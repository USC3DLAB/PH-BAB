//
//  ph.cpp
//  PH-BB
//
//  Created by Semih Atakan on 6/29/15.
//  Copyright (c) 2015 3D Lab. All rights reserved.
//

#include "ph.h"

ph::ph () {}

ph::~ph () {
	if (foutput) ph_output.close();
}

void ph::initialize(SPprob *problem)
{
	prob	= problem;
	
	prob->convert_to_LP();
	
	/* initialize variables */
	inc.resize		(prob->numberOfVar[0]);
	x.resize		(prob->nb_scen);
	lambda.resize	(prob->nb_scen);
	y.resize		(prob->nb_scen);
	
	for (int s=0; s<prob->nb_scen; s++) {
		x[s].resize(prob->numberOfVar[0], 0.0);
		y[s].resize(prob->numberOfVar[1], 0.0);
		
		lambda[s].resize(prob->numberOfVar[0], 0.0);
	}
	
	itr_lim = 1e4;
	tot_nb_itr = 0;
	
	// for display the progress of PH
	verbose = false;
	foutput = false;
	
	if (foutput) ph_output.open("ph_output.log");
}

void ph::load(node *cur_node)
{
	nodeptr = cur_node;
	
	if (nodeptr->nodeid != 0)		{
		// initialize the dual solution
		lambda = nodeptr->dual_soln;
		
		// initialize the primal incumbent solution
		for (int j=0; j<prob->numberOfVar[0]; j++) {
			inc[j] = nodeptr->prim_soln[j];
			
			// projection to the new feasible set
			if (inc[j] > nodeptr->upperbounds[j])	inc[j] = nodeptr->upperbounds[j];
			if (inc[j] < nodeptr->lowerbounds[j])	inc[j] = nodeptr->lowerbounds[j];
		}
	}
	
	if (nodeptr->last_rho == -1)	{ rho = PHBB_rho; coef = 4; }
	else							{ rho = nodeptr->last_rho; coef = nodeptr->last_coef; }
	
	prob->change_t1var_bounds( nodeptr->lowerbounds, nodeptr->upperbounds );
}

void ph::update_inc ()
{
	for (int j=0; j<prob->numberOfVar[0]; j++)	inc[j] = x[0][j] * prob->sceProb[0];
	for (int s=1; s<prob->nb_scen; s++) {
		for (int j=0; j<prob->numberOfVar[0]; j++)
			inc[j] += x[s][j] * prob->sceProb[s];
	}
}

bool ph::isNonanticipative()
{
	for (int j=0; j<prob->numberOfVar[0]; j++) {
		if ( fabs(nodeptr->upperbounds[j] - nodeptr->lowerbounds[j]) > PHBB_Eps )	return false;
	}
	return true;
}

bool ph::optimize_as_mip (double &cutoff)
{
	cout << "Subproblems are solved as MIP... " << flush;
	
#ifdef PHBB_Parallel
	
	// convert to MIP / remove_ph_regularizer / use_linear_objective
	{
		boost::thread_group threads;
		for (int k=0; k<prob->nb_threads; k++)	threads.create_thread( boost::bind(&SPprob::convert_to_MIP, boost::ref(prob), k) );
		threads.join_all();
		for (int k=0; k<prob->nb_threads; k++)	threads.create_thread( boost::bind(&SPprob::remove_ph_regularizer, boost::ref(prob), k ) ) ;
		threads.join_all();
		for (int k=0; k<prob->nb_threads; k++)	threads.create_thread( boost::bind(&SPprob::use_linear_objective, boost::ref(prob), k ) ) ;
		threads.join_all();
	}
	
	vector<double> objs (prob->nb_scen, 0.0), bnodes (prob->nb_scen, 0.0);
	
	// if an incumbent is readily available, solve the root relaxations first
	if (cutoff < PHBB_PosInfty)
	{
		/***** Parallel programming stuff (START) *****/
		boost::asio::io_service io_service;				// create an io_service
		boost::asio::io_service::work work(io_service);	// and some work to stop its run() function from exiting if it has nothing else to do
		boost::thread_group threads;					// start some worker threads
		
		thread_map.clear();
		for (int k=0; k<prob->nb_threads; k++)	{
			boost::thread *new_thread = threads.create_thread(boost::bind(&boost::asio::io_service::run, &io_service));
			thread_map.insert( pair<boost::thread::id, unsigned short> (new_thread->get_id(), k) );
		}
		/***** Parallel programming stuff (END)   *****/
		
		// Solve Subproblems
		for (int s=0; s<prob->nb_scen; s++)
			io_service.post( boost::bind(&ph::optimize_as_mip_root, boost::ref(*this), s, boost::ref(objs[s])) );
		
		/***** Parallel programming stuff (START) *****/
		work.~work();		// let the io_service shutdown by removing the thread
		threads.join_all();
		/***** Parallel programming stuff (END)   *****/
		
		obj_val = 0;
		for (int s=0; s<prob->nb_scen; s++)	obj_val += objs[s];
		
		if (obj_val == PHBB_PosInfty)	{
			// convert to LP
			{
				boost::thread_group threads;
				for (int k=0; k<prob->nb_threads; k++)	threads.create_thread( boost::bind(&SPprob::convert_to_LP, boost::ref(prob), k) );
				threads.join_all();
			}
			cout << endl;
			return false;	// infeasible subproblem
		}
	}
	
	// if the lower bounds already exceed the cutoff value, no need to explore this node any further
	if (cutoff < obj_val)	{
		nodeptr->isOptimal = true; // it is not optimal, yet it will be fathomed in any case
		
		nodeptr->ub				= obj_val;
		nodeptr->lb				= obj_val;
		nodeptr->obj_val		= obj_val;
		
		// convert to LP
		{
			boost::thread_group threads;
			for (int k=0; k<prob->nb_threads; k++)	threads.create_thread( boost::bind(&SPprob::convert_to_LP, boost::ref(prob), k) );
			threads.join_all();
		}
		cout << "cutoff [obj: " << obj_val << "]" << endl;
		
		return true;
	}
	
	// otherwise, continue searching for an integer solution

	/***** Parallel programming stuff (START) *****/
	boost::asio::io_service io_service;				// create an io_service
	boost::asio::io_service::work work(io_service);	// and some work to stop its run() function from exiting if it has nothing else to do
	boost::thread_group threads;					// start some worker threads
	
	thread_map.clear();
	for (int k=0; k<prob->nb_threads; k++)	{
		boost::thread *new_thread = threads.create_thread(boost::bind(&boost::asio::io_service::run, &io_service));
		thread_map.insert( pair<boost::thread::id, unsigned short> (new_thread->get_id(), k) );
	}
	/***** Parallel programming stuff (END)   *****/
	
	// Solve Subproblems
	for (int s=0; s<prob->nb_scen; s++)
		io_service.post( boost::bind(&ph::optimize_as_mip_full, boost::ref(*this), s, boost::ref(bnodes[s]), boost::ref(objs[s])) );
	
	/***** Parallel programming stuff (START) *****/
	work.~work();		// let the io_service shutdown by removing the thread
	threads.join_all();
	/***** Parallel programming stuff (END)   *****/
	
	obj_val = 0;
	double bestbound = 0;
	for (int s=0; s<prob->nb_scen; s++)	{
		obj_val		+= objs[s];
		bestbound	+= bnodes[s];
	}
	
	if (obj_val == PHBB_PosInfty)
	{
		// convert to LP
		{
			boost::thread_group threads;
			for (int k=0; k<prob->nb_threads; k++)	threads.create_thread( boost::bind(&SPprob::convert_to_LP, boost::ref(prob), k) );
			threads.join_all();
		}
		cout << endl;
		return false;
	}
	
	nodeptr->isOptimal = true;
	
	store_solution(nodeptr->isOptimal);
	nodeptr->lb = bestbound;
	nodeptr->ub = obj_val;
	
	// convert to LP
	{
		boost::thread_group threads;
		for (int k=0; k<prob->nb_threads; k++)	threads.create_thread( boost::bind(&SPprob::convert_to_LP, boost::ref(prob), k) );
		threads.join_all();
	}
	
	cout << endl;
	return true;
	
#else
	
	int solver_index = 0;
	
	prob->convert_to_MIP(solver_index);
	prob->remove_ph_regularizer(solver_index);
	prob->use_linear_objective(solver_index);
	
	obj_val = 0;
	
	// if an incumbent is readily available, solve the root relaxations first
	if (cutoff < PHBB_PosInfty)
	{
		obj_val = 0;
		for (int s=0; s<prob->nb_scen; s++) {
			
			prob->setup_scenario_subproblem (s, solver_index);
			
			// Solve this node
			prob->change_t2var_bounds(nodeptr->lowerbounds, nodeptr->upperbounds, s, solver_index);
			
			prob->solver[solver_index].cplex.setParam(IloCplex::NodeLim, 1);
			bool status = prob->solve(solver_index);
			prob->solver[solver_index].cplex.setParam(IloCplex::NodeLim, 2100000000);	// default node limit
			
			if (status)	{
				obj_val += prob->sceProb[s] * prob->solver[solver_index].cplex.getBestObjValue();	// get the best bound
			}
			else {
				cout << prob->solver[solver_index].cplex.getCplexStatus() << endl;
				prob->convert_to_LP(solver_index);
				return false;	// infeasible
			}
		}
	}
	
	// if the lower bounds already exceed the cutoff value, no need to explore this node any further
	if (cutoff < obj_val)	{
		nodeptr->isOptimal = true; // it is not optimal, yet it will be fathomed in any case
		
		nodeptr->ub				= obj_val;
		nodeptr->lb				= obj_val;
		nodeptr->obj_val		= obj_val;
		
		prob->convert_to_LP(solver_index);
		
		cout << "cutoff [obj: " << obj_val << "]" << endl;
		
		return true;
	}
	
	// otherwise, continue finding an integer solution
	
	obj_val=0;
	double bestbound = 0;
	for (int s=0; s<prob->nb_scen; s++)
	{
		prob->setup_scenario_subproblem (s, solver_index);
		
		// Solve of this node
		prob->change_t2var_bounds(nodeptr->lowerbounds, nodeptr->upperbounds, s, solver_index);
		
		bool status = prob->solve(solver_index);
		
		if (status)	{
			prob->getSoln(x[s], 1, solver_index);	// get 1st-stage solution
			prob->getSoln(y[s], 2, solver_index);	// get 2nd-stage solution
			
			obj_val		+= prob->sceProb[s] * prob->solver[solver_index].cplex.getObjValue();
			bestbound	+= prob->sceProb[s] * prob->solver[solver_index].cplex.getBestObjValue();
		}
		else {
			prob->convert_to_LP(solver_index);
			
			return false;	// infeasible
		}
	}
	
	nodeptr->isOptimal = true;
	
	store_solution(nodeptr->isOptimal);
	nodeptr->lb = bestbound;
	nodeptr->ub = obj_val;
	
	prob->convert_to_LP(solver_index);
	
	cout << endl;
	return true;
	
#endif
}

bool ph::isOptimal()
{
	double norm_diff = 0, temp = 0;
	for (int s=0; s<x.size(); s++)
	{
		temp = 0;
		
		for (int j=0; j<x[s].size(); j++)	temp += pow(x[s][j] - inc[j], 2.0);
		norm_diff += sqrt(temp) * prob->sceProb[s];
		
		if (nodeptr->nodeid == 0){
			if(norm_diff > PHBB_NormTol)	return false;
		}
		else {
			if(norm_diff > PHBB_NormTol)	return false;
		}
	}
	return true;
}

void ph::stabilization(double &obj_val, double &prev_obj_val)
{
	double rel_gap = get_rel_gap(obj_val, prev_obj_val);
	if (rel_gap >  0.00) {
		if (direction >= 0)	direction++;
		else	direction = 1;
	}
	else if (rel_gap <= 0.00) {
		if (direction <= 0)	direction--;
		else	direction = -1;
	}
	
	if (direction == 2) {
		rho *= coef;
		//				rho = min(100.0, rho);
		direction = 0;
	}
	else if (direction == -2) {
		rho /= coef;
		//				rho = max(100.0, rho);
		direction = 0;
	}
	
	if (fabs(rel_gap) <= 0.01 && abs(direction) > 1) 	{
		nb_stability++;
		//				cout << direction << endl;
	}
	else						nb_stability = 0;
	
	if (nb_stability == 10)
	{
		coef = (coef + 1)/2.0;
		nb_stability = 0;
	}

}

/******************************************************************************
 Computes the objective function, associated with the minimization subproblems
 i.e., sum_s c_s x_s + d_s y_s + lambda_s x_s
 for a given x, y, and lambda
******************************************************************************/
double ph::compute_dual_function_obj_value (vector< vector<double> > &x,
											vector< vector<double> > &y,
											vector< vector<double> > &lambda,
											int thread)
{
	if (prob->probType.compare("qp") == 0 || prob->probType.compare("miqp") == 0 || prob->probType.compare("miqcp") == 0) {
		cout << "Error: compute_dual_function_obj_value is not suitable for quadratic problems" << endl;
	}
	
	double f = 0;
	for (int s=0; s<prob->nb_scen; s++) {
		// first-stage obj coefs
		for (int i=0; i<prob->numberOfVar[0]; i++)	f += prob->sceProb[s] * (prob->obj_coefs[thread][0][i] + lambda[s][i]) * x[s][i];
		
		// second-stage obj coefs
		for (int i=0; i<prob->numberOfVar[1]; i++)	f += prob->sceProb[s] * prob->obj_coefs[thread][1][i] * y[s][i];
	}
	return f;
}

//double ph::compute_subprob_primal_obj_value (vector< IloNumArray > &x,
//											 vector< IloNumArray > &y)
//{
//	if (prob->probType.compare("qp") == 0 || prob->probType.compare("miqp") == 0 || prob->probType.compare("miqcp") == 0) {
//		cout << "Error: compute_dual_function_obj_value is not suitable for quadratic problems" << endl;
//	}
//	
//	double f = 0;
//	for (int s=0; s<prob->nb_scen; s++) {
//		// first-stage obj coefs
//		for (int i=0; i<prob->numberOfVar[0]; i++)	f += prob->sceProb[s] * prob->obj_coefs[0][i] * x[s][i];
//		
//		// second-stage obj coefs
//		for (int i=0; i<prob->numberOfVar[1]; i++)	f += prob->sceProb[s] * prob->obj_coefs[1][i] * y[s][i];
//	}
//	return f;
//}

/******************************************************************************
 Computes the objective function, associated with the minimization subproblems
 i.e., sum_s c_s x_s + d_s y_s for a given nonanticipative x, and arbitrary y
 ******************************************************************************/
//double ph::compute_primal_obj_value (IloNumArray &x, vector< IloNumArray > &y)
//{
//	if (prob->probType.compare("qp") == 0 || prob->probType.compare("miqp") == 0 || prob->probType.compare("miqcp") == 0) {
//		cout << "Error: compute_primal_obj_value is not suitable for quadratic problems" << endl;
//	}
//
//	double f = 0;
//	for (int s=0; s<prob->nb_scen; s++) {
//		// first-stage obj coefs
//		for (int i=0; i<prob->numberOfVar[0]; i++)	f += prob->sceProb[s] * prob->obj_coefs[0][i] * x[i];
//		
//		// second-stage obj coefs
//		for (int i=0; i<prob->numberOfVar[1]; i++)	f += prob->sceProb[s] * prob->obj_coefs[1][i] * y[s][i];
//	}
//	return f;
//}


/******************************************************************************
 Computes the Lagrangian lower bound, and upper bound (fixed soln)
 ******************************************************************************/
void ph::update_lb_ub ()
{
	
#ifdef PHBB_Parallel
	
	vector<double> lbs (prob->nb_scen, 0.0), ubs (prob->nb_scen, 0.0);
	
	// remove the ph regularizer
	{
		boost::thread_group threads;
		for (int k=0; k<prob->nb_threads; k++)	threads.create_thread( boost::bind(&SPprob::remove_ph_regularizer, boost::ref(prob), k ) ) ;
		threads.join_all();
	}
	
	/***** Parallel programming stuff (START) *****/
	boost::asio::io_service io_service;				// create an io_service
	boost::asio::io_service::work work(io_service);	// and some work to stop its run() function from exiting if it has nothing else to do
	boost::thread_group threads;					// start some worker threads
	
	thread_map.clear();
	for (int k=0; k<prob->nb_threads; k++)	{
		boost::thread *new_thread = threads.create_thread(boost::bind(&boost::asio::io_service::run, &io_service));
		thread_map.insert( pair<boost::thread::id, unsigned short> (new_thread->get_id(), k) );
	}
	/***** Parallel programming stuff (END)   *****/
	
	// Solve Subproblems
	for (int s=0; s<prob->nb_scen; s++)
		io_service.post( boost::bind(&ph::update_lb_ub_scenario_subproblem, boost::ref(*this), s, boost::ref(lbs[s]), boost::ref(ubs[s])) );
	
	/***** Parallel programming stuff (START) *****/
	work.~work();		// let the io_service shutdown by removing the thread
	threads.join_all();
	/***** Parallel programming stuff (END)   *****/

	double lb_candidate = 0.0, ub_candidate = 0.0;
	for (int s=0; s<prob->nb_scen; s++)	{
		lb_candidate += lbs[s];
		ub_candidate += ubs[s];
	}
	if (ub_candidate < nodeptr->ub)	nodeptr->ub = ub_candidate;
	if (lb_candidate > nodeptr->lb)	nodeptr->lb = lb_candidate;

#else 
	
	int solver_index = 0;
	
	bool status = false;
	
	// remove the ph regularizer
	prob->remove_ph_regularizer(solver_index);

	double lb_candidate = 0.0, ub_candidate = 0.0;
	
	for (int s=0; s<prob->nb_scen; s++) {
		
		// change 2nd-stage variable bounds (first-stage is already changed)
		prob->change_t2var_bounds(nodeptr->lowerbounds, nodeptr->upperbounds, s, solver_index);
		
		// setup the scenario subproblem
		prob->setup_scenario_subproblem(s, solver_index);
		
		// begin with the saddle objective
		prob->use_saddle_objective(lambda[s], solver_index);

		// solve the lower-bounding problem, get the bound
		status = prob->solve(solver_index);
		if (!status) {
			cout << "Could not minimize the Lagrangian function for this reason: " << prob->solver[solver_index].cplex.getCplexStatus() << endl;
		}
		else {
			lb_candidate += prob->sceProb[s] * prob->getObjValue(solver_index);
		}
		
		if (false)	// UB CHECK IS DISABLED
		{
			// remove the lagrangian terms
			prob->use_linear_objective(solver_index);
			
			// fix the first-stage solution
			prob->change_t1var_bounds(inc, inc, solver_index);
			
			// solve the upper-bounding problem, get the bound
			bool status = prob->solve(solver_index);
			if (!status) {
				if (prob->solver[solver_index].cplex.getCplexStatus() == IloCplex::Infeasible) ub_candidate = PHBB_PosInfty;
				else {
					cout << "Could not minimize the primal function for this reason: " << prob->solver[solver_index].cplex.getCplexStatus() << endl;
				}
			}
			else {
				ub_candidate += prob->sceProb[s] * prob->getObjValue(solver_index);
			}
			
			// reset the 1st-stage variable bounds
			prob->change_t1var_bounds(nodeptr->lowerbounds, nodeptr->upperbounds, solver_index);
		}
		else {
			ub_candidate = PHBB_PosInfty;
		}
	}

	if (ub_candidate < nodeptr->ub)	nodeptr->ub = ub_candidate;
	if (lb_candidate > nodeptr->lb)	nodeptr->lb = lb_candidate;
	
#endif
}

bool ph::optimize_root()
{

#ifdef PHBB_Parallel
	
	vector<double> objs (prob->nb_scen, 0.0);
	
	/***** Parallel programming stuff (START) *****/
	boost::asio::io_service io_service;				// create an io_service
	boost::asio::io_service::work work(io_service);	// and some work to stop its run() function from exiting if it has nothing else to do
	boost::thread_group threads;					// start some worker threads
	
	thread_map.clear();
	for (int k=0; k<prob->nb_threads; k++)	{
		boost::thread *new_thread = threads.create_thread(boost::bind(&boost::asio::io_service::run, &io_service));
		thread_map.insert( pair<boost::thread::id, unsigned short> (new_thread->get_id(), k) );
	}
	/***** Parallel programming stuff (END)   *****/
	
	// solve subproblems
	for (int s=0; s<prob->nb_scen; s++)	io_service.post( boost::bind(&ph::optimize_root_scenario_subproblem, boost::ref(*this), s, boost::ref(objs[s])) );
	
	/***** Parallel programming stuff (START) *****/
	work.~work();		// let the io_service shutdown by removing the thread
	threads.join_all();
	/***** Parallel programming stuff (END)   *****/
	
	obj_val = 0;
	for (int s=0; s<prob->nb_scen; s++)	obj_val += objs[s];
	
	if (obj_val == PHBB_PosInfty)	return false;
#else
	
	int solver_index = 0;
	
	/** INITIALIZATION **/
	bool cont = true;
	
	/*** ITERATION 0: Root Processing ***/
	// solve scenario subproblems (w/o regularization)

	obj_val = 0;
	for (int s=0; s<prob->nb_scen; s++) {
		prob->setup_scenario_subproblem(s, solver_index);
		
		cont = prob->solve(solver_index);
		
		if (cont) {
			prob->getSoln(x[s], 1, solver_index);	// get 1st stage soln
			prob->getSoln(y[s], 2, solver_index);	// get 2nd stage soln
			
			obj_val += prob->sceProb[s] * prob->getObjValue(solver_index);
		}
		else {
			cout << prob->solver[solver_index].cplex.getStatus() << " " << prob->solver[solver_index].cplex.getCplexStatus() << endl;
			return false;
		}
	}
	
#endif 
	
	update_inc();

	update_dual_soln();
	
	nodeptr->lb = obj_val;
	
	if (verbose)	{
		cout << setw(5) << 0 << fixed << setprecision(6) << " " << obj_val << " " << nodeptr->lb << endl;
	}
	if (foutput)	{
		ph_output << setw(5) << 0 << fixed << setprecision(6) << " " << obj_val << " " << nodeptr->lb << endl;
	}

	bool isOpt = false;
	
	store_solution(isOpt);
	
	tot_nb_itr++;
	
	return true;
}

void ph::update_dual_soln ()
{
	for (int s=0; s<prob->nb_scen; s++) {
		for (int j=0; j<lambda[0].size(); j++) {
			lambda[s][j] += ( rho * ( x[s][j] - inc[j] ) );
		}
	}
}

bool ph::update_primal_soln()
{
	
#ifdef PHBB_Parallel
	
	// add_ph_regularizer
	{
		boost::thread_group threads;
		for (int k=0; k<prob->nb_threads; k++)	threads.create_thread( boost::bind(&SPprob::add_ph_regularizer, boost::ref(prob), boost::ref(rho), k ) ) ;
		threads.join_all();
	}

	vector<double> objs (prob->nb_scen, 0.0);
	
	{
		/***** Parallel programming stuff (START) *****/
		boost::asio::io_service io_service;				// create an io_service
		boost::asio::io_service::work work(io_service);	// and some work to stop its run() function from exiting if it has nothing else to do
		boost::thread_group threads;					// start some worker threads
		
		thread_map.clear();
		for (int k=0; k<prob->nb_threads; k++)	{
			boost::thread *new_thread = threads.create_thread(boost::bind(&boost::asio::io_service::run, &io_service));
			thread_map.insert( pair<boost::thread::id, unsigned short> (new_thread->get_id(), k) );
		}
		/***** Parallel programming stuff (END)   *****/
		
		
		// solve subproblems
		for (int s=0; s<prob->nb_scen; s++)	io_service.post( boost::bind(&ph::update_primal_soln_scenario_subproblem, boost::ref(*this), s, boost::ref(objs[s])) );
		
		/***** Parallel programming stuff (START) *****/
		work.~work();		// let the io_service shutdown by removing the thread
		threads.join_all();
		
		
		/***** Parallel programming stuff (END)   *****/
	}
	obj_val = 0;
	for (int s=0; s<prob->nb_scen; s++)	obj_val += objs[s];
	
	if (obj_val == PHBB_PosInfty)	return false;
	else							return true;
	
#else
	
	int solver_index = 0;
	double constant = 0;
	
	prob->add_ph_regularizer(rho, solver_index);
	
	obj_val = 0;
	for (int s=0; s<prob->nb_scen; s++) {
		
		prob->change_t2var_bounds(nodeptr->lowerbounds, nodeptr->upperbounds, s, solver_index);
		
		// prepare the scenario subproblems
		constant = prob->setup_regularized_subproblem (s, inc, lambda[s], rho, solver_index);
		
		//	add_scenario_cuts(s);
		
		double tmp = get_wall_time();
		bool status = prob->solve(solver_index);
		ph_t += get_wall_time() - tmp;
		
		// if the Barrier optimizer fails, do a reoptimization with numerical emphasis
		if (prob->solver[solver_index].cplex.getCplexStatus() == IloCplex::NumBest) {
			prob->solver[solver_index].cplex.setParam(IloCplex::NumericalEmphasis, 1);
			status = prob->solve(solver_index);
			cout << "Reoptimizing subproblem with numerical emphasis. Result: " << prob->solver[solver_index].cplex.getCplexStatus() << endl;
			prob->solver[solver_index].cplex.setParam(IloCplex::NumericalEmphasis, 0);
		}
		
		if (status)
		{
			prob->getSoln(x[s], 1, solver_index);	// get 1st stage solution
			prob->getSoln(y[s], 2, solver_index);	// get 2nd stage solution
			
			obj_val += prob->sceProb[s] * (prob->getObjValue(solver_index) + constant);
			//			remove_scenario_cuts(s);	// remove scenario-specific cuts from the model
		}
		else {	// infeasible scenario problem
			//			remove_scenario_cuts(s);	// remove scenario-specific cuts from the model
			return false;
		}		
	}
	return true;

#endif
}

void ph::store_solution(bool &isOpt)
{
	// store primal solutions
	if (nodeptr->prim_soln.size() == 0)	nodeptr->prim_soln.resize( prob->numberOfVar[0] + prob->numberOfVar[1] * prob->nb_scen );
	
	int k=0;
	for (; k<inc.size(); k++) nodeptr->prim_soln[k] = inc[k];
	
	for (int s=0; s<prob->nb_scen; s++) {
		for (int j=0; j<y[0].size(); j++)	{
			nodeptr->prim_soln[k] = y[s][j];
			k++;
		}
	}
	
	// store dual solutions
	nodeptr->dual_soln	= lambda;
	
	// store the PH obj value
	nodeptr->obj_val	= obj_val;
	
	// lower and upper bounds are already updated
	// ..
	
	// if optimal, raise the flag. Otherwise, store the parameters for later use.
	if (isOpt) 	{
		nodeptr->isOptimal = true;
	}
	else		{
		nodeptr->nb_of_passes++;
		nodeptr->last_rho	= rho;
		nodeptr->last_coef	= coef;
	}
}

/******************************************************************************
 Check if the upper and lower cutoff limits, and the iteration limit are 
 exceeded. If they are, stop. Otherwise, continue.
 ******************************************************************************/
bool ph::check_termination_criteria ()
{
	if ( nodeptr->lb > cutoff )		return false;
	if ( itr > itr_lim )			return false;
	if ( haveBoundsConverged() )	return false;

	return true;
}

bool ph::haveBoundsConverged() {
	return ( get_rel_gap(nodeptr->ub, nodeptr->lb) < PHBB_NodeRelOptTol );
}

bool ph::optimize_w_bounds(int itr_lim, double cutoff)
{
	// set the limits
	this->itr_lim = itr_lim;
	this->cutoff = cutoff;

	itr = 1;
	
	obj_val = nodeptr->obj_val;
	
	bool cont = true;
	
	while (cont) {
		
		prev_obj_val = obj_val;
		
		/* Primal iteration (wt regularization term) */
		cont = update_primal_soln ();
		if (!cont)	{ return false; }
		
		update_inc();
		
		/* Dual iteration (wt implicit regularization) */
		update_dual_soln();
		
		// get the bounds on the objective value
		if (itr == itr_lim || itr % PHBB_NodeBndItrLim == 0)	update_lb_ub();
		
		if (verbose)	{
			cout << setw(5) << itr << fixed << setprecision(6) << " " << obj_val << " " << nodeptr->lb << " " << nodeptr->ub << " " << setprecision(2) << fabs(get_rel_gap(nodeptr->lb, nodeptr->ub)) * 100 << endl;
		}
		if (foutput)	{
			ph_output << setw(5) << itr << fixed << setprecision(6) << " " << obj_val << " " << nodeptr->lb << " " << nodeptr->ub << " " << setprecision(2) << fabs(get_rel_gap(nodeptr->lb, nodeptr->ub)) * 100 << endl;
		}
		itr++;
		
		cont  = check_termination_criteria();
	
	}
	
	tot_nb_itr += --itr;
	// END OF MAIN LOOP
	
	bool isOpt = haveBoundsConverged();
	store_solution( isOpt );
	
	return true;
}

bool ph::optimize()
{
	itr = 1;
	
	obj_val = nodeptr->obj_val;
	
	direction=0, nb_stability=0;
	
	bool cont = true, isOpt = false;
	
	while (!isOpt) {
		
		prev_obj_val = obj_val;
		
		/* Primal iteration (wt regularization term) */
		cont = update_primal_soln ();
		if (!cont)	{ return false; }
		
		update_inc();
		
		/* Dual iteration (wt implicit regularization) */
		update_dual_soln();
		
		// rho stabilization
//		stabilization(obj_val, prev_obj_val);
		
		if (verbose)	{
			cout << setw(5) << itr << fixed << setprecision(6) << " " << obj_val << " " << rho << " " << coef << endl;
		}
		if (foutput)	{
			ph_output << setw(5) << itr << fixed << setprecision(6) << " " << obj_val << " " << rho << " " << coef << endl;
		}

		itr++;
		
		isOpt = isOptimal();
	}
	
	tot_nb_itr += --itr;
	// END OF MAIN LOOP
	
	nodeptr->lb = obj_val;
	nodeptr->ub = obj_val;
	store_solution(isOpt);
	
	return true;
}

bool ph::heuristic()
{
	ph_t = clock();
	
	// initialize the nodeptr with a dummy node
	nodeptr = new node(prob->lb, prob->ub);
	nodeptr->lb = PHBB_NegInfty;
	nodeptr->ub = PHBB_PosInfty;
	
	// convert the problem to MIP
	prob->convert_to_MIP();
	
	// create a list of obj values to check if the algorithm is repeating itself
	set<double> obj_vals;
	unsigned int count_of_repeated_objs = 0;
	
	itr = 1;
	obj_val = 0.0;
	rho = PHBB_rho;

	cout << "\nRunning progressive hedging with MIP subproblems" << endl;
	
	printf("  Itr        ObjVal   Time\n");
	
	double itr_t = clock();
	optimize_root();
	printf("%5d %13.4f %6.1f\n", itr, obj_val, (clock()-itr_t)/CLOCKS_PER_SEC);
	
	bool isOpt = isOptimal(), cont = true;
	
	itr++;
	
	while (!isOpt) {
		
		itr_t = clock();
		
		/* Primal iteration (wt regularization term) */
		cont = update_primal_soln ();
		if (!cont)	{ return false; }
		
		update_inc();
		
		/* Dual iteration (wt implicit regularization) */
		update_dual_soln();
	
		printf("%5d %13.4f %6.1f\n", itr, obj_val, (clock()-itr_t)/CLOCKS_PER_SEC);

		itr++;
		
		isOpt = isOptimal();
		
		// check if the algorithm starts repeating itself
		if ( (obj_vals.find(round(obj_val * 1e6)/1e6)) == obj_vals.end() ) {
			obj_vals.insert( round(obj_val * 1e6)/1e6 ) ;
			count_of_repeated_objs = 0;
		} else { count_of_repeated_objs++; }
		
		if (count_of_repeated_objs >= 10) {
			cout << "Algorithm is not converging" << endl;
			break;
		}
		
		// check time limit
		if ( (clock()-ph_t)/CLOCKS_PER_SEC > PHBB_TiLim	) {
			cout << "Time limit has been reached" << endl;
			break;
		}
	}
	
	tot_nb_itr = itr-1;
	// END OF MAIN LOOP

	ph_t = clock() - ph_t;

	// round the incumbent to get a feasible solution
	for (int j=0; j<prob->numberOfVar[0]; j++) {
		if (prob->integrality[j])	inc[j] = round(inc[j]);
	}
	
	// get lower and upper bounds (integrality restrictions are still imposed)
	update_lb_ub();
	
	printOutput(nodeptr->lb, nodeptr->ub);
	
	cout << endl << fixed << setprecision(6) << "PH Obj Val: " << obj_val << " LB: " << nodeptr->lb << " UB: " << nodeptr->ub << endl;
	
	// remove the dummy node
	delete nodeptr;
	nodeptr = NULL;
	
	return true;
}

bool ph::optimize_w_bounds_mip_option (int itr_lim, double cutoff)
{
	if (isNonanticipative()) {
		return optimize_as_mip(cutoff);
	}
	else {
		return optimize_w_bounds(itr_lim, cutoff);
	}
}

bool ph::optimize_w_mip_option (double cutoff)
{
	if (isNonanticipative())	{
		return optimize_as_mip(cutoff);
	}
	else  {
		return optimize();
	}
}

double ph::get_rel_gap(const double &compare_this, double &wrt_this) {
	return ((compare_this - wrt_this)/(fabs(wrt_this) + PHBB_Eps));
}

void ph::printOutput (double &lb, double &ub)
{
	ofstream output;
	output.open("output.log", ios::app);
	
	output << currentDateTime() << "\t" << prob->probName << "\t" << obj_val << "\t" << lb << "\t" << ub << "\t" << tot_nb_itr << "\t" << double(ph_t)/CLOCKS_PER_SEC << endl;
}

#ifdef PHBB_Parallel

void ph::update_primal_soln_scenario_subproblem (int s, double &obj)
{
	int solver_index = thread_map[boost::this_thread::get_id()];
	
	prob->change_t2var_bounds(nodeptr->lowerbounds, nodeptr->upperbounds, s, solver_index);
	
	// prepare the scenario subproblems
	double constant = prob->setup_regularized_subproblem (s, inc, lambda[s], rho, solver_index);
	
	double tmp = get_wall_time();
	bool status = prob->solve(solver_index);
	ph_t += get_wall_time() - tmp;

	// if the Barrier optimizer fails, do a reoptimization with numerical emphasis
	if (prob->solver[solver_index].cplex.getCplexStatus() == IloCplex::NumBest) {
		prob->solver[solver_index].cplex.setParam(IloCplex::NumericalEmphasis, 1);
		status = prob->solve(solver_index);
		cout << "Reoptimizing subproblem with numerical emphasis. Result: " << prob->solver[solver_index].cplex.getCplexStatus() << endl;
		prob->solver[solver_index].cplex.setParam(IloCplex::NumericalEmphasis, 0);
	}
	
	if (status)
	{
		prob->getSoln(x[s], 1, solver_index);	// get 1st stage solution
		prob->getSoln(y[s], 2, solver_index);	// get 2nd stage solution
		
		obj = prob->sceProb[s] * (prob->getObjValue(solver_index) + constant);
	}
	else {	// infeasible scenario problem
		obj = PHBB_PosInfty;
	}
}

void ph::update_lb_ub_scenario_subproblem (int s, double &lb_obj, double &ub_obj)
{
	int solver_index = thread_map[boost::this_thread::get_id()];
	
	// change 2nd-stage variable bounds (first-stage is already changed)
	prob->change_t2var_bounds(nodeptr->lowerbounds, nodeptr->upperbounds, s, solver_index);
	
	// setup the scenario subproblem
	prob->setup_scenario_subproblem(s, solver_index);
	
	// begin with the saddle objective
	prob->use_saddle_objective(lambda[s], solver_index);
	
	// solve the lower-bounding problem, get the bound
	bool status = prob->solve(solver_index);
	if (!status) {
		cout << "Could not minimize the Lagrangian function for this reason: " << prob->solver[solver_index].cplex.getCplexStatus() << endl;
		lb_obj = PHBB_NegInfty;
	}
	else {
		lb_obj = prob->sceProb[s] * prob->getObjValue(solver_index);
	}
	
	if (false)	// UB CHECK IS DISABLED TEMPORARILY
	{
		// remove the lagrangian terms
		prob->use_linear_objective(solver_index);
		
		// fix the first-stage solution
		prob->change_t1var_bounds(inc, inc, solver_index);
		
		// solve the upper-bounding problem, get the bound
		bool status = prob->solve(solver_index);
		if (!status) {
			if (prob->solver[solver_index].cplex.getCplexStatus() == IloCplex::Infeasible) {
				ub_obj = PHBB_PosInfty;
			}
			else {
				cout << "Could not minimize the primal function for this reason: " << prob->solver[solver_index].cplex.getCplexStatus() << endl;
				ub_obj = PHBB_PosInfty;
			}
		}
		else {
			ub_obj = prob->sceProb[s] * prob->getObjValue(solver_index);
		}
		
		// reset the 1st-stage variable bounds
		prob->change_t1var_bounds(nodeptr->lowerbounds, nodeptr->upperbounds, solver_index);
	}
	else {
		ub_obj = PHBB_PosInfty;
	}
}

/******************************************************************************
 Optimizes the subproblem as a mixed-integer program, but stops after the
 root relaxation is fully processed. In other words, no branching is allowed.
 ******************************************************************************/
void ph::optimize_as_mip_root (int s, double &obj)
{
	int solver_index = thread_map[boost::this_thread::get_id()];
	
	prob->setup_scenario_subproblem (s, solver_index);
	
	// Solve this node
	prob->change_t2var_bounds(nodeptr->lowerbounds, nodeptr->upperbounds, s, solver_index);
	
	prob->solver[solver_index].cplex.setParam(IloCplex::NodeLim, 1);
	bool status = prob->solve(solver_index);
	prob->solver[solver_index].cplex.setParam(IloCplex::NodeLim, 2100000000);	// default node limit
	
	if (status)	{
		obj = prob->sceProb[s] * prob->solver[solver_index].cplex.getBestObjValue();	// get the best bound
	}
	else {
		obj = PHBB_PosInfty;
	}
}

/******************************************************************************
 Optimizes the subproblem as a mixed-integer program, to optimality.
 (Compare it to optimize_as_mip_root)
 ******************************************************************************/
void ph::optimize_as_mip_full (int s, double &bnode, double &obj)
{
	int solver_index = thread_map[boost::this_thread::get_id()];

	prob->setup_scenario_subproblem (s, solver_index);
	
	// Solve of this node
	prob->change_t2var_bounds(nodeptr->lowerbounds, nodeptr->upperbounds, s, solver_index);
	
	bool status = prob->solve(solver_index);
	
	if (status)	{
		prob->getSoln(x[s], 1, solver_index);	// get 1st-stage solution
		prob->getSoln(y[s], 2, solver_index);	// get 2nd-stage solution
		
		obj		= prob->sceProb[s] * prob->solver[solver_index].cplex.getObjValue();
		bnode	= prob->sceProb[s] * prob->solver[solver_index].cplex.getBestObjValue();
	}
	else {
		obj		= PHBB_PosInfty;
		bnode	= PHBB_PosInfty;
	}
}

void ph::optimize_root_scenario_subproblem (int s, double &obj)
{
	int solver_index = thread_map[boost::this_thread::get_id()];
	
	prob->setup_scenario_subproblem(s, solver_index);
	
	bool cont = prob->solve(solver_index);
	
	if (cont) {
		prob->getSoln(x[s], 1, solver_index);	// get 1st stage soln
		prob->getSoln(y[s], 2, solver_index);	// get 2nd stage soln
		
		obj = prob->sceProb[s] * prob->getObjValue(solver_index);
	}
	else {
		cout << prob->solver[solver_index].cplex.getStatus() << " " << prob->solver[solver_index].cplex.getCplexStatus() << endl;
		obj = PHBB_PosInfty;
	}
}

#endif
