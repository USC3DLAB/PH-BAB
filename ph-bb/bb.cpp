//
//  bb.cpp
//  PH-BB
//
//  Created by Semih Atakan on 6/7/15.
//  Copyright (c) 2015 3D Lab. All rights reserved.
//

#include "bb.h"

bb::bb () {
	// welcome msg
	cout << "Progressive-Hedging embedded Branch-and-Bound version 3.0" << endl << endl;
}

void bb::load(SPprob &prob)
{
	// reset the solver
	reset ();
	
	// initialize problem pointer of the solver
	problem = &prob;
	
	// create the root node
	node *root = new node(problem->lb, problem->ub);
	root->lb = PHBB_NegInfty;
	root->ub = PHBB_PosInfty;
	
	add_node(root);
	
	// initialize incumbent/lowerbound obj values
	inc_obj_val = PHBB_PosInfty;
	lb_obj_val  = PHBB_NegInfty;
	
	// initialize the solver with the problem parameters
	solver.initialize(problem);
}

bool bb::optimize() { return true; }

void bb::add_node(node *cur_node)
{
	// add the node to the set
	nodes.insert( pair<unsigned int, node*> (cur_node->nodeid, cur_node) );
	
	// add its guaranteed lower bound to lb set
	lb_list.insert( pair<double, unsigned int> (cur_node->lb, cur_node->nodeid) );
	
	// add its priority value to the priority list
	priority_list.push( double_uint_pair(cur_node->lb, cur_node->nodeid) );
}

node* bb::get_next_node()
{
	std::tr1::unordered_map<unsigned int, node*>::iterator it;
	
	// find the top-priority element
	while (true) {
		unsigned int lookup_id = (priority_list.top()).nodeid;	// lookup from priority list
		priority_list.pop();
		
		it = nodes.find( lookup_id );		// look up from the node list
		if (it != nodes.end()) break;	// node not already fathomed
	}
	
	node *nodeptr = (it->second);		// get a handle on the node
	
	nodes.erase(it);					// remove the node from the list
	
	// remove its lower bound from the LB list (it will hopefully change)
	multimap<double, unsigned int>::iterator map_it = lb_list.find( nodeptr->lb );
	while ( map_it->second != nodeptr->nodeid ) {
		map_it++;	// the lb must match the node (i.e. avoid multiple nodes with the same lower bound)
	}
	lb_list.erase(map_it);
	
	return nodeptr;
}

double bb::get_rel_abs_gap(const double &compare_this, double &wrt_this) {
	return (fabs(compare_this - wrt_this)/(fabs(wrt_this) + PHBB_Eps));
}

double bb::get_rel_gap(const double &compare_this, double &wrt_this) {
	return (compare_this > PHBB_NegInfty && wrt_this < PHBB_PosInfty) ? ((compare_this - wrt_this)/(fabs(wrt_this) + PHBB_Eps)) : PHBB_NegInfty;
}

void bb::process_int_node (node *cur_node, bool delete_cur_node)
{
	if ( get_rel_gap(cur_node->ub, inc_obj_val) < PHBB_SolnPoolRelOptTol )
	{
		// add to the solution to the pool
		soln_pool.push( solution(cur_node->prim_soln, cur_node->ub) );
		
		// update the incumbent, if we have a better one
		if ( inc_obj_val > cur_node->ub + PHBB_Eps )
		{
			/*********** incumbent update! *************/
			inc_obj_val = cur_node->ub;
			inc_soln	= cur_node->prim_soln;
			/*******************************************/
			
			
			/*************** fathoming! ******************
			 Nodes that have higher lower bound than
			 the incumbent objective are fathomed.
			 ********************************************/
			if ( lb_list.size() > 0 && (lb_list.rbegin())->first > inc_obj_val )
			{
				multimap<double, unsigned int>::iterator it = lb_list.upper_bound(inc_obj_val);
				
				// remove the nodes
				for (multimap<double, unsigned int>::iterator it2 = it; it2 != lb_list.end(); it2++) {
					delete nodes[it2->second];	// kill the node
					nodes[it2->second] = NULL;
					nodes.erase(it2->second);
				}
				
				// remove the lower bounds
				lb_list.erase(it, lb_list.end());
			}
			/******************* end ********************/
			
			
			/********** solution pool clean up! **********
			 Nodes in the solution pool that does not
			 have an objective within SolnPoolRelGap are
			 removed.
			 ********************************************/
			while ( soln_pool.size() > 1
				   && get_rel_gap((soln_pool.top()).obj_val, inc_obj_val) > PHBB_SolnPoolRelOptTol )
			{
				soln_pool.pop();
			}
			/******************* end ********************/
		}
		if (delete_cur_node) {
			delete cur_node;
			cur_node = NULL;
		}
	}
	else {
		fathomed = true;
		
		no_impr_count++;
		delete cur_node;
		cur_node = NULL;
	}
}

void bb::process_int_node (node *cur_node)
{
	process_int_node(cur_node, true);
}

void bb::process_frac_node(node *cur_node)
{
	int frac_index = maxInfeasVar(cur_node->prim_soln);
//	int frac_index = minInfeasVar(cur_node->prim_soln);
//	cout << frac_index << endl;
//	int frac_index = strong_branch(cur_node);
//	cout << frac_index << endl;
	
	// if the node seems promising to enter into the solution pool, otherwise fathom it
	if ( cur_node->lb < inc_obj_val ) {
		split(cur_node, frac_index);
	}
	else {
		fathomed = true;
		delete cur_node;
		cur_node = NULL;
	}
}

void bb::fix_and_branch(node *cur_node)
{
	if ( cur_node->lb < inc_obj_val )
	{
		// split the node
		
		int not_fixed_index = isFixed(cur_node->lowerbounds, cur_node->upperbounds);
		
		if (not_fixed_index != -1)	// start fixing
		{
			// create the node where the current 1st-stage solution values are separated
			node *node_sep = new node (cur_node->lowerbounds, cur_node->upperbounds);
			
			// create the node where the current 1st-stage solution values are fixed
			// fixing:
			vector<int> fixed_var_indices;
			for (int i=0; i<problem->numberOfVar[0]; i++)
			{
				if (problem->integrality[i] && fabs( cur_node->lowerbounds[i] - cur_node->upperbounds[i] ) > PHBB_Eps)
				{
					if( fabs(cur_node->lowerbounds[i] - cur_node->prim_soln[i]) <= PHBB_IntTol )	{
						cur_node->upperbounds[i] = 0;
						fixed_var_indices.push_back(i);
					}
					else if ( fabs(cur_node->upperbounds[i] - cur_node->prim_soln[i]) <= PHBB_IntTol )	{
						cur_node->lowerbounds[i] = 1;
						fixed_var_indices.push_back(i);
					}
				}
			}
			node *node_fix = new node (cur_node->lowerbounds, cur_node->upperbounds);
			
			// first estimates of the primal objective values are passed from the mother node
			node_sep->lb		= cur_node->lb;
			node_sep->ub		= PHBB_PosInfty;
			node_sep->obj_val	= cur_node->obj_val;
			
			node_fix->lb		= cur_node->lb;
			node_fix->ub		= cur_node->ub;		// the fixed guy also gets the upper bound from the mother
			node_fix->obj_val	= cur_node->obj_val;
			
			node_sep->prim_soln = cur_node->prim_soln;	// the primal solution will be projected in the LP solver for warmer start
			node_fix->prim_soln = cur_node->prim_soln;
			
			node_sep->dual_soln = cur_node->dual_soln;	// the dual solution is passed for a warmer start
			node_fix->dual_soln = cur_node->dual_soln;
			
			if (cur_node->cut_indices.size() > 0)	{	// pass the previous cut indices of the nodes
				node_sep->cut_indices = cur_node->cut_indices;
				node_fix->cut_indices = cur_node->cut_indices;
			}
			
			// prepare the L2-feasibility cut
			
			long new_cut_index = problem->derive_L2_feasibility_cut(cur_node->lowerbounds, fixed_var_indices);	// lowerbounds = ub, is the unique first-stage soln
			
			node_sep->cut_indices.push_back( new_cut_index );
			
			delete cur_node;	// delete the mother
			cur_node = NULL;
			
			add_node(node_fix);	// add the children
			add_node(node_sep);
		}
		else 						// continue branching on the 2nd-stage variables
		{
			cout << "I am not supposed to be here in the current implementation" << endl;
			process_frac_node(cur_node);
		}
	}
	else {
		fathomed = true;
		delete cur_node;
		cur_node = NULL;
	}
}

void bb::split(node *cur_node, int &frac_index)
{
	// node down is the child where x <= floor(\barx), node up is x >= ceil(\barx)
	node *node_d = new node (cur_node->lowerbounds, cur_node->upperbounds);
	node *node_u = new node (cur_node->lowerbounds, cur_node->upperbounds);

//	cout << "X" << frac_index+1 << " is either less than " << floor(cur_node->prim_soln[frac_index] + PHBB_Eps) << " or greater than " << ceil(cur_node->prim_soln[frac_index] - PHBB_Eps) << endl;
//	cout << "     ";

	node_d->upperbounds[frac_index] = floor(cur_node->prim_soln[frac_index] + PHBB_Eps);
	node_u->lowerbounds[frac_index] = ceil(cur_node->prim_soln[frac_index] - PHBB_Eps);

	// first estimates of the primal objective values are passed from the mother node
	node_d->lb = cur_node->lb;
	node_u->lb = cur_node->lb;
	
	node_d->ub = PHBB_PosInfty;
	node_u->ub = PHBB_PosInfty;
	
	node_d->obj_val = cur_node->obj_val;
	node_u->obj_val = cur_node->obj_val;
	
	node_d->prim_soln = cur_node->prim_soln;	// the primal solution will be projected in the LP solver for warmer start
	node_u->prim_soln = cur_node->prim_soln;

	node_d->dual_soln = cur_node->dual_soln;	// the dual solution is passed for a warmer start
	node_u->dual_soln = cur_node->dual_soln;

	node_d->cut_indices = cur_node->cut_indices;	// pass the node specific cuts
	node_u->cut_indices = cur_node->cut_indices;
	
//	cout << cur_node->getId() << " " << node_d->getId() << " " << node_u->getId() << endl;
//	cout << cur_node->lb_obj_val << endl;
	
	add_node(node_d);	// children are born
	add_node(node_u);
	
	delete cur_node;	// delete the mother
	cur_node = NULL;
}

int bb::isFixed(vector<double> &lb, vector<double> &ub)
{
	/*	Checks if the vector have all upper and lower bounds equal.
		Returns the variable index that violates this first.
		Returns -1 if the solution is fixed.
	 */
	
	for (int i=0; i<problem->numberOfVar[0]; i++)		// iterate only over the 1st-stage variables
	{
		if (fabs(ub[i] - lb[i]) > PHBB_IntTol)	return i;
	}
	
	return -1;
}

int bb::isIntegral(vector<double> &x)
{
	/*	Checks if the vector conforms the integrality restrictions,
		Returns the variable index that violates first.
		Returns -1 if the solution is feasible
	 */
	
	for (int i=0; i<x.size(); i++)
	{
		if (problem->integrality[i] &&
			fabs(round(x[i]) - x[i]) > PHBB_IntTol){
			return i;
		}
	}
	
	return -1;
}

int bb::strong_branch(node *cur_node)
{
	int picked_var = -1;
	
	double obj_diff = IloInfinity;
	
	for (int i=0; i<problem->numberOfVar[0]; i++)
	{
		if (problem->integrality[i] &&
			fabs(round(cur_node->prim_soln[i]) - cur_node->prim_soln[i]) > PHBB_IntTol)
		{
			// fractional solution found
			
			// create two nodes
			node *node_d = new node (cur_node->lowerbounds, cur_node->upperbounds);
			node *node_u = new node (cur_node->lowerbounds, cur_node->upperbounds);
			
			node_d->upperbounds[i] = 0.0;
			node_u->lowerbounds[i] = 1.0;
			
			node_d->lb = cur_node->lb;
			node_u->lb = cur_node->lb;
			
			node_d->ub = PHBB_PosInfty;
			node_u->ub = PHBB_PosInfty;
			
			node_d->obj_val = cur_node->obj_val;
			node_u->obj_val = cur_node->obj_val;
			
			node_d->prim_soln = cur_node->prim_soln;	// the primal solution will be projected in the LP solver for warmer start
			node_u->prim_soln = cur_node->prim_soln;
			
			node_d->dual_soln = cur_node->dual_soln;	// the dual solution is passed for a warmer start
			node_u->dual_soln = cur_node->dual_soln;

			solver.load(node_d);
			solver.optimize_w_bounds(1, inc_obj_val);
			
			if (obj_diff > fabs(cur_node->obj_val - node_d->obj_val))	{
				obj_diff = fabs(cur_node->obj_val - node_d->obj_val);
				picked_var = i;
			}

			solver.load(node_u);
			solver.optimize_w_bounds(1, inc_obj_val);
			
			if (obj_diff > fabs(cur_node->obj_val - node_u->obj_val))	{
				obj_diff = fabs(cur_node->obj_val - node_d->obj_val);
				picked_var = i;
			}

			delete node_d;
			node_d = NULL;
			
			delete node_u;
			node_u = NULL;
		}
	}
	
	if (picked_var == -1) {
 		picked_var = maxInfeasVar(cur_node->prim_soln);
	}
	
	return picked_var;
}

int bb::minInfeasVar (vector<double> &x)
{
	/*	Returns the variable with the minimum infeasibility.
			(*) 1-st stage variables have the highest priority
	 */
	
	int picked_var = -1;
	double mininfeasibility = 1.0, temp;
	for (int i=0; i<problem->numberOfVar[0]; i++)
	{
		if (problem->integrality[i]) {
			temp = fabs(round(x[i]) - x[i]);
			if ( temp > PHBB_IntTol && temp < mininfeasibility ) {
				mininfeasibility = temp;
				picked_var = i;
			}
		}
	}

	if ( mininfeasibility == 1.0 )	// 1st stage is integral
	{
		for (int i=problem->numberOfVar[0]; i<problem->numberOfVar[1] * problem->nb_scen; i++)
		{
			if (problem->integrality[i]) {
				temp = fabs(round(x[i]) - x[i]);
				if ( temp > PHBB_IntTol && temp < mininfeasibility ) {
					mininfeasibility = temp;
					picked_var = i;
				}
			}
		}
	}
	return picked_var;
}

int bb::maxInfeasVar (vector<double> &x)
{
	/*	Returns the variable with the maximum infeasibility.
	 (*) 1-st stage variables have the highest priority
	 */
	
	int picked_var = -1;
	double maxinfeasibility = 0.0;
	for (int i=0; i<problem->numberOfVar[0]; i++)
	{
		if (problem->integrality[i]) {
			if ( fabs(round(x[i]) - x[i]) > maxinfeasibility ) {
				maxinfeasibility = fabs(round(x[i]) - x[i]);
				picked_var = i;
			}
		}
	}
	
	if ( maxinfeasibility < PHBB_IntTol )	// 1st stage is integral
	{
		for (int i=problem->numberOfVar[0]; i<problem->numberOfVar[0] + problem->numberOfVar[1] * problem->nb_scen; i++)
		{
			if (problem->integrality[i]) {
				if ( fabs(round(x[i]) - x[i]) > maxinfeasibility ) {
					maxinfeasibility = fabs(round(x[i]) - x[i]);
					picked_var = i;
				}
			}
		}
		cout << "Branching on second stage" << endl;
//		cout << "Branching on second-stage: " << problem->vars[1][ picked_var % problem->nb_scen ].getName() << " (no: " << picked_var << ")" << endl;;
	}
	return picked_var;
}

void bb::printSolution()
{
	string output_fname = problem->probName + ".sol";
	ofstream output;
	output.open(output_fname.c_str());

	if ( soln_pool.size() == 0 )	{
		output << "No solution" << endl;
		output.close();
		return;
	}

	// copy the solutions (in case they are requested later)
	priority_queue<solution, vector<solution>> soln_pool_copy = soln_pool;
	
	// get the optimal solution (optimal solution is on the bottom of the list, because fathoming worse-objective solutions from the top is faster)
	while ( soln_pool_copy.size() > 1 )	soln_pool_copy.pop();
	
	solution opt_soln = soln_pool_copy.top();

	// print first-stage variables
	output << "STAGE 1" << endl;
	for (int k=0; k<problem->numberOfVar[0]; k++) {
		output << problem->vars[0][0][k].getName() << "\t" << opt_soln.vars[k] << endl;
	}

	// print second-stage variables
	output << "STAGE 2" << endl;
	for (int j=0; j<problem->numberOfVar[1]; j++)
	{
		output << problem->vars[0][1][j].getName();
		for (int s=0; s<problem->nb_scen; s++)
		{
			output << "\t" << opt_soln.vars[ problem->numberOfVar[0] + s * problem->numberOfVar[1] + j];
		}
		output << endl;
	}
	output.close();
}

void bb::displaySolutions()
{
	if ( soln_pool.size() == 0 )	{
		cout << "No solution" << endl;
		return;
	}

	priority_queue<solution, vector<solution>> soln_pool_copy = soln_pool;
	
	// store the obj values first
	vector<double> objs;
	while (!soln_pool_copy.empty()) {
		objs.push_back( soln_pool_copy.top().obj_val );
		soln_pool_copy.pop();
	}
	std::reverse( objs.begin(), objs.end() );
	
	printf("Soln          ObjVal\n");
	for (int k=0; k<objs.size(); k++)	printf("%4d  %14.6f\n", k, objs[k]);
}

void bb::add_node_cuts (node *cur_node)
{
	for (int c=0; c<cur_node->cut_indices.size(); c++)	problem->add_L2_feasibility_cut( cur_node->cut_indices[c] );
}

void bb::remove_node_cuts (node *cur_node)
{
	for (int c=0; c<cur_node->cut_indices.size(); c++)	problem->remove_L2_feasibility_cut(cur_node->cut_indices[c]);
}

void bb::setUpperCutoff(double val)
{
	inc_obj_val = val;
}

void bb::reset ()
{
	// empty the containers
	nodes.clear();
	lb_list.clear();

	while (priority_list.size() > 0)	priority_list.pop();
	while (soln_pool.size() > 0)		soln_pool.pop();
	
	// reset the node counter
	node::nodecount = 0;
}

double bb::getObjValue()
{
	return inc_obj_val;
}

double bb::getBestLowerBound()
{
	return lb_obj_val;
}

/******************************************************************************
 Prints algorithm performance statistics into the log file
 ******************************************************************************/
void bb::printLog ()
{
	ofstream output;
	output.open("output.log", ios::app);
	output << currentDateTime() << "\t" << problem->probName << "\t" << setprecision(3) << fixed << lb_obj_val << "\t" << inc_obj_val << "\t" << nb_nodes_processed << "\t" << solver.tot_nb_itr << "\t" << wall_t << "\t" << cpu_t << endl;
	output.close();
}

/******************************************************************************
 Displays small scripts if the optimization takes a bit long.
 ******************************************************************************/
void bb::displayMotivatingMessage ()
{
	switch ( rand() % 7 ) {
		case 0:
			cout << "Hang in there...";
			break;
		case 1:
			cout << "This seems like a tough problem!";
			break;
		case 2:
			cout << "Still working on it...";
			break;
		case 3:
			cout << "Give me a few more minutes, I'll be done soon.";
			break;
		case 4:
			cout << "This appears to be taking some time. Why don't you put a cup of tea?";
			break;
		case 5:
			cout << "I guess this problem is NP-Hard?";
			break;
		default:
			cout << "I'm almost there.";
			break;
	}
}

/******************************************************************************
 Displays algorithm performance statistics.
 Designed to be used after optimization ends.
 ******************************************************************************/
void bb::displayFinalMessage()
{
	printf("Obj value:             %.4f\n", inc_obj_val);
	printf("Wall / CPU time:       %.3f / %.3f\n\n", wall_t, cpu_t);
	
	printf("Nb of nodes processed: %d\n", nb_nodes_processed);
	printf("Nb of PH iterations:   %d\n", solver.tot_nb_itr);
	printf("Nb of solutions:       %d\n\n", (int)soln_pool.size());
}
