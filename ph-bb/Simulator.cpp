//
//  Simulator.cpp
//  ph-bb
//
//  Created by Semih Atakan on 10/6/17.
//  Copyright © 2017 3D Lab. All rights reserved.
//

#include "Simulator.hpp"

Simulator::Simulator (string probname)
{
	prob = new SPprob (probname);
	
	mem_alloc();
}

Simulator::Simulator (string probname, string sto_filename)
{
	prob = new SPprob (probname);
	
	prob->probName = sto_filename;	// override stochastic info
	prob->readSTO();
	
	mem_alloc();
}

Simulator::~Simulator ()
{
	delete prob;
}

void Simulator::mem_alloc()
{
	y.resize( prob->nb_scen );
	for (int s=0; s<prob->nb_scen; s++) y[s].resize(prob->numberOfVar[1]);
	
	obj_vals.resize( prob->nb_scen );
}

void Simulator::setFirstStageSoln (string sol_filename)
{
	/* Read solution */
	
	// open file
	ifstream input;
	sol_filename = sol_filename + ".sol";
	input.open(sol_filename.c_str());
	if (input.fail()) {
		cout << "Could not open "<< sol_filename << endl;
		exit(1);
	}

	// parse file
	string temp_str;
	double temp_dbl;
	
	input >> temp_str >> temp_dbl;	// the first tokens must be "STAGE 1"
	while (!input.eof()) {
		input >> temp_str >> temp_dbl;
		
		if (temp_str.compare("STAGE") != 0)	{
			firstStageSoln.insert( pair<string, double> (temp_str, temp_dbl) );
		}
		else {
			break;	// end of first-stage solution
		}
	}
	input.close();
	
	cout << "SOL file has been read successfully." << endl;
	
	/* fix the solution */
	for (int j=0; j<prob->vars[0][0].getSize(); j++)
	{
		temp_str = prob->vars[0][0][j].getName();
		
		for (int k=0; k<prob->nb_threads; k++) {
			prob->vars[k][0][j].setBounds(firstStageSoln[temp_str], firstStageSoln[temp_str]);
		}
	}
	
	/* ** ** **/
	// in case you want to fix some second-stage variables, write their solutions as if they are first-stage variables, and execute this loop
	for (int j=0; j<prob->vars[0][1].getSize(); j++)
	{
		temp_str = prob->vars[0][1][j].getName();
		
		if ( firstStageSoln.find(temp_str) == firstStageSoln.end() ) {
			// if this soln is not to be fixed, continue...
			continue;
		}
		else {
			// fix the second-stage variable
			for (int k=0; k<prob->nb_threads; k++) {
				prob->vars[k][1][j].setBounds(firstStageSoln[temp_str], firstStageSoln[temp_str]);
			}
		}
	}
	/* ** ** ** */
	
	cout << "1st-stage solution is fixed" << endl;
}

#ifdef PHBB_Parallel

void Simulator::simulate ()
{
	cout << "Simulation begins..." << endl;

	double			start_time = get_wall_time();
	
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
	for (int s=0; s<prob->nb_scen; s++) {
		io_service.post( boost::bind(&Simulator::optimize_subproblem_as_mip, boost::ref(*this), s) );
	}
	
	
	/***** Parallel programming stuff (START) *****/
	work.~work();		// let the io_service shutdown by removing the thread
	threads.join_all();
	/***** Parallel programming stuff (END)   *****/
	
	// Compute the objective value
	double obj_val = 0;
	for (int s=0; s<prob->nb_scen; s++)	{
		obj_val	+= prob->sceProb[s] * obj_vals[s];
	}

	cout << "Objective value= " << setprecision(4) << fixed << obj_val << endl;
	cout << "Simulation is completed [" << setprecision(1) << fixed << get_wall_time() - start_time << "s]" << endl;
}

/******************************************************************************
 Optimizes the subproblem as a mixed-integer program, to optimality.
 ******************************************************************************/
void Simulator::optimize_subproblem_as_mip (int s)
{
	unsigned short thread_no = thread_map[boost::this_thread::get_id()];
	
	prob->setup_scenario_subproblem (s, thread_no);
		
	bool status = prob->solve(thread_no);
	
//	if ((s+1)%5 == 0)	{ cout << "* Scenario " << s+1 << " is solved..." << endl; }
	
	if (status)	{
		obj_vals[s]	= prob->getObjValue(thread_no);
		prob->getSoln(y[s], 2, thread_no);
	}
	else {
		obj_vals[s]	= PHBB_PosInfty;
	}
}

#else

void Simulator::simulate ()
{
	cout << "Simulation begins..." << endl;
	
	unsigned short	thread_no = 0;
	bool			status = true;
	double			start_time = get_wall_time();

	for (int s=0; s<prob->nb_scen; s++)
	{
		prob->setup_scenario_subproblem (s, thread_no);
		
		status = prob->solve(thread_no);
		
		if (status)	{
			obj_vals[s] = prob->getObjValue(thread_no);
			prob->getSoln(y[s], 2, thread_no);
		}
		else {
			obj_vals[s] = PHBB_PosInfty;
			cout << "Infeasible subproblem" << endl;
			break;
		}
	}
	
	// Compute the objective value
	double obj_val = 0;
	for (int s=0; s<prob->nb_scen; s++)	{
		obj_val	+= prob->sceProb[s] * obj_vals[s];
	}
	
	cout << "Simulated objective value= " << setprecision(4) << fixed << obj_val << endl;
	cout << "Simulation is completed [" << setprecision(1) << fixed << get_wall_time() - start_time << "s]" << endl;
}

#endif

void Simulator::printSimulationSoln()
{
	string output_fname = prob->probName + "_simulation.sol";
	ofstream output;
	output.open(output_fname.c_str());
	
	// print first-stage variables
	output << "STAGE 1" << endl;
	for (std::map<string,double>::iterator it=firstStageSoln.begin(); it!=firstStageSoln.end(); ++it) {
		output << it->first << "\t" << it->second << endl;
	}
	
	// print second-stage variables
	output << "STAGE 2" << endl;
	for (int j=0; j<prob->numberOfVar[1]; j++)
	{
		output << prob->vars[0][1][j].getName();
		for (int s=0; s<prob->nb_scen; s++) {
			output << "\t" << y[s][j];
		}
		output << endl;
	}
	output.close();
}
