//
//  main.cpp
//  PH-BB
//
//  Created by Semih Atakan on 6/6/15.
//  Copyright (c) 2015 3D Lab. All rights reserved.
//

#include <iostream>
#include <vector>

using namespace std;

#include "SPproblem.h"
#include "bb_approx.hpp"
#include "bb_exact.hpp"
#include "ph.h"
#include "Simulator.hpp"

string probname = "./sslp/sslp_5_50_100", sim_filename = "", algorithm = "phbab";

bool read_program_inputs(int argc, const char * argv[])
{
	// read program inputs
	for (int a=1; a<argc; a++) {
		if (strcmp(argv[a], "-f") == 0) {
			probname = argv[++a];
		} else if (strcmp(argv[a], "-s") == 0) {
			sim_filename = argv[++a];
		} else if (strcmp(argv[a], "-a") == 0) {
			algorithm = argv[++a];
		}
	}

	// error check
	if (algorithm.compare("phbab") != 0 && algorithm.compare("phbab-apx") != 0) {
		cout << "Error: Wrong algorithm name is specified" << endl;
		return false;
	}
	
	// read out loud the program inputs
	cout << "Problem:               " << probname << endl;
	if (sim_filename.compare("") != 0)	cout << "Simulation Stoch File: " << sim_filename << endl;
	cout << "Algorithm:             " << algorithm << endl << endl;;
	
	return true;
}

int main(int argc, const char * argv[])
{
	/* read program inputs */
	bool status = read_program_inputs(argc, argv);
	if (!status) {
		return -1;
	}
	
	/* Create the Stochastic Problem */
	SPprob problem (probname);

	/* Optimize the Stochastic Problem */
	if (algorithm.compare("phbab") == 0) {
		// PHBAB: PH solves are partially completed. Lower bounding LPs are solved to determine
		// a guaranteed lower-bound on the optimal objective value of the node.
		bb_exact solver;
		
		// load the problem
		solver.load( problem );
		
		// optimize the problem
		solver.optimize();
		
		// print the solution
		solver.printSolution();
	}
	else if (algorithm.compare("phbab-apx") == 0) {
		// PHBAB-Apx: PH solves are done up to a numerical tolerance and the objective value is assumed to be optimal
		// Warning: Loose numerical leads to inexact-pruning of the branch-and-bound tree.
		bb_approx solver;
		
		// load the problem
		solver.load( problem );
		
		// optimize the problem
		solver.optimize();
		
		// print the solution
		solver.printSolution();
	}
	
	/* Simulate a Solution */
	if (sim_filename.compare("") != 0) {
		// Simulator: By providing an additional "validation" STOC file (ideally involving thousands of scenarios), one can simulate
		// the "optimal" first-stage decisions, obtained above.
		Simulator simulator (probname, sim_filename);
		
		// read first-stage solution from the solution file (printed above)
		simulator.setFirstStageSoln (probname);
		
		// execute the simulation
		simulator.simulate();
		
		// print simulation solution
		simulator.printSimulationSoln();
	}

	return 0;
}


