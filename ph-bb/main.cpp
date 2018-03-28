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


int main(int argc, const char * argv[])
{
	// initializations
	string probname = "./sslp/sslp_5_25_50";
	
	/* Create the Stochastic Problem */
	SPprob problem (probname);

	/* Optimize the Stochastic Problem */
	// PHBAB-Apx: PH solves are done up to a numerical tolerance and the objective value is assumed to be optimal
	// Warning: Loose numerical leads to inexact-pruning of the branch-and-bound tree.
	// bb_approx solver;
	
	// PHBAB: PH solves are partially completed. Lower bounding LPs are solved to determine
	// a guaranteed lower-bound on the optimal objective value of the node.
	bb_exact solver;
	
	// load the problem
	solver.load( problem );
	
	// optimize the problem
	solver.optimize();
	
	// print the solution
	solver.printSolution();
	
	/* Simulate a Solution *
	// Simulator: By providing an additional "validation" STOC file (ideally involving thousands of scenarios), one can simulate
	// the "optimal" first-stage decisions, obtained above.	
	Simulator simulator (probname, sto_filename);
	
	// read first-stage solution from the solution file (printed above)
	simulator.setFirstStageSoln (probname);
	
	// execute the simulation
	simulator.simulate();
	
	// print simulation solution
	simulator.printSimulationSoln();
	*/
	
	return 0;
}


