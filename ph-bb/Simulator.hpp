//
//  Simulator.hpp
//  ph-bb
//
//  Created by Semih Atakan on 10/6/17.
//  Copyright © 2017 3D Lab. All rights reserved.
//

#ifndef Simulator_hpp
#define Simulator_hpp

#include <iostream>
#include <fstream>
#include <string>

#include "SPproblem.h"

#ifdef PHBB_Parallel	// Boost Library is linked for multi-threading
#include <boost/thread/thread.hpp>
#include <boost/bind.hpp>
#include <boost/asio.hpp>
#endif


class Simulator {
	
public:
	Simulator (string probname);
	Simulator (string probname, string sto_filename);
	~Simulator ();
	
	void	simulate ();
	void	setFirstStageSoln (string sol_filename);
	void	printSimulationSoln ();
	
private:
	SPprob				*prob;
	map<string, double>	firstStageSoln;		// first-stage solution retrieved from SOL file

	vector< vector<double> >	y;			// second-stage solution from the simulator
	vector< double >			obj_vals;	// scenario subproblem objective value
	
	void	mem_alloc ();
	
#ifdef PHBB_Parallel
	void	optimize_subproblem_as_mip (int s);
	map<boost::thread::id, unsigned short> thread_map;
#endif

};

#endif /* Simulator_hpp */
