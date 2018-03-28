//
//  cpx.hpp
//  ph-bb
//
//  Created by Semih Atakan on 3/13/17.
//  Copyright © 2017 3D Lab. All rights reserved.
//

#ifndef cpx_hpp
#define cpx_hpp

#include <string>
#include <vector>

#include "config.h"
#include "ilcplex/ilocplex.h"

using namespace std;

/* Class for CPLEX objects */
class cpx {
	
	friend class SPprob;
	friend class Simulator;
	friend class ph;

public:
	cpx ();
	
	void initialize				();
	void importModel			(string fname);
	void extractModel			();
	void convertToLP			();
	void convertToMIP			(vector<bool> &integrality);
	void end					();								// must call to release CPX memory
	void setDefaultParameters	();
	
	string getProblemType		();
	
	

private:
	IloEnv		env;
	IloModel	model;
	IloCplex	cplex;
	
	IloObjective	obj;		// objective function
	IloNumVarArray	var;		// variables
	IloRangeArray	rng;		// constraints

	IloConversion	ConvToLP, ConvToMIP;	// for converting the problem to MIP or LP
};


#endif /* cpx_hpp */
