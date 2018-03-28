//
//  cpx.cpp
//  ph-bb
//
//  Created by Semih Atakan on 3/13/17.
//  Copyright © 2017 3D Lab. All rights reserved.
//

#include "cpx.hpp"

cpx::cpx() { }

void cpx::end() {
	env.end();
}

void cpx::initialize()
{
	env	  = IloEnv();
	model = IloModel (env);
	cplex = IloCplex (env);
	
	obj = IloObjective	 (env);
	var = IloNumVarArray (env);
	rng = IloRangeArray  (env);

	setDefaultParameters();
}

void cpx::importModel(string fname)
{
	try {
		cplex.importModel(model, fname.c_str(), obj, var, rng);
	}
	catch (IloException &e) {
		cout << e << endl;
		exit(1);
	}
}

string cpx::getProblemType()
{
	// get problem type
	if (cplex.isMIP())	{
		if (cplex.isQC())		return "miqcp";
		else if (cplex.isQO())	return "miqp";
		else					return "milp";
	}
	else {
		if (cplex.isQC())		return "qcp";
		else if (cplex.isQO())	return "qp";
		else					return "lp";
	}
	return "error";
}

void cpx::extractModel()
{
	try {
		cplex.extract(model);		// prepare the solver
	}
	catch (IloException &e) {
		cout << e << endl;
		exit(1);
	}
}

void cpx::setDefaultParameters()
{
	/* Default Parameters */
	cplex.setOut(env.getNullStream());
	cplex.setWarning(env.getNullStream());
	
	cplex.setParam(IloCplex::RootAlg,	IloCplex::Algorithm::Dual);
	cplex.setParam(IloCplex::EpGap,		PHBB_MIPRelOptTol);
	cplex.setParam(IloCplex::EpAGap,	PHBB_MIPAbsOptTol);
	cplex.setParam(IloCplex::TiLim,		PHBB_TiLim);
	cplex.setParam(IloCplex::Threads,	1);
}

void cpx::convertToLP()
{
	// remove any old conversion first (may not exist, and hence, we need try-catch)
	try {
		model.remove(ConvToMIP);
		ConvToMIP.end();
	}
	catch (...) {}
	
	// create and add the new conversion
	ConvToLP = IloConversion(env, var, ILOFLOAT);
	model.add( ConvToLP );
}

void cpx::convertToMIP(vector<bool> & integrality)
{
	// remove any old conversion first (may not exist, and hence, we need try-catch)
	try {
		model.remove(ConvToLP);
		ConvToLP.end();
	}
	catch (...) {}
	
	// create and add the new conversion
	IloNumVarArray temp (env);
	for (int j=0; j<var.getSize(); j++) {
		if (integrality[j])	temp.add( var[j] );
	}
	
	ConvToMIP = IloConversion (env, temp, ILOINT);
	model.add( ConvToMIP );
}

