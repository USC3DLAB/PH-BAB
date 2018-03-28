//
//  randVar.hpp
//  ph-bb
//
//  Created by Semih Atakan on 8/17/17.
//  Copyright © 2017 3D Lab. All rights reserved.
//

#ifndef randVar_hpp
#define randVar_hpp

#include <string>
#include <vector>
#include "ilcplex/ilocplex.h"
using namespace std;

/* Class for random variables */
class RandVar {

	friend class SPprob;
	friend class Simulator;
	friend class cut_gen;

public:
	
	RandVar ();
	
	void set_cumulative_probs ();	// Calculates the cumulative probabilites (invoke after reading the STOCH file)
	
	bool isInRHS ();
	bool isInObj ();
	bool isInMat ();
	
	IloInt getColIndex();
	IloInt getRowIndex();
	double getValue(int &s);
	
private:
	
	string colName;								// variable's col name
	string rowName;								// variable's row name
	IloInt rowIndex = -1;						// variable's row index (-1 if in the objective)
	IloInt colIndex = -1;						// variable's col index (-1 if in the right-hand-side)
	
	vector<double> Value;						// realizations
	vector<double> Prob;						// probabilities
	vector<double> cuProb;						// cumulative probabilities
	
	enum rVarLoc { obj, rhs, mat };				// location of randomness
	rVarLoc loc;
};

#endif /* randVar_hpp */
