//
//  Node.h
//  PH-BB
//
//  Created by Semih Atakan on 6/6/15.
//  Copyright (c) 2015 3D Lab. All rights reserved.
//

#ifndef __ph_bb__Node__
#define __ph_bb__Node__

#include <stdio.h>
#include <vector>
#include "config.h"
using namespace std;


/* Class node
 
 Element of a partition of the initial variable bounds bounds. The node is defined by upper & lower bounds of the integral variables.
 
 */

class node {
	
	friend class bb;
	
public:
	
	node	();
	node	(vector<double> lbs, vector<double> ubs);
	~node	();
	
	int getId	();
	
	unsigned int nodeid;	// identification number of the node
	
	vector<double> lowerbounds;
	vector<double> upperbounds;
	
	vector<int> cut_indices;	// the IDs of node-specific cuts
	
	vector<double> prim_soln;				// the final primal solution associated with the node
	vector< vector<double> > dual_soln;		// the final dual solution associated with the node
	
	double lb		= PHBB_NegInfty;	// a guaranteed lower bound
	double ub		= PHBB_PosInfty;	// a guaranteed upper bound
	double obj_val	= 0.0;				// PH objective value
		
	
	bool isOptimal;				// guaranteed optimality flag
	
	unsigned int nb_of_passes;	// count how many times the node has been processed
	
	double last_rho, last_coef;
	
private:
	static int nodecount;	// a counter that is incrementally increased whenever the node constructor is called

};

#endif /* defined(__ph_bb__Node__) */
