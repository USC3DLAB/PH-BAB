//
//  Node.cpp
//  PH-BB
//
//  Created by Semih Atakan on 6/6/15.
//  Copyright (c) 2015 3D Lab. All rights reserved.
//

#include "Node.h"
#include <iostream>

int node::nodecount = 0;	// initializer for the static variable

node::node()
{
	nodeid = nodecount++;
	isOptimal = false;
	nb_of_passes = 0;
	
	last_rho = -1;
}

node::~node() {}

node::node(vector<double> lbs, vector<double> ubs)
{
	nodeid = nodecount++;
	lowerbounds = lbs;
	upperbounds = ubs;
	isOptimal = false;
	nb_of_passes = 0;
	
	last_rho = -1;
	
	//	cout << "Node " << nodeid << " created." << endl;
}

int node::getId()
{
	return nodeid;
}
