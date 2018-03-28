//
//  randVar.cpp
//  ph-bb
//
//  Created by Semih Atakan on 8/17/17.
//  Copyright © 2017 3D Lab. All rights reserved.
//

#include "RandVar.hpp"

RandVar::RandVar() {}

void RandVar::set_cumulative_probs() {
	cuProb.resize(Prob.size());
	cuProb[0] = Prob[0];
	for (int i=1; i<Prob.size()-1; i++)    {
		cuProb[i] = cuProb[i-1] + Prob[i];
	}
	cuProb[ cuProb.size()-1 ] = 1;
}

bool RandVar::isInRHS() {
	return loc == RandVar::rhs;
}

bool RandVar::isInObj() {
	return loc == RandVar::obj;
}

bool RandVar::isInMat() {
	return loc == RandVar::mat;
}

IloInt RandVar::getColIndex() {
	return colIndex;
}

IloInt RandVar::getRowIndex() {
	return rowIndex;
}

double RandVar::getValue(int &s) {
	return Value[s];
}
