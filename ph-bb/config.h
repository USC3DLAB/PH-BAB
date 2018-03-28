//
//  config.h
//  PH-BB
//
//  Created by Semih Atakan on 6/7/15.
//  Copyright (c) 2015 3D Lab. All rights reserved.
//

#ifndef ph_bb_config_h
#define ph_bb_config_h

#include <limits>

#define PHBB_PosInfty			std::numeric_limits<double>::infinity()		// + infinity
#define PHBB_NegInfty			-std::numeric_limits<double>::infinity()	// - infinity
#define PHBB_Eps				1e-10		// epsilon for numerical corrections

#define PHBB_SolnPoolRelOptTol	1e-1		// solns with this percent optimality gap will be kept

#define PHBB_IntTol				1e-5		// tolerance for integrality checks
#define PHBB_RelOptTol			1e-4		// tolerance for relative optimality checks
#define PHBB_AbsOptTol			1e-6		// tolerance for absolute optimality checks

#define PHBB_MIPRelOptTol		1e-4		// relative optimality tolerance for MIP solves
#define	PHBB_MIPAbsOptTol		1e-6		// absolute optimality tolerance for MIP solves 

#define PHBB_NormTol			1e-1
#define PHBB_NodeRelOptTol		1e-6
#define PHBB_NodeItrLim			1
#define PHBB_NodeBndItrLim		1			// per how many iterations the lower/upper bounding procedures for a node should be invoked

#define PHBB_rho				1

#define PHBB_TiLim				3600		// time limit

#define PHBB_Parallel						// comment out if serial implementation is desired

#ifdef PHBB_Parallel
#define PHBB_Threads			4
#else 
#define PHBB_Threads			1
#endif

#endif
