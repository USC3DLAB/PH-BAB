//
//  ph.h
//  PH-BB
//
//  Created by Semih Atakan on 6/29/15.
//  Copyright (c) 2015 3D Lab. All rights reserved.
//

#ifndef __ph_bb__ph__
#define __ph_bb__ph__

#include <stdio.h>
#include <iomanip>

#include "commons.h"
#include "SPproblem.h"
#include "Node.h"
#include "config.h"
#include <fstream>
#include <vector>
#include <string>
#include <set>

#ifdef PHBB_Parallel
// the boost library is added for multi-threading
#include <boost/thread/thread.hpp>
#include <boost/bind.hpp>
#include <boost/asio.hpp>

#endif


class ph {
	
public:
	
	ph ();
	~ph();
	
	void initialize (SPprob *problem);

	bool optimize_root					();
	bool optimize_w_bounds_mip_option	(int itr_lim, double cutoff);
	bool optimize_w_bounds				(int itr_lim, double cutoff);
	bool optimize_w_mip_option			(double cutoff);
	bool optimize						();
	
	bool heuristic	();
	
	bool isNonanticipative();			// true if the 1st stage bounds are all equality
	bool isOptimal	();					// check the violation in the norm ||x_s - x||^2
	bool haveBoundsConverged();			// check the convergence of the bounds

	void load		(node *cur_node);	// load the solver with a node
	
	bool verbose, foutput;				// verbose=true displays output, foutput=true prints output
	
	unsigned int tot_nb_itr;			// total number of iterations

	double ph_t;
	
private:
	vector<double> inc;
	vector< vector<double> > x, y;
	vector< vector<double> > lambda;
	
	double opt_obj_val, obj_val_for_L2_cut;

	double obj_val, prev_obj_val;
	
	bool optimize_as_mip	(double &cutoff);
	bool update_primal_soln ();
	void update_dual_soln	();
	void update_lb_ub		();
	void update_inc			();

	// objects that are required for parallelization
	void update_primal_soln_scenario_subproblem (int s, double &obj);
	void update_lb_ub_scenario_subproblem		(int s, double &lb_obj, double &ub_obj);
	void optimize_root_scenario_subproblem		(int s, double &obj);
	void optimize_as_mip_full					(int s, double &bnode, double &obj);
	void optimize_as_mip_root					(int s, double &obj);
	
	bool check_termination_criteria ();

	double cutoff;		// usually the incumbent obj value, it is used to stop a solver from spending unnecessary time
	
	void printOutput	(double &lb, double &ub);
	void store_solution (bool &isOpt);
	
	ofstream ph_output;
	
	void stabilization (double &obj_val, double &prev_obj_val);
	int direction, nb_stability;
	
	double compute_dual_function_obj_value (vector< vector<double> > &x,
											vector< vector<double> > &y,
											vector< vector<double> > &lambda,
											int thread);
	double compute_dual_function_obj_value (vector< vector<double> > &x,
											vector< vector<double> > &y,
											vector< vector<double> > &lambda,
											vector< vector<double> > &prev_x,
											int thread);

	double compute_primal_obj_value	(IloNumArray &x, vector< IloNumArray > &y);
	double compute_subprob_primal_obj_value (vector< IloNumArray > &x, vector< IloNumArray > &y);
	
	double get_rel_gap	(const double &compare_this, double &wrt_this);

	
	double rho, prim_rho, coef;
	unsigned int itr, itr_lim;
	
	node	*nodeptr;	
	SPprob	*prob;
	
#ifdef PHBB_Parallel
	map<boost::thread::id, unsigned short> thread_map;
#endif

};

#endif /* defined(__dualMILP__ph__) */
