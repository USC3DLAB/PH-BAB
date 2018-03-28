//
//  bb.h
//  PH-BB
//
//  Created by Semih Atakan on 6/7/15.
//  Copyright (c) 2015 3D Lab. All rights reserved.
//

#ifndef __ph_bb__bb__
#define __ph_bb__bb__

#include <stdio.h>
#include <vector>
#include <queue>
#include <ctime>
#include <tr1/unordered_map>
#include <algorithm>
#include <map>

#include "commons.h"
#include "SPproblem.h"
#include "ph.h"
#include "Node.h"
#include "config.h"

/******************* Preliminaries ************************/

// pairs of some values and the associated node IDs
struct double_uint_pair {
	double value;
	unsigned int nodeid;
 
	double_uint_pair (double val, unsigned int no)
	{
		value = val;	nodeid = no;
	}
};

// comparison objects
class larger_value_wins {
	public:
	bool operator() (const node &lhs, const node &rhs) const
	{
		return lhs.ub < rhs.ub;
	}
};

class smaller_value_wins {
	public:
	bool operator() (const double_uint_pair &lhs, const double_uint_pair &rhs) const {
		return lhs.value > rhs.value;
	}
};

/******************* end of preliminaries ************************/


class bb {
	
public:
	
	bb();
	
	void load				(SPprob &prob);
	void reset				();
	virtual bool optimize	();
	
	void setUpperCutoff		(double val);	// sets the inc_obj_val so that anything worse can be fathomed

	void printSolution				();
	void displaySolutions			();
	
	double getObjValue			();
	double getBestLowerBound	();
	
	double inc_obj_val;
	vector<double> inc_soln;
	
	double lb_obj_val;
	vector<double> lb_soln;
	
protected:
	SPprob *problem;	// the problem object
	
	ph solver;		// the solver object
	
	void split (node *cur_node, int &frac_index);
	
	int isIntegral		(vector<double> &x);
	int isFixed			(vector<double> &lb, vector<double> &ub);
	
	int minInfeasVar	(vector<double> &x);
	int maxInfeasVar	(vector<double> &x);
	int strong_branch	(node *cur_node);
	
	double get_rel_gap		(const double &compare_this, double &wrt_this);
	double get_rel_abs_gap	(const double &compare_this, double &wrt_this);
	
	void add_node			(node *cur_node);
	void process_int_node	(node *cur_node, bool delete_cur_node);
	void process_int_node	(node *cur_node);
	void process_frac_node	(node *cur_node);
	void fix_and_branch		(node *cur_node);

	bool rounding_heuristic (node *cur_node);
	
	node *get_next_node		();
	
	void printLog					();
	void displayMotivatingMessage	();
	void displayFinalMessage		();

	unsigned int no_impr_count;
	unsigned int nb_nodes_processed;
	
	void add_node_cuts		(node *cur_node);	// these store the L2 cuts which we use in the specialized PH-BAB algorithm
	void remove_node_cuts	(node *cur_node);
	
	/*********************************** Node lists ***********************************
	 - nodes:			contains nodes waiting for (further) processing
	 - soln_pool:		contains primal solutions within specified opt tolerances
	 
	 - priority_list:	sorts the node IDs according to some branching rule
	 - lb_list:			contains a list of lower bounds and ID of nodes claiming ownership	*/
	
	tr1::unordered_map<unsigned int, node *> nodes;
	
	struct solution {
		vector<double> vars;
		double obj_val;

		solution (vector<double> &vars, double &obj_val) {
			this->vars		= vars;
			this->obj_val	= obj_val;
		}
		
		bool operator < (solution other) const {
			return obj_val < other.obj_val;
		}
		bool operator > (solution other) const {
			return obj_val > other.obj_val;
		}
	};
	priority_queue<solution, vector<solution>> soln_pool;
	
	priority_queue<double_uint_pair, vector<double_uint_pair>, smaller_value_wins>	priority_list;
	multimap<double, unsigned int>	lb_list;
	
	/******************************* end of node lists *********************************/
		
	// logging
	bool fathomed, integral, infeasible, unknown;
	
	// timing
	double wall_t, cpu_t;
};
#endif /* defined(__ph_bb__bb__) */
