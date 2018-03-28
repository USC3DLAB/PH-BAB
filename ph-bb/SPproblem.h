//
//  SPproblem.h
//  PH-BB
//
//  Written by Yifan Liu on 02/05/15,
//	Revised by Semih Atakan on 06/28/15
//	Revised by Semih Atakan on 09/13/15
//
//  Copyright (c) 2015 3D Lab. All rights reserved.
//

#ifndef __ph_bb__SPproblem__
#define __ph_bb__SPproblem__

#define nb_stages 2

#include <stdio.h>
#include <string>
#include <vector>
#include <map>

#include "config.h"
#include "RandVar.hpp"
#include "cpx.hpp"

#include "ilcplex/ilocplex.h"
#include "commons.h"
using namespace std;

class SPprob {
	
	friend class Simulator;
public:
	
	SPprob	(string filename);
	~SPprob ();

	bool	solve			(int solver_index);
	double	getObjValue		(int solver_index);
	void	getSoln		(vector<double> &soln, int stage, int solver_index);		// extract the solution of variables in a particular stage from the corresponding solver
	
	void	convert_to_LP ();
	void	convert_to_LP (int k);
	void	convert_to_MIP();
	void	convert_to_MIP(int k);

	void	use_linear_objective	(int solver_index);		// non-regularized objective with no dual multipliers (for the 1st itr of PH)
	void	use_saddle_objective	(vector<double> &lambda, int solver_index);	// non-regularized objective with dual multipliers
	
	void	add_ph_regularizer		(double &rho, int solver_index);	// first-stage regularization terms are added
	void	remove_ph_regularizer	(int solver_index);					// first-stage regularization terms are removed

	double	setup_regularized_subproblem	(int &s, vector<double> &inc, vector<double> &lambda, double &rho, int solver_index);
	void	setup_scenario_subproblem		(int &s, int solver_index);
	void	set_random_rhs					(int &s, int solver_index);
	
	void	change_t1var_bounds (vector<double> &new_lbs, vector<double> &new_ubs);
	void	change_t1var_bounds (vector<double> &new_lbs, vector<double> &new_ubs, int &solver_index);
	void	change_t2var_bounds (vector<double> &new_lbs, vector<double> &new_ubs, int &s, int solver_index);
	
	long	derive_L2_feasibility_cut	(vector<double> &soln, vector<int> &fixed_var_indices);
	void	add_L2_feasibility_cut		(int &index);
	void	remove_L2_feasibility_cut	(int &index);
	
	void	displayProblemStats	();
	
	/* Problem Attributes */
	string probName;				// problem name
	string probType;				// problem type (lp, qp, mip)

	vector<double>	lb;				// variable lower bounds
	vector<double>	ub;				// variable upper bounds
	vector<bool>	integrality;	// true if the variable is integer
	vector<char>	con_type;		// e: equality, g: greater or equal, l: less or equal

	vector< double > sceProb;					// a vector for storing scenario probability

	int		nb_scen;
	long	nb_vars;		// total nb of variables in each scenario (first + second stages)
	int		nb_rvars;		// nb of random variables
	int		nb_rvars_obj, nb_rvars_mat, nb_rvars_rhs;
	int		nb_threads;

	vector< IloInt > numberOfVar;			// number of variables at each stage
	
	vector< RandVar > randomVars;			// independent random variables from stoch file
	
	// CPX solver(s)
	vector<cpx> solver;

	vector< vector<IloNumArray> >		obj_coefs;	// stores linear obj coefs for each solver and each stage
	vector< vector<IloNumVarArray> >	vars;		// variables are partitioned into 1st and 2nd stages (for ease of access)
	vector< IloRangeArray >				node_cuts;	// storage for node specific cuts
	
private:

	/* READ and STORE INPUT */
	void initialize();						// initialize the enviroment and populate all necessary data
	
	void readCOR ();						// read .cor file
	void readTIM ();						// read .tim file
	void readSTO ();						// read .sto file
	
	void store_constraint_info	();			// sets con_type
	void store_variable_info	();			// sets lb, ub, integrality vectors
	void store_objective_info	();			// sets obj_coef

	void get_linear_obj_coefs	(int k);	// collects the linear objective function coefficients
	void assign_vars_to_stages	();			// ex: 2nd stage variables will be assigned to vars[2]
	
	map<string, long>			var2index;	// returns the index of the desired variable
	map<string, long>			con2index;	// returns the index of the desired constraint
	
	vector< vector<string> >	time_info;	// store the information of the time 
};

#endif /* defined(__ph_bb__SPproblem__) */
