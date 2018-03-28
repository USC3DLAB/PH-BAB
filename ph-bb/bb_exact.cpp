//
//  bb_exact.cpp
//  ph-bb
//
//  Created by Semih Atakan on 10/7/16.
//  Copyright © 2016 3D Lab. All rights reserved.
//

#include "bb_exact.hpp"

bool bb_exact::optimize()
{
	wall_t	= get_wall_time();
	cpu_t	= get_cpu_time();
	
	cout << "             Nodes" << endl;
	cout << "     Node    Left     Incumbent\t     Bound\t    Rel Gap" << endl;
	
	no_impr_count = 0;
	nb_nodes_processed = 0;
	
	unsigned int cur_node_id = -1, time_interval = 1;
	bool status;
	
	/******************************************************************
	 Process the root node
	 
	 Currently, it solves the subproblems with no non-anticipativity, 
	 gets a lower bound, and stops. 
	 
	 Future work: Valid inequalities can be used for strengthening.
	 ******************************************************************/
	{
		node *root = get_next_node();
		solver.load(root);
		
		status = solver.optimize_root();
		if (!status) {		// infeasible root relaxation
			infeasible = true;
			delete root;
			root = NULL;
		}
		else {
			add_node(root);	// add it back for processing
		}
	}
	
	/******************************************************************
	 Branch-and-bound loop
	 ******************************************************************/
	while (!nodes.empty())
	{
		fathomed=false;
		infeasible=false;
		integral=false;
		unknown=false;
		
		// get the next node (node selection rule is determined in the data structure)
		node *cur_node = get_next_node();
		
		cur_node_id = cur_node->nodeid;
				
		// prepare the solver
		solver.load(cur_node);

		// add node-specific cuts
		add_node_cuts(cur_node);
		
		// solve the node
		status = solver.optimize_w_bounds_mip_option(PHBB_NodeItrLim, inc_obj_val);
		
		// remove node-specific cuts
		remove_node_cuts(cur_node);
		
		// process the node
		nb_nodes_processed++;

		/******************************/
		/*** For debugging purposes ***
		ofstream log;
		log.open("alg_progress.log", ios::app);

		log << "--" << endl;
		log << cur_node->getId() << " " << cur_node->lb << " " << cur_node->ub << " " << cur_node->obj_val << endl;
		for ( auto it = nodes.begin(); it != nodes.end(); ++it ) {
			log << it->first << " " << it->second->lb << " " << it->second->ub << " " << it->second->obj_val << endl;
		}
		log.close();
		 ******************************/

		if (!status)	// infeasible node
		{
			infeasible = true;
			delete cur_node;
			cur_node = NULL;
		}
		else			// feasible node
		{
			// check integrality
			int frac_index = isIntegral(cur_node->prim_soln);
			
			if (frac_index >= 0)	// fractional solution
			{
				if (frac_index < problem->numberOfVar[0])	// 1st-stage is still fractional
				{
					process_frac_node(cur_node);
				}
				else										// 2nd-stage is still fractional
				{
					// create two branches:
					//  (1) fix 1st-stage solns,
					//  (2) separate 1st-stage solns.
					fix_and_branch(cur_node);
				}
			}
			else					// integer-feasible solution
			{
				// if the node is optimal, it will be fathomed
				bool delete_cur_node = cur_node->isOptimal;
				
				// if the node's ub is not computed, PH objective value will be used
				if (cur_node->ub == INFINITY)	cur_node->ub = cur_node->obj_val;
				
				// process the node
				process_int_node(cur_node, delete_cur_node);
				
				integral = true;
				
				if (!fathomed && !delete_cur_node)	// if not optimal, this node needs to be reprocessed
				{
					unknown=true;
					add_node(cur_node);
				}
			}
		}
		
		// update the lower bound
		if (lb_list.size() > 0)		lb_obj_val = (lb_list.begin())->first;
		else						lb_obj_val = inc_obj_val;

		// Display log
		if (fathomed)			cout << "     ";
		else if (integral)		cout << "  *  ";
		else if (unknown)		cout << "  ?  ";
		else if (infeasible)	cout << "  x  ";
		else					cout << "     ";
		
		// current node, nb of nodes left
		cout << left << setw(8) << cur_node_id << setw(8) << nodes.size() << right;
		
		cout << fixed << setw(10) << setprecision(3);
		
		if (inc_obj_val == PHBB_PosInfty)	cout << "-" << "\t";
		else								cout << inc_obj_val << "\t";
		
		cout << fixed << setw(10) << setprecision(3);
		
		cout << lb_obj_val << "\t" << setw(10) << setprecision(2) << fixed;
		
		double rel_gap = min(100.0, fabs(lb_obj_val - inc_obj_val)/(fabs(inc_obj_val + PHBB_Eps))*100);
		if (inc_obj_val != PHBB_PosInfty && lb_obj_val != PHBB_NegInfty)	cout << rel_gap << "%" << endl;
		else																cout << "-" << endl;
				
		if (rel_gap/100.0 < PHBB_RelOptTol || fabs(lb_obj_val - inc_obj_val) < PHBB_AbsOptTol || get_wall_time()-wall_t > PHBB_TiLim )	break;
		
		// SPEAK!
		if ( get_wall_time() - wall_t > double(time_interval) * 120.0 ) {
			time_interval = int(get_wall_time() - wall_t)/120.0 + 1;
			cout << "Elapsed time = " << (get_wall_time() - wall_t) << ", ";
			displayMotivatingMessage();
			cout << endl;
		}
	}
	
	wall_t	= get_wall_time() - wall_t;
	cpu_t	= get_cpu_time() - cpu_t;
	
	cout << endl;
	
	displayFinalMessage();
	displaySolutions();
	
	cout << endl;
	
	printLog();
	
	return 0;
}
