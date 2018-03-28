//
//  SPproblem.cpp
//  PH-BB
//
//  Created by Semih Atakan on 6/28/15.
//  Copyright (c) 2015 3D Lab. All rights reserved.
//

#include "SPproblem.h"

SPprob::SPprob (string fileName)
{
	probName = fileName;
	
	initialize();
}

SPprob::~SPprob ()
{
	for (int k=0; k<nb_threads; k++)	solver[k].end();
}

void SPprob::initialize()
{
	nb_threads = PHBB_Threads;

	// initialize CPX objects
	solver.resize(nb_threads);
	for (int k=0; k<nb_threads; k++)	solver[k].initialize();
	
	// read the SMPS file
	readCOR();
	readTIM();
	readSTO();

	// miscellaneous initializations
	store_constraint_info	();
	store_variable_info		();
	store_objective_info	();

	// add variables to the formulation (so that you can extract them even if they are presolved to be 0)
	for (int k=0; k<nb_threads; k++)	solver[k].model.add( solver[k].var );

	// initialize node-specific cut container for each solver
	node_cuts.resize(nb_threads);
	for (int k=0; k<nb_threads; k++)	node_cuts[k] = IloRangeArray ( solver[k].env );

	probType = solver[0].getProblemType();
	
	// extract the model and ready the LP/MIP solver
	for (int k=0; k<nb_threads; k++)	solver[k].extractModel();

	displayProblemStats();
	cout << "CPLEX version " << solver[0].cplex.getVersion() << endl;
}

void SPprob::store_objective_info ()
{
	// store linear obj coefs

	obj_coefs.resize(nb_threads);
	for (int k=0; k<nb_threads; k++)
	{
		obj_coefs[k].resize(nb_stages);
		for (int t=0; t<nb_stages; t++) obj_coefs[k][t] = IloNumArray( solver[k].env, numberOfVar[t] );
		
		get_linear_obj_coefs (k);
	}
}

void SPprob::assign_vars_to_stages()
{
	vars.resize(nb_threads);
	for (int k=0; k<nb_threads; k++) {
		vars[k].resize(nb_stages);
		for (int t=0; t<nb_stages; t++)	vars[k][t] = IloNumVarArray ( solver[k].env );
		
		for (int i=0; i<nb_vars; i++) {
			if ( var2index[solver[k].var[i].getName()] < var2index[time_info[1][0]] )	{
				vars[k][0].add( solver[k].var[i] );
			}
			else {
				vars[k][1].add( solver[k].var[i] );
			}
		}
	}
	
	numberOfVar.resize(2);
	numberOfVar[0] = vars[0][0].getSize();
	numberOfVar[1] = vars[0][1].getSize();
}

void SPprob::readCOR()
{
	string corFileName = probName + ".cor";
	
	// initialize one cplex object per thread
	for (int k=0; k<nb_threads; k++)	solver[k].importModel(corFileName);
	
	nb_vars = solver[0].var.getSize();
	
	cout << "CORE file has been read successfully." << endl;
}

void SPprob::store_constraint_info()
{
	/* Store Constraint Types *
	 
	 Constraints have the form LB <= Ax <= UB,
	 
	 if LB = UB are not infinity, then,
	 we have Ax = LB (= UB),
	 else if (only) UB is not infinity, then,
	 we have Ax <= UB,
	 else if (only) LB is not infinity, then,
	 we have Ax >= LB.
	 */
	
	con_type = vector<char>	( solver[0].rng.getSize() );
	
	for (int c=0; c<solver[0].rng.getSize(); c++) {
		if (solver[0].rng[c].getUB() != IloInfinity && solver[0].rng[c].getLB() != -IloInfinity) {
			con_type[c] = 'e';
		}
		else if (solver[0].rng[c].getUB() != IloInfinity) {
			con_type[c] = 'l';
		}
		else {
			con_type[c] = 'g';
		}
	}
}

void SPprob::readTIM()
{
	string timFileName = probName + ".tim";
	ifstream timFileStream;
	timFileStream.open(timFileName.c_str());
	if (timFileStream.fail()) {
		cout << "Could not open "<< timFileName << endl;
		exit(1);
	}

	string line;
	int numberOfStage = 0;
	while (getline(timFileStream, line)) {
		if (line[0] == '*' || line[0] == '\r' || (line[0] >= '0' && line[0] <= 'Z'))
			continue; // skip comment, empty line and basically all other stuff
		
		/* if we are here, chances are that we meet the line
		 * which stores stage decomposition information
		 */
		stringstream lineStream(line);		// cast this string into a stream so that we can obtain its component
		string field1, field2;				// field1 stores the column's name
		lineStream >> field1 >> field2;		// field2 stores the row's name
		
		time_info.push_back( vector<string> (2) );
		
		time_info[numberOfStage][0] = field1;
		time_info[numberOfStage][1] = field2;

		numberOfStage++;							  // keep track of the number of stage
	}
	
	if (numberOfStage != nb_stages) {
		cout << "Error reading "<< timFileName <<"! Two stage needed but "<< numberOfStage << " stage found";
	}
	
	timFileStream.close();
	cout << "TIME file has been read successfully." << endl;
}

void SPprob::readSTO()
{
	/* STO file is read in SCENARIOS format */

	/* Open STOCH file */
	string stoFileName = probName + ".sto";
	ifstream stoFileStream;
	stoFileStream.open(stoFileName.c_str());	
	if (stoFileStream.fail()) {
		cout << "Could not open "<< stoFileName << endl;
		exit(1);
	}
	
	// reset containers and parameters
	randomVars.clear();
	sceProb.clear();
	
	nb_scen = 0;		// remember, if you don't initialize this, you will loose 6 hours of sleep.
	nb_rvars = 0;
	nb_rvars_obj = 0;
	nb_rvars_rhs = 0;
	nb_rvars_mat = 0;
	
	int var_index = 0;

	double sumOfProbs=0, prob=0;
	string line;
	while ( !safeGetline(stoFileStream, line).eof() ) {
		if (line[0] == '*' || line[0] == '\r' || (line[0] >= '0' && line[0] <= 'Z')) {
			continue; // skip comment, empty line and basically all other stuff
		}
		
		/* if we are here, chances are that we meet the line
		 * which stores stochastic information
		 */
		
		stringstream lineStream(line);
		lineStream >> line;
		if (line.compare("SC") == 0) {	// we are getting a new scenario
			nb_scen++;
			
			lineStream >> line >> line >> line;
			prob = atof(line.c_str());
			
			sumOfProbs += prob;
			sceProb.push_back(prob);
			
			var_index = 0;
		}
		else {							// we are getting the contents of the new scenario
			
			if (fabs(sumOfProbs - prob) <= PHBB_Eps) { // then this is the 1st scenario, the random vars are created here
				nb_rvars++;
				randomVars.push_back(RandVar());
				
				randomVars[nb_rvars - 1].colName = line;
				lineStream >> randomVars[nb_rvars - 1].rowName >> line;
				
				randomVars[nb_rvars - 1].Value.push_back(atof(line.c_str()));
				randomVars[nb_rvars - 1].Prob.push_back(prob);
			}
			else {
				lineStream >> line >> line;
				
				randomVars[ var_index ].Value.push_back( atof(line.c_str() ) );
				randomVars[ var_index ].Prob.push_back(prob);
			}
			var_index++;
		}
	}
	
	if ( sumOfProbs < 1.0 - PHBB_Eps ) {
		cout << "Probabilities of the scenarios do not add up to 1.0 (" << sumOfProbs << ")" << endl;
		exit(1);
	}
	
	/* Find the random variable locations, Set their indices, Count RVs location-wise */
	int nb_deleted_rvar = 0;
	for (int i = 0; i < nb_rvars; i++) {
		if ( randomVars[i].colName.compare("RHS") == 0 )
		{
			// the Rvar is in the right-hand-side
			randomVars[i].loc = RandVar::rhs;
			
			// the column index is -1
			randomVars[i].colIndex = -1;
			
			// find the row index
			int c;
			for (c=0; c<solver[0].rng.getSize(); c++)	{
				if ( randomVars[i].rowName.compare(solver[0].rng[c].getName()) == 0 ) {
					randomVars[i].rowIndex = c;
					break;
				}
			}
			
			// remove the RV that cannot be located
			if (c == solver[0].rng.getSize()) {
				if (nb_deleted_rvar <= 5) {
					cout << "Warning: Random variable " << randomVars[i].colName << " " << randomVars[i].rowName << " is not in the model" << endl;
				}
				randomVars.erase( randomVars.begin() + i );
				nb_rvars--;		// nb_rvars decreases by 1
				i--;			// iterate should not proceed
				nb_deleted_rvar++;
				continue;
			}
			
			nb_rvars_rhs++;
		}
		else if ( randomVars[i].rowName.compare(solver[0].obj.getName()) == 0 )
		{
			// the Rvar is in the objective
			randomVars[i].loc = RandVar::obj;
			
			// the row index is -1
			randomVars[i].rowIndex = -1;
			
			// find the column index
			for (int j=0; j<solver[0].var.getSize(); j++) {
				if ( randomVars[i].colName.compare(solver[0].var[j].getName()) == 0 ) {
					randomVars[i].colIndex = j;
					break;
				}
			}
			nb_rvars_obj++;
		}
		else {
			cout << "Error: This code cannot handle random variables affecting the coefficient matrices!";
			
			// the Rvar is in the matrix
			randomVars[i].loc = RandVar::mat;
			
			// find the row index
			for (int c=0; c<solver[0].rng.getSize(); c++)	{
				if ( randomVars[i].rowName.compare(solver[0].rng[c].getName()) == 0 ) {
					randomVars[i].rowIndex = c;
					break;
				}
			}
			
			// find the column index
			for (int j=0; j<solver[0].var.getSize(); j++) {
				if ( randomVars[i].colName.compare(solver[0].var[j].getName()) == 0 ) {
					randomVars[i].colIndex = j;
					break;
				}
			}
			nb_rvars_mat++;
		}
	}
	
	if (nb_deleted_rvar > 5) {
		cout << "..." << endl << "Warning: " << nb_deleted_rvar << " random variables are not in the model" << endl;
	}
	
	// in the end, let's get the cumulative probability distribution
	for (int i = 0; i < nb_rvars; i++) {
		randomVars[i].set_cumulative_probs();
	}
	
	cout << "STOCH file has been read successfully." << endl;
}

void SPprob::store_variable_info ()
{
	// record variable indices
	for (int i=0; i<nb_vars; i++) {
		var2index.insert( pair<string, long> (solver[0].var[i].getName(), i) );
	}
	
	// split var vector into stage-many pieces
	assign_vars_to_stages();
	
	lb.resize( vars[0][0].getSize() + vars[0][1].getSize() * nb_scen );
	ub.resize( vars[0][0].getSize() + vars[0][1].getSize() * nb_scen );
	
	integrality.resize( vars[0][0].getSize() + vars[0][1].getSize() * nb_scen );
	
	for (int j=0; j<vars[0][0].getSize(); j++)	{
		lb[j] = vars[0][0][j].getLB();
		ub[j] = vars[0][0][j].getUB();

		if ( vars[0][0][j].getType() != 2 )	integrality[j] = true;		// 	enum Type { Int=1, Float=2, Bool=3 };
	}
	
	long index = -1;
	for (int s=0; s<nb_scen; s++) {
		for (int j=0; j<vars[0][1].getSize(); j++)
		{
			index = vars[0][0].getSize() + s * vars[0][1].getSize() + j;
			
			lb[index] = vars[0][1][j].getLB();
			ub[index] = vars[0][1][j].getUB();

			if ( vars[0][1][j].getType() != 2 )	integrality[index] = true;
		}
	}
}

void SPprob::change_t1var_bounds( vector<double> &new_lbs, vector<double> &new_ubs)
{
	for (int k=0; k<nb_threads; k++)	change_t1var_bounds(new_lbs, new_ubs, k);
}

void SPprob::change_t1var_bounds( vector<double> &new_lbs, vector<double> &new_ubs, int &k )
{
	IloNumArray ilolb (solver[k].env, numberOfVar[0]), iloub (solver[k].env, numberOfVar[0]);
	for (int i=0; i<numberOfVar[0]; i++) {
		ilolb[i] = new_lbs[i];
		iloub[i] = new_ubs[i];
	}
	vars[k][0].setBounds(ilolb, iloub);
	ilolb.end(); iloub.end();
}

void SPprob::change_t2var_bounds( vector<double> &new_lbs, vector<double> &new_ubs, int &s, int k )
{
	IloNumArray ilolb (solver[k].env, numberOfVar[1]), iloub (solver[k].env, numberOfVar[1]);
	
	long st_ind = numberOfVar[0] + numberOfVar[1] * s;
	
	for (int i=0; i<vars[k][1].getSize(); i++) {
		ilolb[i] = new_lbs[st_ind+i];
		iloub[i] = new_ubs[st_ind+i];
	}
	vars[k][1].setBounds(ilolb, iloub);
	ilolb.end(); iloub.end();
}

bool SPprob::solve(int solver_index)
{
	return solver[solver_index].cplex.solve();
}

double SPprob::getObjValue (int solver_index)
{
	return solver[solver_index].cplex.getObjValue();
}

/********************************************************
 Tracks and stores the obj coefs using variable IDs
 Assumption: two-stage problem
 *******************************************************/
void SPprob::get_linear_obj_coefs (int k)	// k: solver index
{
	for (IloExpr::LinearIterator it = IloExpr( solver[k].obj.getExpr() ).getLinearIterator(); it.ok(); ++it)    {
		if ( var2index[(it.getVar()).getName()] < var2index[time_info[1][0]] )
			obj_coefs[k][0][ var2index[(it.getVar()).getName()] ] = it.getCoef();
		else
			obj_coefs[k][1][ var2index[(it.getVar()).getName()] - numberOfVar[0] ] = it.getCoef();
	}
}

void SPprob::getSoln(vector<double> &soln, int stage, int k)
{
	IloNumArray cpx_vals ( solver[k].env );
	
	solver[k].cplex.getValues( cpx_vals , vars[k][stage-1] );
	for (int i=0; i<soln.size(); i++)	soln[i] = cpx_vals[i];
	
	cpx_vals.end();
}

/********************************************************
 Assumption: Quadratic terms cannot be present in the first-stage
 *******************************************************/
void SPprob::use_linear_objective(int k)
{
	solver[k].obj.setLinearCoefs(vars[k][0], obj_coefs[k][0]);
}

/********************************************************
 Assumption: Other quadratic terms cannot be present in the first-stage
 Relax it in the future...
 *******************************************************/
void SPprob::remove_ph_regularizer (int k)
{
	for (int j=0; j<numberOfVar[0]; j++)	solver[k].obj.setQuadCoef(vars[k][0][j], vars[k][0][j], 0.0);
}

/********************************************************
 Assumption: Other quadratic terms cannot be present in the first-stage
 Relax it in the future...
 *******************************************************/
void SPprob::add_ph_regularizer (double &rho, int k)
{
	for (int j=0; j<numberOfVar[0]; j++)	solver[k].obj.setQuadCoef(vars[k][0][j], vars[k][0][j], rho/2.0);
}

void SPprob::use_saddle_objective(vector<double> &lambda, int k)
{
	for (int j=0; j<numberOfVar[0]; j++)	solver[k].obj.setLinearCoef(vars[k][0][j], obj_coefs[k][0][j] + lambda[j]);
}

/********************************************************
 sets the random right-hand-side of the s^th subproblem.
 *******************************************************/
void SPprob::set_random_rhs (int &s, int solver_index)
{
	IloInt c = 0;	// constraint index

	for (int i=0; i<nb_rvars; i++)
	{
		c = randomVars[i].rowIndex;
		
		if (con_type[c] == 'e')	{
			solver[solver_index].rng[c].setBounds(randomVars[i].Value[s], randomVars[i].Value[s]);
		}
		else if (con_type[c] == 'l') {
			solver[solver_index].rng[c].setUB(randomVars[i].Value[s]);
		}
		else {
			solver[solver_index].rng[c].setLB(randomVars[i].Value[s]);
		}
	}
}

double SPprob::setup_regularized_subproblem (int &s, vector<double> &inc, vector<double> &lambda, double &rho, int k)
{	
	// set the right-hand-sides that contain random variables
	set_random_rhs(s, k);
	
	// prepare the objective + regularizer coefficients
	IloNumArray t1_new_obj_coefs (solver[k].env, numberOfVar[0]);
	
	for (int j=0; j<numberOfVar[0]; j++)	t1_new_obj_coefs[j] = obj_coefs[k][0][j] + lambda[j] - rho * inc[j];
	
	solver[k].obj.setLinearCoefs(vars[k][0], t1_new_obj_coefs);

	t1_new_obj_coefs.end();

	// compute the constant terms in the objective function (due to the primal regularization term)
	double constant = 0;
	for (int j=0; j<numberOfVar[0]; j++) constant += pow(inc[j], 2.0);
	constant *= rho/2.0;
	
	return constant;
}

void SPprob::setup_scenario_subproblem(int &s, int solver_index)
{
	set_random_rhs(s, solver_index);
	// other random components (e.g., obj function, technology matrix) should come here
}

void SPprob::convert_to_LP(int k)
{
	solver[k].convertToLP();
}

void SPprob::convert_to_LP()
{
	for (int k=0; k<nb_threads; k++)	solver[k].convertToLP();
}

void SPprob::convert_to_MIP()
{
	for (int k=0; k<nb_threads; k++)	solver[k].convertToMIP(integrality);
}

void SPprob::convert_to_MIP(int k)
{
	solver[k].convertToMIP(integrality);
}

long SPprob::derive_L2_feasibility_cut (vector<double> &soln, vector<int> &fixed_var_indices)
{
	int nb_ones;
	
	for (int k=0; k<nb_threads; k++)
	{
		IloExpr delta ( solver[k].env );
		
		nb_ones = 0;
		for (int i=0; i<fixed_var_indices.size(); i++) {
			if ( soln[ fixed_var_indices[i] ] >= 1-PHBB_IntTol )	{
				delta += vars[k][0][ fixed_var_indices[i] ];
				nb_ones++;
			}
			else	delta -= vars[k][0][ fixed_var_indices[i] ];
		}
		IloRange feas_cut ( solver[k].env, 1, nb_ones - delta);
		delta.end();
		
		node_cuts[k].add(feas_cut);
	}
	
	return node_cuts[0].getSize()-1;
}

void SPprob::add_L2_feasibility_cut(int &index)
{
	for (int k=0; k<nb_threads; k++)	solver[k].model.add( node_cuts[k][index] );
}

void SPprob::remove_L2_feasibility_cut(int &index)
{
	for (int k=0; k<nb_threads; k++)	solver[k].model.remove( node_cuts[k][index] );
}

void SPprob::displayProblemStats()
{
	cout << "Problem Type: " << probType << endl;
	printf("Nb Vars: %ld %s, %ld %s, %ld %s\n", numberOfVar[0], "(T1)", numberOfVar[1], "(T2)", solver[0].var.getSize(), "(Total)");
//	printf("Nb Cons: %ld %s, %ld %s, %ld %s\n", numberOfCon[0], "(T1)", numberOfCon[1], "(T2)", solver[0].rng.getSize(), "(Total)");
	
	cout << "Nb Rand: ";
	if ( nb_rvars_obj > 0 ) cout << nb_rvars_obj << " (obj)  ";
	if ( nb_rvars_rhs > 0 ) cout << nb_rvars_rhs << " (rhs)  ";
	if ( nb_rvars_mat > 0 ) cout << nb_rvars_mat << " (mat)  ";
	
	int maxOfVars = max(nb_rvars_rhs, nb_rvars_obj);
	maxOfVars = max(maxOfVars, nb_rvars_mat);
	if ( maxOfVars != nb_rvars ) cout << nb_rvars_mat << " (Total) ";
	cout << endl;
}
