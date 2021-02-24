/**
* Part of the this code is derived from ZMDP: http://www.cs.cmu.edu/~trey/zmdp/
* ZMDP is released under Apache License 2.0
* The rest of the code is released under GPL v2
*/


#include <stdlib.h>
#ifdef _MSC_VER
#else
#include <unistd.h>
#endif

#include <stdio.h>
#include <assert.h>

#include <iostream>
#include <fstream>

#include "MOMDP.h"
#include "SARSOP.h"     // for SARSOPAlphaPlaneTuple
//#include "BlindLBInitializer.h"
//#include "AlphaPlanePool.h"
#include "AlphaPlanePoolSet.h"   // for modifying LB in adaptive management
#include "FastInfUBInitializer.h"
#include "FullObsUBInitializer.h"
#include "BeliefValuePairPoolSet.h"

using namespace std;
//using namespace MatrixUtils;

//#define DEBUGSYL_290908 1
//#define DEBUGSYL_160908a 1

namespace momdp
{

FastInfUBInitializer::FastInfUBInitializer( SharedPointer<MOMDP> problem, BeliefValuePairPoolSet* _bound)
{
	pomdp = ( SharedPointer<MOMDP>) problem;
	bound = _bound;

	// pomdp->XStates->size()
	//alphasByState.resize(pomdp->XStates->size());   // SYL260809 commented out
	//upperBoundBVpair = _upperBoundBVpair;

}

FastInfUBInitializer::FastInfUBInitializer(SharedPointer<MOMDP> problem)
{
	pomdp = ( SharedPointer<MOMDP>) problem;

	// pomdp->XStates->size()
	//alphasByState.resize(pomdp->XStates->size());   // SYL260809 commented out
	//upperBoundBVpair = _upperBoundBVpair;

}

void FastInfUBInitializer::initialize(double targetPrecision)
{
	if(pomdp->XStates->size() != 1 && pomdp->hasPOMDPMatrices())
	{
		// only does this if convert fast is called to produce pomdp version of the matrices
		initFIB_unfac(targetPrecision, false); // need pomdp matrix
	}
	else   // goes there in my case
	{
		initFIB(targetPrecision, false);
	}

}

void FastInfUBInitializer::getFIBsolution(double targetPrecision)
{
	if(pomdp->XStates->size() != 1 && pomdp->hasPOMDPMatrices())
	{
		// only does this if convert fast is called to produce pomdp version of the matrices
		initFIB_unfac(targetPrecision, true); // need pomdp matrix
	}
	else
	{
		initFIB(targetPrecision, true);
	}

	//if (pomdp->XStates->size() == 1)
	//	initFIB(targetPrecision);
	//else
	//	initFIB_unfac(targetPrecision);
}


void FastInfUBInitializer::initMDP_unfac(double targetPrecision)
{
	// set alpha to be the mdp upper bound
	FullObsUBInitializer m;
	m.valueIteration_unfac(pomdp, targetPrecision);
	DenseVector calpha;
	copy(calpha, m.alpha);

	alphas.clear();
	alphas.push_back(calpha);

	// **** ADDED PRINTOUTS TO SCREEN HERE ******** 6 APRIL 2009 *******
	/*  cout << pomdp->YStates->size() << endl;
	FOR (i, pomdp->YStates->size()) {
	cout << calpha(i) << " ";
	}
	cout << endl;
	*/

}

void FastInfUBInitializer::initFIB_unfac(double targetPrecision, bool getFIBvectors)
{
	// calculates the fast informed bound (Hauskrecht, JAIR 2000)
	std::vector< DenseVector > al(pomdp->actions->size());
	std::vector< DenseVector > nextAl(pomdp->actions->size());
	DenseVector tmp, beta_aoi, beta_ao, diff;
	double maxResidual;
	alpha_vector backup;

	initMDP_unfac(MDP_RESIDUAL);

	// TODO:: bound->elapsed = bound->heurTimer.elapsed();
	// TODO:: printf("%.2fs initMDP(MDP_RESIDUAL) done\n", bound->elapsed);

	// initialize al array with weak MDP upper bound
	FullObsUBInitializer m;
	m.pomdp = pomdp;
	alpha_vector& alpha = alphas[0];
	copy(m.alpha, alpha);

	FOR (a, pomdp->actions->size())
	{
		//al[a].resize(pomdp->YStates->size());
		//nextAl[a].resize(pomdp->YStates->size());
		al[a].resize(pomdp->XStates->size() * pomdp->YStates->size());
		nextAl[a].resize(pomdp->XStates->size() * pomdp->YStates->size());
		m.nextAlphaAction_unfac(al[a], a);
	}

	// iterate FIB update rule to approximate convergence
	do
	{
		FOR (a, pomdp->actions->size())
		{
			DenseVector& beta_a = nextAl[a];

			set_to_zero( beta_a );


			FOR (o, pomdp->observations->size())
			{
				FOR (i, pomdp->actions->size())
				{
					emult_column( tmp, *(*(pomdp->pomdpO))[a], o, al[i] );
					//	  emult_column( tmp, pomdp->pomdpO[a], o, al[i] );
					mult( beta_aoi, tmp, *(*(pomdp->pomdpTtr))[a] );
					//	  mult( beta_aoi, tmp, pomdp->pomdpTtr[a] );
					if (0 == i)
					{
						beta_ao = beta_aoi;
					}
					else
					{
						max_assign( beta_ao, beta_aoi );
					}
				}
				beta_a += beta_ao;
			}

			beta_a *= pomdp->discount;
			copy_from_column( tmp, (*(pomdp->pomdpR)), a );
			//     copy_from_column( tmp, pomdp->pomdpR, a );
			beta_a += tmp;
		}

		maxResidual = 0;
		FOR (a, pomdp->actions->size())
		{
			diff = nextAl[a];
			diff -= al[a];
			maxResidual = std::max( maxResidual, diff.norm_inf() );

			al[a] = nextAl[a];
		}

	}
	while ( maxResidual > targetPrecision );

	if (!getFIBvectors)
	{
		DenseVector dalpha;
		FOR (a, pomdp->actions->size())
		{
			if (0 == a)
			{
				dalpha = al[a];
			}
			else
			{
				max_assign(dalpha, al[a]);
			}
		}

		// post-process: make sure the value for all terminal states
		// is exactly 0, since that is how the ubVal field of terminal
		// nodes is initialized.
		FOR (i, pomdp->XStates->size()*pomdp->YStates->size())
		{

			// convert i to x and y values
			unsigned int x = (unsigned int) i/(pomdp->YStates->size());
			unsigned int y = i % (pomdp->YStates->size());

			if (pomdp->isPOMDPTerminalState[x][y])
			{

				dalpha(i) = 0.0;
			}
		}

		// **** ADDED PRINTOUTS TO SCREEN HERE ******** 6 APRIL 2009 *******
		/*  cout << pomdp->YStates->size() << endl;
		FOR (i, pomdp->YStates->size()) {
		cout << dalpha(i) << " ";
		}
		cout << endl;
		*/
		// write out result
		//  bound->points.clear();
		//  copy(upperBoundBVpair->cornerPoints, dalpha);

		// write dalpha entries into dalphas[state_idx] entries
		std::vector<DenseVector> dalphas(pomdp->XStates->size());

		FOR (s, pomdp->XStates->size() * pomdp->YStates->size())
		{
			// convert i to x and y values
			unsigned int x = (unsigned int) s/(pomdp->YStates->size());
			unsigned int y = s % (pomdp->YStates->size());

			if (y==0) dalphas[x].resize(pomdp->YStates->size());

			dalphas[x](y) = dalpha(s);
		}

		// write out result - do it for each of the BoundsSet
		FOR (state_idx, pomdp->XStates->size())
		{
			bound->set[state_idx]->points.clear();
			copy(bound->set[state_idx]->cornerPoints,dalphas[state_idx]);

			/*cout << "state_idx  " << state_idx << endl;
			dalphas[state_idx].write(std::cout);
			cout << endl; */

		}
	}
	else
	{
		actionAlphaByState.resize(pomdp->actions->size());

		FOR (a, pomdp->actions->size())
		{
			FOR (i, pomdp->XStates->size()*pomdp->YStates->size())
			{

				// convert i to x and y values
				unsigned int x = (unsigned int) i/(pomdp->YStates->size());
				unsigned int y = i % (pomdp->YStates->size());

				if (pomdp->isPOMDPTerminalState[x][y])
				{

					al[a](i) = 0.0;
				}

				if (x==0) actionAlphaByState[a].resize(pomdp->XStates->size());
				if (y==0) actionAlphaByState[a][x].resize(pomdp->YStates->size());
				actionAlphaByState[a][x](y) = al[a](i);
			}
		}

	}
}

void FastInfUBInitializer::initMDP(double targetPrecision, FullObsUBInitializer& m) //called
// SYL260809 prevly:	void FastInfUBInitializer::initMDP(double targetPrecision)
{
	// set alpha to be the mdp upper bound
	// FullObsUBInitializer m;		// SYL260809 commented out
	m.valueIteration(pomdp, targetPrecision); // this puts the MDP vector at m.alphaByState

	//alphasByState.clear();

	//alphasByState = m.alphaByState; // SYL260809 commented out

	/* #if DEBUGSYL_160908
	//### check alphasByState contents compared to m.alphaByState contents

	cout << "m.alphaByState" << endl;
	for (unsigned int stateidx=0; stateidx < m.alphaByState.size(); stateidx++)
	{
	m.alphaByState[stateidx].write(std::cout);
	cout << endl;
	}
	cout << "alphasByState" << endl;
	for (unsigned int stateidx=0; stateidx < alphasByState.size(); stateidx++)
	{
	alphasByState[stateidx].write(std::cout);
	cout << endl;
	}
	#endif */



	//std::vector<DenseVector> calphas;
	//copy(calpha, m.alpha);


	//alphas.clear();
	//alphas.push_back(calpha);

	/*
	if (zmdpDebugLevelG >= 1) {
	cout << "initUpperBoundMDP: alpha=" << sparseRep(alphas[0]).c_str() << endl;
	cout << "initUpperBoundMDP: val(b)=" << inner_prod(alphas[0], pomdp->initialBelief) << endl;
	}
	*/
}

void FastInfUBInitializer::initFIB(double targetPrecision, bool getFIBvectors)   //called
{
	// calculates the fast informed bound (Hauskrecht, JAIR 2000)
	//std::vector< DenseVector > al(pomdp->actions->size());
	//std::vector< DenseVector > nextAl(pomdp->actions->size());

	std::vector< std::vector<DenseVector> > al(pomdp->actions->size());
	std::vector< std::vector<DenseVector> > nextAl(pomdp->actions->size());

	DenseVector tmp, tmp2, diff, beta_ao; //beta_aoi, beta_aoiSum
	DenseVector beta_aoXn_unweighted, beta_aoXn;
	double maxResidual;
	//alpha_vector backup;

	FullObsUBInitializer M;  // SYL260809 added

	initMDP(MDP_RESIDUAL, M); // SYL260809 modified

	// SYL260809 commented out - the MDP upper bound is already in M.alphaByState

	FOR (a, pomdp->actions->size())
	{
		al[a].resize(pomdp->XStates->size());
		nextAl[a].resize(pomdp->XStates->size());

		M.nextAlphaAction(al[a], a); // SYL260809 prevly: m.nextAlphaAction(al[a], a);
	}

	if (pomdp->adaptiveManagement){ //evaluate each corner's optimal policy on all corners
		// to generate good LB alpha planes.
		vector< vector< int > > bestActions(pomdp->YStates->size()); // best action in each X and Y.

		M.findBestActions(bestActions);
		FOR (Yc, pomdp->YStates->size()){
			vector< alpha_vector > alpha(pomdp->XStates->size());
			if (Yc == 0) {
				FOR (Xc, pomdp->XStates->size()){
					cout << "best action for state " << Xc << " : " << bestActions[Yc][Xc] <<endl;
				}

               }
			M.evaluateCornerPolicies(pomdp, targetPrecision, bestActions[Yc], alpha);




			FOR (Xc, pomdp->XStates->size()){
				SharedPointer<AlphaPlane> plane (new AlphaPlane());
//			cornerAlphaPlanes[Xc].resize(pomdp->YStates->size());

//				cornerAlphaPlanes[Xc][Yc] = plane;

				// Create plane and store in lower bound through member 'lbSet'
				copy(*plane->alpha, alpha[Xc]);
				plane->action = bestActions[Yc][Xc];
				plane->sval = Xc;
				plane->setTimeStamp(0);


				SARSOPAlphaPlaneTuple *tempTuple = (SARSOPAlphaPlaneTuple *)plane->solverData;
				tempTuple->certed = 0; //init certed count to 0

				if (Xc == pomdp->uniqueInitialX){
//				if (Xc == pomdp->XStates->size() - 2){
//				if (Xc == 3){
//				cout << endl << endl<< " X= " << Xc << " Y= " << Yc << endl ;
//					FOR (Ycc, pomdp->YStates->size())
//						cout<< " Y= " << Ycc<< " - value =  " << alpha[Xc](Ycc);
				}



				bound->lbSet->set[Xc]->addAlphaPlane(plane);

			}
		}
//		int yy; cin >> yy;

	}


	bool	valid_o_and_Xn; // SYL260809 to allow skipping operation for an observation o, if there are
	// no valid Xn values (for given action a) associated with it
	int steps = 0;
	do
	{
		steps++;
		// SYL260809 keep track of maximum residual across actions
		maxResidual = 0;

		// the following loops illustrate well the complexity of one FIB upate: A^2.S^2.O.
		FOR (a, pomdp->actions->size())
		{
			std::vector< DenseVector >& beta_a = nextAl[a];

			FOR (Xc, pomdp->XStates->size())
			{
				beta_a[Xc].resize(pomdp->YStates->size());

				FOR (o, pomdp->observations->size())
				{
					beta_ao.resize(pomdp->YStates->size());
					valid_o_and_Xn = false; // SYL260809

					// only iterate over possible X states
					const vector<int>& possibleXns = pomdp->XTrans->getMatrix(a, Xc)->nonEmptyColumns();
					FOREACH (int, XnIt, possibleXns)
					{
						int Xn = *XnIt;
						if (!(pomdp->obsProb->getMatrix(a, Xn)->isColumnEmpty(o)))
						{
							valid_o_and_Xn = true;	 // SYL260809

							FOR (i, pomdp->actions->size())
							{
								emult_column( tmp, *pomdp->obsProb->getMatrix(a, Xn), o, al[i][Xn] );
								mult( tmp2, *pomdp->YTrans->getMatrix(a, Xc, Xn), tmp); // SYL270809
								// SYL270809 mult(vector, matrix, vector) is faster than mult(vector, vector, matrix)
								//mult( tmp2, tmp, *pomdp->XYTrans->getMatrixTr(a, Xc) );

								if (0 == i)
								{
									beta_aoXn_unweighted = tmp2;
								}
								else
								{
									max_assign( beta_aoXn_unweighted, tmp2 );
								}

							}

							emult_column(beta_aoXn, *pomdp->XTrans->getMatrix(a, Xc), Xn, beta_aoXn_unweighted);
							beta_ao += beta_aoXn;
						}
					}
					// SYL260809
					if (valid_o_and_Xn)	// false means there's no valid Xn (given the
						//action and Xc) for this observation o
						beta_a[Xc] += beta_ao; //beta_ao;

				}

				beta_a[Xc] *= pomdp->discount;

				copy_from_column( tmp, *pomdp->rewards->getMatrix(Xc), a );
				beta_a[Xc] += tmp;

				// SYL260809 keep track of maximum residual
				diff = beta_a[Xc];
				diff -= al[a][Xc];
				maxResidual = std::max( maxResidual, diff.norm_inf() );
			}
		}

		// SYL260809  assign the FIB vectors for next iteration
		al = nextAl;

		/*
		if (zmdpDebugLevelG >= 1) {
		cout << ".";
		cout.flush();
		}
		*/
	}
	while ( maxResidual > targetPrecision );

//	cout <<  endl << "steps: " << steps;
//	cout << endl << targetPrecision<< "  / "  << diff.norm_inf() << "  / " << maxResidual << endl;
//	int pp ; cin >> pp;


//

	// Prepare to return the best alpha-vectors for each state:
//	cout << endl << "Display alpha vectors: " << endl;
//	FOR (s, pomdp->getBeliefSize())
//	{
//		FOR (Xc, al[0].size())
//		if (0!= pomdp->initialBeliefX->data[Xc] )
//		{
//			cout << endl << endl << "x: "  << Xc << ", y: " << s  ;
//			double maxAction = 0;
//			int bestAction = 0;
//			FOR (a, pomdp->actions->size())
//			{
//				if (al[a][Xc](s) > maxAction)
//				{
//					maxAction = al[a][Xc](s);
//					bestAction = a;
//				}
//			}
//			cout << endl << "Best action:  " << bestAction << ". Values: ";
//			FOR (ss, pomdp->getBeliefSize())
//			{
//				cout << al[bestAction][Xc](ss) << " / ";
//			}
//		}
//	}
//
//	int qq ;
//	cin >> qq;

	//cout << "targetPrecision : " << targetPrecision << endl;


	/*
	if (zmdpDebugLevelG >= 1) {
	cout << endl;
	}
	*/
#if DEBUGSYL_290908
	cout << "targetPrecision : " << targetPrecision << endl;
	cout << "maxResidual : " << maxResidual << " iterationCount : " << iterationCount << endl;
#endif

#if DEBUGSYL_160908
	cout << "targetPrecision : " << targetPrecision << endl;
	cout << "maxResidual : " << maxResidual << " iterationCount : " << iterationCount << endl;
	cout << "After iteration, al " << endl;
	FOR (a, pomdp->actions->size())
	{
		cout << "a : " << a << endl;

		for (unsigned int stateidx=0; stateidx < alphasByState.size(); stateidx++)
		{
			cout << "stateidx : " << stateidx << endl;
			al[a][stateidx].write(std::cout);
			cout << endl;
		}
	}
#endif

	if (!getFIBvectors) // where the upper bound is created. Up to now they is alpha-vector per
		// action, state x and state y. Now only the upper bound over y is kept (max_assign).
	{
		// SYL260809 we only need one dalpha of size YStates.size at a time
		DenseVector dalpha;
		//std::vector<DenseVector> dalpha(pomdp->XStates->size());
		FOR (state_idx, pomdp->XStates->size())
		{

			FOR (a, pomdp->actions->size())
			{
				if (0 == a)
				{
					dalpha = al[a][state_idx];
				}
				else
				{
					max_assign(dalpha, al[a][state_idx]);
				}
			}

			// at this point, the vector at dalpha[state_idx] contains the highest value, across actions

			// post-process: make sure the value for all terminal states
			// is exactly 0, since that is how the ubVal field of terminal
			// nodes is initialized.
			FOR (i, pomdp->YStates->size())
			{
				if (pomdp->isPOMDPTerminalState[state_idx][i])
				{
					dalpha(i) = 0.0;
				}
			}

			// at this point, the vector at dalpha[state_idx] has taken into account terminal state
			bound->set[state_idx]->points.clear();
			copy(bound->set[state_idx]->cornerPoints,dalpha);

		}

	}
	else      // output one vector for each action
	{

		FOR (state_idx, pomdp->XStates->size())
		{

			FOR (a, pomdp->actions->size())
			{
				// post-process: make sure the value for all terminal states
				// is exactly 0
				FOR (i, pomdp->YStates->size())
				{
					if (pomdp->isPOMDPTerminalState[state_idx][i])
					{
						al[a][state_idx](i) = 0.0;
					}
				}
			}

		}

		actionAlphaByState = al;

	}
}

}; // namespace zmdp


