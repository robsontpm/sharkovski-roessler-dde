/* ========================================================================
 * Roessler DDE system
 * Anna Gierzkiewicz and Robert Szczelina
 * ========================================================================
 * Program: Roessler_Sharkovskii.cpp
 * Date: 18 IV 2024
 * ========================================================================
 * The program contains the proofs of existence of contracting grids
 * in the Roessler system ...
 * ========================================================================*/ 

#include <iostream>
using namespace std;

#include "capd/capdlib.h"
using namespace capd;
using namespace capd::alglib;
using namespace capd::matrixAlgorithms;

#include "utils.h"

// Main function

int main()
{
  	cout.precision(16);
	cout << boolalpha;  
	try
	{			
		system3d roessler525(interval(5.25));						// The Roessler system with a=5.25 defined 					
		///===================== variables used in Procedure 1:  =====================
		HSet2D grid3 (IVector ({-6.38401, 0.0327544}) , IMatrix ({{-1., 0.000656767}, {-0.000656767, -1.}}) , DVector({3.63687,0.0004}));			// Attractor's container		
		vector<HSet2D> c3(3);
		c3[0] = HSet2D(IVector ({-3.46642, 0.0346316}) , IMatrix ({{-1., 0.000656767}, {-0.000656767, -1.}}) , DVector({0.072,0.00048}));			// cube 1
		c3[1] = HSet2D(IVector ({-6.26401, 0.0326544}) , IMatrix ({{-1., 0.000656767}, {-0.000656767, -1.}}) , DVector({0.162,0.00066}));			// cube 2
		c3[2] = HSet2D(IVector ({-9.74889, 0.0307529}) , IMatrix ({{-1., 0.000656767}, {-0.000656767, -1.}}) , DVector({0.036,0.00072}));			// cube 3
		
		roessler525.makeHistory(c3[0]);

		const int CUTS_Y = 1000;
		const int CUTS_Z = 3;
		GridSet gridSet(2);
		grid3.gridSet(gridSet,CUTS_Y,CUTS_Z);
		for(auto i = gridSet.begin();i!=gridSet.end();++i)
		{
			Set3d = C0Rect2Set(expand(*i), );
			roessler525.makeHistoryD(grid3, CUTS_Y, CUTS_Z);
		}
		return 0;

		cout << "===========================================================" << endl;
		cout << "|| Roessler system, a = 5.25" << endl;
		cout << "===========================================================" << endl;
	
		cout << "P(C1)<C2? ... " << roessler525.inside(c3[0],c3[1],3,1) << endl;		// true
		cout << "----------------------------------------" <<  endl;
		cout << "P(C2)<C3? ... " << roessler525.inside(c3[1],c3[2],20,1) << endl;		// true
		cout << "----------------------------------------" << endl;
		cout << "P(C3)<C1? ... " << roessler525.inside(c3[2],c3[0],40,1) << endl;		// true
		cout << "----------------------------------------" <<  endl;
		cout << "Is the grid G3 forward-invariant? ... " << roessler525.inside(grid3,grid3,CUTS_Y,CUTS_Z) << endl;		// check if P(grid) < grid divided into 500x3 pieces
		cout << "----------------------------------------" << endl;				//~ {	
	}
	catch(exception& e)
  	{
    	cout << "\n\nException caught: "<< e.what() << endl;
  	}
  return 0;
} 
