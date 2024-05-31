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

#define DDES_ALLOW_SYSTEM

#include <iostream>
using namespace std;

#include "capd/capdlib.h"
using namespace capd;
using namespace capd::alglib;
using namespace capd::matrixAlgorithms;

#include "utils.h"

// Main function
int main(int argc, char* argv[])
{
	capd::ddeshelper::ArgumentParser parser(argc, argv);
	int iy = 0, iz = 0;
	int CUT_Y = 1, CUT_Z = 1;
	parser.parse("cuty=", CUT_Y, "how many pieces in y");
	parser.parse("cutz=", CUT_Z, "how many pieces in z");
	parser.parse("iy=", iy, "piece number in y coordinate");
	parser.parse("iz=", iz, "piece number in z coordinate");
	if (parser.isHelpRequested()){
		std::cout << parser.getHelp() << endl;
		return 0;
	}


  	cout.precision(16);
	cout << boolalpha;  
	try
	{			
		interval a = 5.25;
		interval epsi = interval(0.0);
		epsi = interval(0.0001);
		//epsi = interval(0.005);
		interval tau = 0.5;
		system3d roessler525(tau, epsi, a);
		///===================== variables used in Procedure 1:  =====================
		HSet2D grid3 (IVector ({-6.38401, 0.0327544}) , IMatrix ({{-1., 0.000656767}, {-0.000656767, -1.}}) , DVector({3.63687,0.0004}));			// Attractor's container		
		vector<HSet2D> c3(3);
		c3[0] = HSet2D(IVector ({-3.46642, 0.0346316}) , IMatrix ({{-1., 0.000656767}, {-0.000656767, -1.}}) , DVector({0.072,0.00048}));			// cube 1
		c3[1] = HSet2D(IVector ({-6.26401, 0.0326544}) , IMatrix ({{-1., 0.000656767}, {-0.000656767, -1.}}) , DVector({0.162,0.00066}));			// cube 2
		c3[2] = HSet2D(IVector ({-9.74889, 0.0307529}) , IMatrix ({{-1., 0.000656767}, {-0.000656767, -1.}}) , DVector({0.036,0.00072}));			// cube 3

		bool result = roessler525.inside_piece(grid3, grid3, CUT_Y, CUT_Z, iy, iz);
		ostringstream oss; oss << boolalpha;
		oss << iy << "/" << CUT_Y << " " << iz << "/" << CUT_Z << ": " << result;
		cout << oss.str() << endl;
	}
	catch(exception& e)
  	{
		cout << iy << " " << iz << ": Exception caught: "<< e.what() << endl;
  	}
  return 0;
} // main()

