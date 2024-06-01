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
	std::string outname = "test";
	std::string str_src = "grid3", str_dst = "grid3";
	parser.parse("cuty=", CUT_Y, "how many pieces in y");
	parser.parse("cutz=", CUT_Z, "how many pieces in z");
	parser.parse("iy=", iy, "piece number in y coordinate");
	parser.parse("iz=", iz, "piece number in z coordinate");
	parser.parse("src=", str_src, "id of the lhs set in covering relation: 'c3[0]', 'c3[1]', 'c3[2]', or 'grid3'");
	parser.parse("dst=", str_dst, "id of the rhs set in covering relation: 'c3[0]', 'c3[1]', 'c3[2]', or 'grid3'");
	parser.parse("out=", outname, "prefix of the name of the files generated with this program");
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

		auto choice_set = [&](std::string const& name) -> HSet2D& {
			if (name == "c3[0]") return c3[0];
			if (name == "c3[1]") return c3[1];
			if (name == "c3[2]") return c3[2];
			return grid3;
		};

		auto choice_mid = [&](std::string const& name) {
			string filename = "grid3_x0.ivector.bin";
			if (name == "c3[0]") filename = "c3_0_x0.ivector.bin";
			if (name == "c3[1]") filename = "c3_1_x0.ivector.bin";
			if (name == "c3[2]") filename = "c3_2_x0.ivector.bin";
			IVector mid_point(roessler525.M());
			capd::ddeshelper::readBinary(filename, mid_point);
			return mid_point;
		};
		IVector rel_r0(roessler525.M());
		IVector rel_Xi(roessler525.d * roessler525.p);
		capd::ddeshelper::readBinary("grid3_new_r0.ivector.bin", rel_r0);
		capd::ddeshelper::readBinary("grid3_new_Xi.ivector.bin", rel_Xi);

		HSet2D& src_set = choice_set(str_src);
		auto src_mid = choice_mid(str_src);
		HSet2D& dst_set = choice_set(str_dst);
		auto dst_mid = choice_mid(str_dst);

		IVector Phull(roessler525.M()), PXi(roessler525.p * roessler525.d);
		bool result = roessler525.inside_piece(src_set, src_mid, dst_set, dst_mid, CUT_Y, CUT_Z, iy, iz, Phull, PXi);
		ostringstream oss; oss << boolalpha;
		oss << str_src << "=>" << str_dst << " " << iy << "/" << CUT_Y << " " << iz << "/" << CUT_Z << ": " << result;
		cout << oss.str() << endl;

		ostringstream prefix;
		prefix << outname << "-iy" << iy << "-" << iz;
		string dirpath = "Pimages/";
		capd::ddeshelper::mkdir_p(dirpath);
		ofstream out(dirpath + prefix.str() + ".ivector");
		out.precision(16);
		out << Phull << endl;
		out.close();

		// compute the size relative to grid3 on the tail and output data suitable for drawing!
		for (int i = 3; i < roessler525.M(); i++)
			Phull[i] /= rel_r0[i].rightBound();

		interval hullPx = Phull[3];
		for (int i = 4; i < Phull.dimension(); i++)
			hullPx = intervalHull(hullPx, (interval)Phull[i]);

		IVector rel_Xi_rad = rel_Xi;
		capd::vectalg::split(rel_Xi, rel_Xi_rad);
		PXi -= rel_Xi;
		for (int i = 0; i < PXi.dimension(); i++)
			PXi[i] /= rel_Xi_rad[i].rightBound();

		interval hullXi = PXi[0];
		for (int i = 1; i < PXi.dimension(); i++)
			hullXi = intervalHull(hullXi, (interval)PXi[i]);

		ofstream dat(dirpath + prefix.str() + ".dat");
		dat.precision(16);

		dat << capd::ddeshelper::to_dat(Phull[1]) << " ";
		dat << capd::ddeshelper::to_dat(Phull[2]) << " ";
		dat << capd::ddeshelper::to_dat(hullPx) << " ";
		dat << capd::ddeshelper::to_dat(hullXi) << " ";

		dat.close();
	} catch(exception& e){

		cout << iy << " " << iz << ": Exception caught: "<< e.what() << endl;
  	}
  return 0;
} // main()

