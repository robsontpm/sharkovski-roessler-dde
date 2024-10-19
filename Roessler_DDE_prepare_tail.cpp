/* ========================================================================
 * Roessler DDE system
 * Anna Gierzkiewicz and Robert Szczelina
 * ========================================================================
 * Program: Roessler_Sharkovskii.cpp
 * Date: 18 IV 2024
 * ========================================================================
 * The program contains the procedure to estimate the tails (far and close)
 * for the proof in the accompanied manuscript.
 * ========================================================================*/ 
//
//#define DDES_ALLOW_SYSTEM

//#include <iostream>
//using namespace std;
//
//#include "capd/capdlib.h"
//using namespace capd;
//using namespace capd::alglib;
//using namespace capd::matrixAlgorithms;

#include "setup.h"

// Main function
int main(int argc, char* argv[])
{
	capd::ddeshelper::ArgumentParser parser(argc, argv);
	int iy = 0, iz = 0;
	int CUTS_Y = 1, CUTS_Z = 1;
	string WD = "./tmp/";
	bool verbose = false;
	interval EPSI = 0.;
	parser.parse("epsi=", EPSI, "epsi value for the computation.");
	parser.parse("cuty=", CUTS_Y, "how many pieces in y");
	parser.parse("cutz=", CUTS_Z, "how many pieces in z");
	parser.parse("iy=", iy, "piece number in y coordinate");
	parser.parse("iz=", iz, "piece number in z coordinate");
	parser.parse("wd=", WD, "the working directory of the program, i.e. where to look for data and to put results. Should end with a '/'.");
	verbose = parser.parse("verbose", std::string("print more output")) || parser.parse("debug", std::string("same as verbose"));
	if (parser.isHelpRequested()){
		std::cout << parser.getHelp() << endl;
		return 0;
	}

	// sanitize input
	if (WD == "") WD = "./";
	if (WD.back() != '/') WD += "/";
	capd::ddeshelper::mkdir_p(WD);

  	cout.precision(16);	// precision of nonrigorous output
	cout << boolalpha;  // outputting bools as true/false
	SHA_DEBUG("# " << parser.getCommandLine()) << endl;
	SHA_DEBUG("# The working Directory is '" << WD << "'") << endl;
	bool result = false;
	try
	{			
		// get the configuration of the Roessler system with the selected parameters.
		Config local_config(EPSI);
		// the names are inherited from paper [GZ2022]
		auto& roessler525 = config.roessler525;
		auto& grid3 = config.grid3;
		auto Mdim = roessler525.M();
		int p = roessler525.p;
		int d = roessler525.d;

		IMatrix C(Mdim, Mdim);
		IMatrix invC(Mdim, Mdim);
		IVector x0(Mdim);
		IVector r0(Mdim);
		IVector Xi(d*p);
		capd::ddeshelper::readBinary(WD + "G_x0.ivector.bin", x0);
		capd::ddeshelper::readBinary(WD + "G_M.imatrix.bin", C);
		capd::ddeshelper::readBinary(WD + "G_invM.imatrix.bin", invC);
		capd::ddeshelper::readBinary(WD + "G_r0.ivector.bin", r0);
		capd::ddeshelper::readBinary(WD + "G_Xi.ivector.bin", Xi);
		SHA_DEBUG("x0:      " << ddeshelper::slice(x0)) << endl;
		SHA_DEBUG("M:       ") << endl;
		SHA_DEBUG("         " << ddeshelper::slice(C.row(0))) << endl;
		SHA_DEBUG("         " << ddeshelper::slice(C.row(1))) << endl;
		SHA_DEBUG("         " << ddeshelper::slice(C.row(2))) << endl;
		SHA_DEBUG("M^{-1}:  ") << endl;
		SHA_DEBUG("         " << ddeshelper::slice(invC.row(0))) << endl;
		SHA_DEBUG("         " << ddeshelper::slice(invC.row(1))) << endl;
		SHA_DEBUG("         " << ddeshelper::slice(invC.row(2))) << endl;
		SHA_DEBUG("r0:      " << ddeshelper::slice(r0)) << endl;
		SHA_DEBUG("Xi:      " << ddeshelper::slice(Xi)) << endl;

		// compute images
		IVector new_r0 = r0;
		IVector new_Xi(p * d);
		result = roessler525.inside_piece(
			grid3, x0,
			grid3, x0,
			r0, Xi,
			C, invC,
			CUTS_Y, CUTS_Z, iy, iz,
			new_r0, new_Xi
		);
		new_r0 = capd::vectalg::intervalHull(new_r0, r0);
		new_r0[0] = 0.;

		string dirpath = WD + "estimate_piece/";
		capd::ddeshelper::mkdir_p(dirpath);

		ostringstream prefix;
		prefix << "_" << iy << "_" << iz;

		std::ofstream human_readable_new_r0(dirpath + "human_readable_new_radius" + prefix.str() + ".txt");
		human_readable_new_r0 << "r0:" << endl << new_r0 << endl;
		human_readable_new_r0 << "Xi:" << endl << new_Xi << endl;
		human_readable_new_r0.close();

		capd::ddeshelper::saveBinary(dirpath + "G_new_r0" + prefix.str() + ".ivector.bin", new_r0);
		capd::ddeshelper::saveBinary(dirpath + "G_new_Xi0" + prefix.str() + ".ivector.bin", new_Xi);
	}
	catch(exception& e)
  	{
    	cout << "\n\nException caught: "<< e.what() << endl;
  	}
	std::cout << (result ? "OK" : "KEEP LOOKING") << " at " << iy << " " << iz << endl;
	if (!result) int ignore = system((std::string("touch ") + WD + "KEEP_LOOKING.txt").c_str());
    return result;
} // main()

