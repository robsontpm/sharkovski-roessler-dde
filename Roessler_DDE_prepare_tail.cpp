/* ========================================================================
 * Roessler DDE system
 * Anna Gierzkiewicz and Robert Szczelina
 * ========================================================================
 * Program: Roessler_Sharkovskii.cpp
 * Date: 18 IV 2024
 * ========================================================================
 * The program contains the procedure to estimate the tails (far and close)
 * for the proof in the accompanied manuscript. It is used to get the
 * initial data for the proofs.
 * ========================================================================*/ 

#include "setup.h"

// Main function
int main(int argc, char* argv[])
{
	// this just reads the arguments to the program
	// you can run program with -h or --help to see the help.
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

		// read the data from files, unfortunately now the names are fixed
		// the names are such that they match names generated with Roessler_DDE_gencoords
		IMatrix C(Mdim, Mdim);
		IMatrix invC(Mdim, Mdim);
		IVector x0(Mdim);
		IVector r0(Mdim);
		IVector Xi(d*p);
		// we read in binary to have exact representaton.
		// This could be changed to reading the data form text files,
		// but then all number would be in machine precision range around read numbers.
		capd::ddeshelper::readBinary(WD + "G_x0.ivector.bin", x0);
		capd::ddeshelper::readBinary(WD + "G_M.imatrix.bin", C);
		capd::ddeshelper::readBinary(WD + "G_invM.imatrix.bin", invC);
		capd::ddeshelper::readBinary(WD + "G_r0.ivector.bin", r0);
		capd::ddeshelper::readBinary(WD + "G_Xi.ivector.bin", Xi);
		// just print some info if 'verbose' or 'debug' flag is presented
		// useful when preparing data of your own proofs.
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

		// compute image of the piece
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
		// we expand old data if necessary
		// (we leave the Xi as the new values, it is important for the first step, when the Xi = [0,0]^{p,d}
		// maybe we need to change that in future. But it works now.
		new_r0 = capd::vectalg::intervalHull(new_r0, r0);
		new_r0[0] = 0.;

		// save data for future use
		string dirpath = WD + "estimate_piece/";
		capd::ddeshelper::mkdir_p(dirpath);

		ostringstream prefix;
		prefix << "_" << iy << "_" << iz;

		// save in human readable form (for debug)
		std::ofstream human_readable_new_r0(dirpath + "human_readable_new_radius" + prefix.str() + ".txt");
		human_readable_new_r0 << "r0:" << endl << new_r0 << endl;
		human_readable_new_r0 << "Xi:" << endl << new_Xi << endl;
		human_readable_new_r0.close();

		// save in exact precision
		capd::ddeshelper::saveBinary(dirpath + "G_new_r0" + prefix.str() + ".ivector.bin", new_r0);
		capd::ddeshelper::saveBinary(dirpath + "G_new_Xi0" + prefix.str() + ".ivector.bin", new_Xi);
	}
	catch(exception& e)
  	{
    	cout << "\n\nException caught: "<< e.what() << endl;
  	}
	std::cout << (result ? "OK" : "KEEP LOOKING") << " at " << iy << " " << iz << endl;
	if (!result) int ignore = system((std::string("touch ") + WD + "KEEP_LOOKING.txt").c_str());
    return 0;
} // main()

