/* ========================================================================
 * Roessler DDE system
 * Anna Gierzkiewicz and Robert Szczelina
 * ========================================================================
 * Program: Roessler_Sharkovskii.cpp
 * Date: 18 IV 2024
 * ========================================================================
 * The program contains the proofs of existence of contracting grids
 * in the Roessler system. It works on chunks of pieces from the original
 * proof in finite dimensions from [GZ2022] (see utils.h for reference)
 *
 * The program is scripted for the setup of 3-periodic orbit. It should be
 * easy to change it to work also for other sets of coverings from [GZ2022]
 * ========================================================================*/ 

#include "setup.h"

// Main function
int main(int argc, char* argv[])
{
	// this is a helper class to parse arguments from command line
	// we will prepare a lot of commands and run them in parallel to get the
	// speed boost. Each computation will check the conditions of the small
	// chunk of the original set. The arguments from comand line are nice for
	// controlling that.
	// you can run program with -h or --help to see the help.
	capd::ddeshelper::ArgumentParser parser(argc, argv);
	int iy = 0, iz = 0;
	int CUT_Y = 1, CUT_Z = 1;
	std::string outname = "test";
	std::string str_src = "G", str_dst = "G";
	string WD = "./";
	bool verbose = false;
	parser.parse("cuty=", CUT_Y, "how many pieces in y");
	parser.parse("cutz=", CUT_Z, "how many pieces in z");
	parser.parse("iy=", iy, "piece number in y coordinate");
	parser.parse("iz=", iz, "piece number in z coordinate");
	parser.parse("src=", str_src, "id of the lhs set in covering relation: 'C[i]' - i-th set, or 'G' - the grid");
	parser.parse("dst=", str_dst, "id of the rhs set in covering relation: 'C[i]' - i-th set, or 'G' - the grid");
	parser.parse("out=", outname, "prefix of the name of the files generated with this program");
	parser.parse("wd=", WD, "the working directory of the program, i.e. where to look for data and to put results. Should end with a '/'.");
	verbose = parser.parse("verbose", std::string("print more output")) || parser.parse("debug", std::string("same as verbose"));
	if (parser.isHelpRequested()){
		std::cout << parser.getHelp() << endl;
		return 0;
	}

	// sanitize input
	if (WD == "") WD = "./";
	if (WD.back() != '/') WD += "/";

  	cout.precision(16);	// precision of nonrigorous output
	cout << boolalpha;  // outputting bools as true/false
	SHA_DEBUG("# " << parser.getCommandLine()) << endl;
	SHA_DEBUG("# The working Directory is '" << WD << "'") << endl;
	try
	{			
		// for shorter writing later
		auto& roessler525 = config.roessler525;
		auto& C = config.c3;
		auto& G = config.grid3;

		// this is a helper function to decide which C_i set to choose.
		// we use lambda functions to easily link with external variables, like C, G, etc.
		auto get_C_number = [&](std::string const& name) -> int {
			if (name[0] != 'C') throw std::logic_error("Bad set name - should star either with C or G.");
			if (name[1] != '_') throw std::logic_error("Bad set name: no underscore in C_i!");
			std::istringstream iss(name.substr(2));
			if (name[2] < '0' || '9' < name[2]) throw std::logic_error("Bad set name: should contain number.");
			int i = 0; iss >> i;
			if (i < 0 || i > C.size()) throw std::logic_error("Bad set name: selected set outside of the range.");
			return i;
		};
		// this is a helper function to select right set.
		// it also checks if the data given by the user is well-formed
		// this is just technical, and has nothing to do with the actual proof.
		auto choice_set = [&](std::string const& name) -> HSet2D& {
			if (name[0] == 'G') {
				SHA_DEBUG("using set G") << endl;
				return G;
			}
			int i = get_C_number(name);
			SHA_DEBUG("using set C_" << i) << endl;
			return C[i];
		};
		// similar to above, but this one returns mid point of the set
		// to be used. Mid point need to be in the full space of the DDE
		// so we cannot get it only from ODE finite dimensional data.
		auto choice_mid = [&](std::string const& name) {
			string filename = "G_x0.ivector.bin";
			if (name[0] != 'G'){
				int i = get_C_number(name);
				ostringstream oss; oss << "C_" << i << "_x0.ivector.bin";
				filename = oss.str();
			}
			IVector mid_point(roessler525.M());
			capd::ddeshelper::readBinary(WD + filename, mid_point);
			SHA_DEBUG("using midpoint from '" << filename << "'") << endl;
			return mid_point;
		};

		// the r0 and Xi (tail) part would be the same for all sets
		// except maybe at the head z(x_0). We will set the head from the
		// selected finite dimensional sets for ODE
		IVector rel_r0(roessler525.M());
		IVector rel_Xi(roessler525.d * roessler525.p);
		capd::ddeshelper::readBinary(WD + "G_r0.ivector.bin", rel_r0);
		capd::ddeshelper::readBinary(WD + "G_Xi.ivector.bin", rel_Xi);
		IMatrix M(roessler525.M(), roessler525.M()), invM(roessler525.M(), roessler525.M());
		capd::ddeshelper::readBinary(WD + "G_M.imatrix.bin", M);
		capd::ddeshelper::readBinary(WD + "G_invM.imatrix.bin", invM);

		// now we select input data using helper functions
		HSet2D& src_set = choice_set(str_src);
		auto src_mid = choice_mid(str_src);
		HSet2D& dst_set = choice_set(str_dst);
		auto dst_mid = choice_mid(str_dst);

		SHA_DEBUG("src_mid: " << ddeshelper::slice(src_mid)) << endl;
		SHA_DEBUG("dst_mid: " << ddeshelper::slice(dst_mid)) << endl;
		SHA_DEBUG("r0:      " << ddeshelper::slice(rel_r0)) << endl;
		SHA_DEBUG("Xi:      " << ddeshelper::slice(rel_Xi)) << endl;

		// in those two variables we store the image of the selected piece after the poincare map
		IVector Phull(roessler525.M()), PXi(roessler525.p * roessler525.d);
		// this is where we check the inclusion! Check details there.
		bool result = roessler525.inside_piece(
			src_set, src_mid,		// left side of the covering relation
			dst_set, dst_mid,		// right side of the covering relation
			rel_r0, rel_Xi,			// the common data (the tail)
			M, invM, 				// the coordinates
			CUT_Y, CUT_Z, iy, iz,	// the definition of the piece to be computed
			Phull, PXi				// the image of the piece, transformed to good coordinates
		);

		// output information on the check, if true, we have proven the inclusion
		cout << str_src << "=>" << str_dst << " " << iy << "/" << CUT_Y << " " << iz << "/" << CUT_Z << ": " << result << endl;

		// what follows is just for pictures, and is not important from the proof perspective.
		// save data for making pictures
		ostringstream prefix;
		prefix << "iy-" << iy << "--iz-" << iz;
		string dirpath = outname + "-Pimages/";
		capd::ddeshelper::mkdir_p(dirpath);
		auto filepath = dirpath + prefix.str() + ".ivector";
		ofstream out(filepath);
		out.precision(16);
		out << Phull << endl;
		out.close();
		SHA_DEBUG("P(set) saved in '" << filepath << "'") << endl;

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

		filepath = dirpath + prefix.str() + ".dat";
		ofstream dat(filepath);
		dat.precision(16);

		// just to have info what we have in the file
		// the y coordinate is more like in the dominant direction, computed during
		// the preparation procedure (see paper and/or docs in Roesler_DDE_gencoords for more explanation)
		dat << "# y_mid y_radius z_mid z_radius tail_hull_mid tail_hull_radius Xi_hull_mid Xi_hull_radius\n";
		dat << capd::ddeshelper::to_dat(Phull[1]) << " ";
		dat << capd::ddeshelper::to_dat(Phull[2]) << " ";
		dat << capd::ddeshelper::to_dat(hullPx) << " ";
		dat << capd::ddeshelper::to_dat(hullXi) << " ";

		dat.close();
		SHA_DEBUG("Max radii of head, tail and Xi saved in '" << filepath << "'") << endl;
	} catch(exception& e){
		// output information about exception and the piece indices
		cout << iy << " " << iz << ": Exception caught: "<< e.what() << endl;
		cout << str_src << "=>" << str_dst << " " << iy << "/" << CUT_Y << " " << iz << "/" << CUT_Z << ": " << false << endl;
  	}
    return 0;
} // main()

