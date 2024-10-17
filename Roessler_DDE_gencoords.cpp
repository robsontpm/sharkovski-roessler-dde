/* ========================================================================
 * Roessler DDE system
 * Anna Gierzkiewicz and Robert Szczelina
 * ========================================================================
 * Program: Roessler_DDE_gencoords.cpp
 * ========================================================================
 * The program contains the method to generate structure of the initial
 * segments of the DDE problem to be feed into machinery later.
 * As the structure of initial set is basically infinite dimensional
 * in a theoretical sense and of very large dimension in rigorous numerics
 * it is not feasible to provide initial data by hand. This is why we need
 * an automatic procedure.
 *
 * The default values of parameters to this program, ale set to match
 * the method used in the accompanying paper. The reader can experiment with
 * other settings, if they wishes.
 * ========================================================================*/ 

#include "setup.h"

/**
 * the function to generate good enough coordinates empirically for a given initial data
 * (point on ection in the ODE. The coordiates are for the segment in the DDE system.
 */
void computeCoordsForward(
		bool verbose,
		system3d& system, DVector const& ivp,
		HSet2D& hset, int CUTS_Y, int CUTS_Z,
		std::vector<DVector>& attractor_out, DSolution& solution_out,
		IVector& x0_out, IMatrix& C_out, IVector& r0_out);

/**
 * This saves a representation of a box in good coordinates to a dat file
 * (nonrigorous). This is just to use them in ploting pretty figures
 * later.
 */
void save_box(std::string const& filepath, IVector const& ibox);

/**
 * This doeas kind of nonrigorous visualization
 * of the set given in the good coords.
 */
void nonrig_box_plot(
		std::string plotdir,
		std::vector<DVector> const& attractor,
		system3d& system, HSet2D& hset, int CUTS_Y, int CUTS_Z,
		DVector const& ref, DVector const& ydir, DVector const& zdir);

// Main function
int main(int argc, char* argv[])
{
	// this is a helper class to parse arguments from command line
	// we will prepare a lot of commands and run them in parallel to get the
	// speed boost. Each computation will check the conditions of the small
	// chunk of the original set.
	capd::ddeshelper::ArgumentParser parser(argc, argv);
	int CUTS_Y = 100, CUTS_Z = 3;
	string WD = "./tmp/";
	bool verbose = false;
	// we use setup with epsi=0 in this step.
	// the generated coordinates should be good enough for small epsilons.
	double EPSI = 0.;
	// initial value for the nonrigorous computations.
	// It should be such that the trajectory is quite dense in the attractor.
	// Might be tricky for periodic-windows parameter values of the system.
	DVector IVP {0.,-5,0.03};
	parser.parse("ivp=", IVP, "initial value for the simulation. ");
	parser.parse("epsi=", EPSI, "epsi value for the computation.");
	parser.parse("cuty=", CUTS_Y, "how many pieces in y");
	parser.parse("cutz=", CUTS_Z, "how many pieces in z");
	parser.parse("wd=", WD, "the working directory of the program, i.e. where to look for data and to put results. Should end with a '/'.");
	verbose = parser.parse("verbose", std::string("print more output")) || parser.parse("debug", std::string("same as verbose"));
	if (parser.isHelpRequested()){
		std::cout << parser.getHelp() << endl;
		return 0;
	}

	// sanitize input
	if (WD == "") WD = "./";
	if (WD.back() != '/') WD += "/";
	// here we will save plots / gnuplot files for figures
	std::string plotdir = WD + "plots/";
	capd::ddeshelper::mkdir_p(plotdir);

  	cout.precision(16);	// precision of nonrigorous output
	cout << boolalpha;  // outputting bools as true/false
	SHA_DEBUG("# " << parser.getCommandLine()) << endl;
	SHA_DEBUG("# The working Directory is '" << WD << "'") << endl;

	try
	{
		// get the configuration of the Roessler system with the selected parameters.
		Config local_config(EPSI);
		// the names are inherited from paper [GZ2022]
		auto& roessler525 = local_config.roessler525;
		auto& c3 = local_config.c3;
		auto& grid3 = local_config.grid3;
		auto& grid = roessler525.grid;
		auto& order = roessler525.order;

		// here we will store the resulting coordinate system
		IVector I_x0(roessler525.M());
		IMatrix I_C(roessler525.M(), roessler525.M());
		IVector I_r0(roessler525.M());
		// here we will hold the attractor segments on the section, for visualization
		std::vector<DVector> attractor;
		// this will hold the whole solution for generating the attractor
		// i.e. it will look like the Roesler attractor in 3D.
		DSolution solution = roessler525.makeDSegment(DVector(roessler525.M()));

		// this is a point that produces reasonably dense trajectory in the ODE for selected parameters
		// and based on that it computes a reasonably good coordinates for the sets.
		// it will save the coordinates in various files.
		cout << "GENERATING coords... " << std::endl;
		computeCoordsForward(
			verbose,
			roessler525,
			IVP, grid3, CUTS_Y, CUTS_Z,
			attractor, solution,
			I_x0, I_C, I_r0
		);
		// do some fancy plots
		splotMany(roessler525, attractor, "forward-history", plotdir);
		std::string filename = "trajectory";
		capd::ddeshelper::plot_value(plotdir + filename, solution.pastTime(), solution.currentTime(), roessler525.h.leftBound() / 10, solution, false);
		std::ostringstream plot;
		plot << "'" << filename << "ddes-plot.dat' u 3:5:7 with lines notitle";
		capd::ddeshelper::splot_many(plotdir, { plot.str() }, false, "trajectory");
		std::ostringstream cmd; cmd << "cd '" << plotdir << "' && gnuplot 'trajectory.gp'";
		capd::ddeshelper::runSystemCommand(cmd.str());
		// end of fancy plots
		cout << "DONE GENERATING coords" << std::endl;

		// Now we compute translation of the mid points of the sets c[i].
		auto ydir = I_C.column(1);
		auto zdir = I_C.column(2);
		// make a nonrig plot of the set
		nonrig_box_plot(
				plotdir, attractor, roessler525,
				grid3, CUTS_Y, CUTS_Z,
				capd::vectalg::midObject<DVector>(I_x0),
				capd::vectalg::midObject<DVector>(ydir),
				capd::vectalg::midObject<DVector>(zdir));

		cout << "GENERATING C_i sets... " << std::endl;
		// TODO: rethink!
		// TODO: now that I have changed the head to match the old hset grid3,
		// TODO: we should just now hmmm... what to do?
		// Attractor's container (TODO: rethink!) This should be just grid3 ??
		//HSet2D gridAI (IVector ({I_x0[1], I_x0[2]}), IMatrix ({{I_C[1][1], I_C[1][2]}, {I_C[2][1], I_C[2][2]}}), DVector({3.63687,0.0004}));
		//auto& refgrid = gridAI;
		//auto& refgrid = grid3;
		HSet2D refgrid(grid3.get_I_x(), grid3.get_I_B(), grid3.get_r());

		// transform the reference points of other sets (C_i's) to the good coordinate frame
		SHA_DEBUG("refgrid.B = " << refgrid.get_I_B() << endl << "B^{-1} = " << refgrid.get_I_invB() << endl);
		auto ref = refgrid.get_I_x();
		auto invC = refgrid.get_I_invB();
		std::vector<IVector> I_c3;
		for (int i = 0; i < c3.size(); i++){
			auto c3_mid = c3[i].get_I_x();
			auto diff = invC * (c3_mid - ref);
			SHA_DEBUG("C[" << i << "].mid in coordinates: " << diff << endl);
			auto new_c3 = capd::vectalg::midObject<DVector>(I_x0 + ydir * diff[0] + zdir * diff[1]);
			SHA_DEBUG("head: ");
			new_c3[0] = 0.; // x = 0 on section
			IVector v2d {new_c3[1], new_c3[2]};
			SHA_DEBUG(v2d << endl << "diff:" << v2d - c3_mid << endl);

			std::ostringstream fileprefix;
			fileprefix << "C_" << i;
			capd::ddeshelper::saveBinary(WD + fileprefix.str() + "_x0.ivector.bin", IVector(new_c3));
			I_c3.push_back(IVector(new_c3));

			save_box(plotdir + "BOX_" + fileprefix.str() + ".dat", diff + c3[i].box());
			cout << "C_" << i << " transformation done" << endl;
		}
		cout << "DONE C_i sets... " << std::endl;

		cout << "SAVING coords... " << std::flush;
		// save data to do pretty pictures
		save_box(plotdir + "BOX_G.dat", grid3.box());
		// save machine representation of numbers
		capd::ddeshelper::saveBinary(WD + "G_x0.ivector.bin", I_x0);
		capd::ddeshelper::saveBinary(WD + "G_M.imatrix.bin", I_C);
		capd::ddeshelper::saveBinary(WD + "G_r0.ivector.bin", I_r0);
		// save human-readable representations (it might get machine precision loss)
		std::ofstream human_readable_coords("human_readable_coords.txt");
		human_readable_coords.precision(16);
		human_readable_coords << "G_x0: " << endl;
		human_readable_coords << capd::vectalg::midObject<DVector>(I_x0) << endl;
		for (int i = 0; i < I_c3.size(); ++i){
			human_readable_coords << "C_" << i << "_x0: " << endl;
			human_readable_coords << capd::vectalg::midObject<DVector>(I_c3[i]) << endl;
		}
		human_readable_coords << "M: " << endl;
		human_readable_coords << capd::vectalg::midObject<DMatrix>(I_C) << endl;
		human_readable_coords << "r0: " << endl;
		human_readable_coords << capd::vectalg::midObject<DVector>(I_r0) << endl;
		std::cout << "DONE. Data saved in '" << WD << "'" << endl;
		std::cout << "Please run from '" << WD << "' folder the following: " << endl;
		std::cout << endl;
		std::cout << "    {capdDDEsDIR}/bin/convmatrix ";
		std::cout << roessler525.M() << " bin G_M.imatrix.bin ";
		std::cout << "inv bin " << "G_invM.imatrix.bin" << endl;
		std::cout << endl;
		std::cout << "NOTE: you need to determine the relative/absolute path of convmatrix script w.r.t. '" << WD << "' path!" << endl;
	}
	catch(exception& e)
  	{
    	cout << "\n\nException caught: "<< e.what() << endl;
  	}
  return 0;
} // main()


void computeCoordsForward(
		bool verbose,
		system3d& system, DVector const& ivp,
		HSet2D& hset, int CUTS_Y, int CUTS_Z,
		std::vector<DVector>& attractor_out, DSolution& solution_out,
		IVector& x0_out, IMatrix& C_out, IVector& r0_out){

	// renaming for convenience
	int& p = system.p;
	int& order = system.order;
	int d = system.d;
	auto& attractor = attractor_out;
	DGrid &grid = system.dgrid;

	// make the non-rigorous integrator for DDE
	DEq rhs(system.a.leftBound(), 0.2, system.eps.mid().leftBound());
	DDDEq vf(rhs, grid(p));
	DSection section(3, 0, 0.);
	DSolver solver(vf, order * 3);
	DPoincare P(solver, section, poincare::MinusPlus);

	// make data needed to propagate the initial fanction, which is constant.
	DSolution constIVP(grid, -grid(p), grid(0), order, ivp);
	DSolution Px(grid, -grid(p), grid(0), order, ivp);
	solution_out = constIVP;

	// generate a lot of points (segments) on the section
	double tp;
	int NUM_ITERS = 100; // TODO: this should be a parameter

	// please note that, we do not need to exclude few initial
	// iterates, as the selected point should be on the attractor (or close)
	// since we are (by default) doing iterates for epsi=0
	// then the tail will be automatically right, after first iterate
	// as it is just a solution to ODE (back in time from the point on the section)
	for (int i = 0; i < NUM_ITERS; ++i){
		P(solution_out, Px, tp);
		DVector pv = Px;
		pv[0] = 0.0;
		attractor.push_back(pv);
	}
	// the solution_out will contain the whole trajectory of the point after 100 iterates of the poincare map!

	// now we start to generate good coordinate frame
	DVector zerov(system.M()); // make a zero vector, as I use that frequently to initialize data

	// this is a simple function to compute average
	auto mean = [zerov](std::vector<DVector> items){
		return items.size() ?
			std::accumulate(items.begin(), items.end(), zerov) / double(items.size())
			: zerov;
	};

	// get the representation of the head set from the finite
	// dimensional set on the section. NoteL it it 2D vector!
	DMatrix fin_C = hset.get_B();
	DVector fin_mid_main = expand(hset.get_x());

	// expand them with the x direction.
	DVector fin_ydir = expand((DVector)fin_C.column(0));
	DVector fin_zdir = expand((DVector)fin_C.column(1));
	SHA_DEBUG(fin_ydir << endl << fin_zdir << endl);

	// sort data with respect to y-coordinate
	std::sort(
		attractor.begin(), attractor.end(),
		[](DVector const& a, DVector const& b){
			return a[1] < b[1];
		}
	);

	// this seems like bad idea to compute average direction
	// instead, i just connect two farthest points and this seems ok-ish
	SHA_DEBUG("Computing ydir" << endl);
	DVector ydir = attractor.front() - attractor.back();
	double norm = abs(ydir[1]); ydir /= norm;

	// we compute mid points among all the pairs
	// from the beginning and the end of the sorted attractor
	// that way, we should have a point that is close to the
	// original hset reference point.
	SHA_DEBUG("Computing mid point" << endl);
	std::vector<DVector> mid_candidates;
	for (int i = 0; i < attractor.size() / 3; ++i)
		for (int j = attractor.size() -1; j >= 2 * attractor.size() / 3; --j)
			mid_candidates.push_back(mean({attractor[i], attractor[j]}));

	// now, we try to find the point, among selected
	// that is as close as possible to the original reference
	// point of the set
	capd::vectalg::SumNorm<DVector, DMatrix> selnorm;
	std::sort(mid_candidates.begin(), mid_candidates.end(), [&fin_mid_main, &selnorm](DVector const& a, DVector const& b){
		DVector aa{a[1], a[2]};
		DVector bb{b[1], b[2]};
		DVector cc{fin_mid_main[1], fin_mid_main[2]};

		return selnorm(aa-cc) < selnorm(bb-cc);
	});
	DVector mid_point = mid_candidates[0];

	// output some debug information, when needed
	auto compheads = [&](DVector const& fin_v, DVector const& full_v){
		DVector proj_v {full_v[0], full_v[1], full_v[2]};
		SHA_DEBUG("old: " << fin_v << endl);
		SHA_DEBUG("new: " << proj_v << endl);
		SHA_DEBUG("dif: " << proj_v - fin_v);
	};
	SHA_DEBUG("midpoint:" << endl); compheads(fin_mid_main, mid_point); SHA_DEBUG(endl);
	SHA_DEBUG("ydir:    " << endl); compheads(fin_ydir, ydir); SHA_DEBUG(endl);

	// we select simple z direction, to be the same as in the original, 3-dimensional set
	DVector zdir(system.M());
	zdir[1] = fin_zdir[1];
	zdir[2] = fin_zdir[2];

	// change the values on head, so they match the original hset
	ydir[1] = fin_ydir[1];
	ydir[2] = fin_ydir[2];
	mid_point[1] = fin_mid_main[1];
	mid_point[2] = fin_mid_main[2];
	// TODO: test those changes!

	// save coordinates to be used later
	x0_out = IVector(mid_point);
	C_out = IMatrix(system.M(), system.M());
	C_out.setToIdentity();
	C_out.column(1) = IVector(ydir);
	C_out.column(2) = IVector(zdir);
	r0_out = IVector(system.M()); r0_out *= 0.;
	// TODO: estimate r) based on the attractor?
}

void save_box(std::string const& filepath, IVector const& ibox){
	ostringstream prefix;
	ofstream dat(filepath);
	dat.precision(16);
	dat << capd::ddeshelper::to_dat(ibox[0]) << " ";
	dat << capd::ddeshelper::to_dat(ibox[1]) << " ";
	dat << capd::ddeshelper::to_dat(interval(-1, 1)) << " ";
	dat << capd::ddeshelper::to_dat(interval(-1, 1)) << " ";
	dat.close();
}

void nonrig_box_plot(
		std::string plotdir,
		std::vector<DVector> const& attractor,
		system3d& system, HSet2D& hset, int CUTS_Y, int CUTS_Z,
		DVector const& ref, DVector const& ydir, DVector const& zdir){

	auto& order = system.order;
	auto& d = system.d;
	auto& p = system.p;
	GridSet gridSet(2);
	hset.gridSet(gridSet, CUTS_Y, CUTS_Z);
	IVector box = expand(gridSet.box());
	std::vector<DVector> the_set_chunks;
	DMatrix gridInvCoords = capd::vectalg::midObject<DMatrix>(hset.get_I_invB());
	for(auto chunk = gridSet.begin(); chunk != gridSet.end(); ++chunk)
	{
		DVector mid_chunk = capd::vectalg::midObject<DVector>(*chunk);
		DVector diff = gridInvCoords*(mid_chunk - hset.get_x());
		DVector v = ref + diff[0] * ydir + diff[1] * zdir;
		the_set_chunks.push_back(v);
	}
	splotMany(system, the_set_chunks, "the_set_chunks_forward", plotdir);

	// for each time point draw (a 3d?) picture of all the points on the
	// generated set and all the points of the attractor
	// we expect that points of the attractor should more or less lie inside the
	// set points
	{
		auto point3d = [&](DVector const& item, int i, int k){
			int bi = (i == 0 ? 0 : d * (1 + (i-1) * (order+1) + k)); // baseindex
			DVector v { item[bi], item[bi+1], item[bi+2] };
			return v;
		};
		auto savedat = [&](std::string const& path, std::vector<DVector> const& items, int i, int k){
			ofstream datpoints(path);
			for (auto const& item: items){
				auto p3d = point3d(item, i, k);
				for (int j = 0; j < d; ++j)
					datpoints << p3d[j] << " ";
				datpoints << "\n";
			}
			datpoints.close();
		};

		std::string plotpath = plotdir + "scatter/";
		capd::ddeshelper::mkdir_p(plotpath);
		for (int k = 0; k <= order; ++k){
			for (int i = 0; i < p; ++i){
				if (i == 0 && k > 0) continue;

				std::ostringstream idprefix;
				idprefix << "-k" << k << "-t" << i;

				savedat(plotpath + "attractor" + idprefix.str() + ".dat", attractor, i, k);
				savedat(plotpath + "gridset" + idprefix.str() + ".dat", the_set_chunks, i, k);

				ostringstream splot_attractor;
				splot_attractor << "'attractor" << idprefix.str() << ".dat" << "' with points pt 7 ps 2 lc rgb 'red'";
				ostringstream splot_gridset;
				splot_gridset << "'gridset" << idprefix.str() << ".dat" << "' with points pt 7 ps 1 lc rgb 'black'";


				capd::ddeshelper::splot_many(
					plotpath + "live" + idprefix.str(),
					{
						splot_attractor.str(),
						splot_gridset.str(),
					},
					true, ""
				);
				capd::ddeshelper::splot_many(
					plotpath + "all" + idprefix.str(),
					{
						splot_attractor.str(),
						splot_gridset.str(),
					},
					false, ""
				);
				// force plot, as xplot is precompiled without the flag...
				std::ostringstream cmd; cmd << "cd '" << plotpath << "' && gnuplot 'all" << idprefix.str() << ".gp" << "'";
				capd::ddeshelper::runSystemCommand(cmd.str());
			}
		}

	}
}
