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

/** this is just for debug purposes and should be deprecated later */
void printHistory(typename system3d::HistoryType const& history);
/** this is just for debug purposes and should be deprecated later */
void printBase(DVector const& x0, std::vector<DVector> const& coords, int i);
/** TODO: docs */
void splotMany(system3d& roessler525, std::vector<DVector>& x0_all, std::string const& name);
/** TODO: docs */
void computeCoordsByC1(system3d&, HSet2D&, int, int);
/** TODO: docs */
void computeCoordsSimple(system3d&, HSet2D&, int, int);
///** TODO: docs */
//void computeCoordsForward(system3d&, DVector const& ivp);
void computeCoordsForward(
		system3d& system, DVector const& ivp,
		HSet2D& hset, int CUTS_Y, int CUTS_Z,
		std::vector<DVector> fin_mid_others = {});

// Main function

int main()
{
  	cout.precision(16);
	cout << boolalpha;  
	try
	{			
		interval a = 5.25;
		interval epsi = interval(0.0);
		//interval epsi = interval(0.0001);
		interval tau = 0.5;
		system3d roessler525(tau, epsi, a);
		///===================== variables used in Procedure 1:  =====================
		HSet2D grid3 (IVector ({-6.38401, 0.0327544}) , IMatrix ({{-1., 0.000656767}, {-0.000656767, -1.}}) , DVector({3.63687,0.0004}));			// Attractor's container
		vector<HSet2D> c3(3);
		c3[0] = HSet2D(IVector ({-3.46642, 0.0346316}) , IMatrix ({{-1., 0.000656767}, {-0.000656767, -1.}}) , DVector({0.072,0.00048}));			// cube 1
		c3[1] = HSet2D(IVector ({-6.26401, 0.0326544}) , IMatrix ({{-1., 0.000656767}, {-0.000656767, -1.}}) , DVector({0.162,0.00066}));			// cube 2
		c3[2] = HSet2D(IVector ({-9.74889, 0.0307529}) , IMatrix ({{-1., 0.000656767}, {-0.000656767, -1.}}) , DVector({0.036,0.00072}));			// cube 3
		

		const int CUTS_Y = 100;
		const int CUTS_Z = 3;

		DVector ivp {0.,-5,0.03};
		// TODO: UNCOMMENT FOR LATER
		computeCoordsForward(roessler525, ivp, grid3, CUTS_Y, CUTS_Z);

		// compute translation of the mid points of the sets c[i]:
		IVector I_x0(roessler525.M());
		IMatrix I_C(roessler525.M(), roessler525.M());
		capd::ddeshelper::readBinary("grid3_x0.ivector.bin", I_x0);
		capd::ddeshelper::readBinary("grid3_C.imatrix.bin", I_C);

		auto ydir = I_C.column(1);
		auto zdir = I_C.column(2);


		HSet2D gridAI (IVector ({I_x0[1], I_x0[2]}) , IMatrix ({{I_C[1][1], I_C[1][2]}, {I_C[2][1], I_C[2][2]}}) , DVector({3.63687,0.0004}));			// Attractor's container

		// save data to do pretty pictures
		IVector box = grid3.box();
		ostringstream prefix;
		prefix << "BOX_grid3";
		string dirpath = "plots/";
		capd::ddeshelper::mkdir_p(dirpath);
		ofstream dat(dirpath + prefix.str() + ".dat");
		dat.precision(16);
		dat << capd::ddeshelper::to_dat(box[0]) << " ";
		dat << capd::ddeshelper::to_dat(box[1]) << " ";
		dat << capd::ddeshelper::to_dat(interval(-1, 1)) << " ";
		dat << capd::ddeshelper::to_dat(interval(-1, 1)) << " ";
		dat.close();

		auto& refgrid = gridAI;
		//auto& refgrid = grid3;

		cout << refgrid.get_I_B() << endl << refgrid.get_I_invB() << endl;
		auto ref = refgrid.get_I_x();
		auto invC = refgrid.get_I_invB();
		for (int i = 0; i < 3; i++){
			auto c3_mid = c3[i].get_I_x();
			auto diff = invC * (c3_mid - ref);
			cout << "c3[" << i << "].mid in coordinates: " << diff << endl;
			auto new_c3 = capd::vectalg::midObject<DVector>(I_x0 + ydir * diff[0] + zdir * diff[1]);
			cout << "head: " << new_c3[0] << " " ;
			new_c3[0] = 0.; // x = 0 on section
			IVector v2d {new_c3[1], new_c3[2]};
			cout << v2d << endl;
			cout << "diff:" << v2d - c3_mid << endl;

			std::ostringstream filename;
			filename << "c3_" << i << "_x0.ivector.bin";
			capd::ddeshelper::saveBinary(filename.str(), IVector(new_c3));


			// save data to do pretty pictures
			IVector box = diff + c3[i].box();
			ostringstream prefix;
			prefix << "BOX_c3_" << i;
			ofstream dat(dirpath + prefix.str() + ".dat");
			dat.precision(16);
			dat << capd::ddeshelper::to_dat(box[0]) << " ";
			dat << capd::ddeshelper::to_dat(box[1]) << " ";
			dat << capd::ddeshelper::to_dat(interval(-1, 1)) << " ";
			dat << capd::ddeshelper::to_dat(interval(-1, 1)) << " ";
			dat.close();

			cout << endl;
		}



		return 0;
	}
	catch(exception& e)
  	{
    	cout << "\n\nException caught: "<< e.what() << endl;
  	}
  return 0;
} // main()


// Auxiliary Functions

void printHistory(typename system3d::HistoryType const& history){
	for (auto jet: history){
		for (auto coeff: jet)
			cout << coeff.first << "\n" << coeff.second << "\n\n";
		cout << "NEXT" << endl;
	}
}

void printBase(DVector const& x0, std::vector<DVector> const& coords, int i){
	std::cout << x0 << endl << endl;
	std::cout << "x0[" << i << "]: " << x0 << endl << endl;
	for (int j = 0; j < 3; ++j)
		std::cout << "C[" << i << "]_*" << j << ": " << coords[j] << endl << endl;
}

void splotMany(system3d& roessler525, std::vector<DVector>& x0_all, std::string const& name){
	// plot all the solutions in the Ambient space
	const std::string plotpath = "./plots/" + name + "/";
	capd::ddeshelper::mkdir_p(plotpath);
	std::vector<std::string> gplots;
	std::cerr << "Plotting " << x0_all.size() << " segments... " << std::flush;
	for (int i = 0; i < x0_all.size(); ++i){
		DVector& v = x0_all[i];
		auto seg = roessler525.makeDSegment(v);
		std::ostringstream filename; filename << "seg-" << i;
		std::string fileprefix = plotpath + filename.str();
		// cout << v << endl;
		// cout << "plotting " << i << " " << fileprefix << endl;
		// capd::ddeshelper::plot_value(fileprefix.str(), seg, false);
		capd::ddeshelper::plot_value(fileprefix, seg.pastTime(), seg.currentTime(), roessler525.h.leftBound() / 10, seg, false);
		std::ostringstream plot;
		plot << "'" << filename.str() << "ddes-plot.dat' u 3:5:7 with lines notitle";
		gplots.push_back(plot.str());
	}
	std::cerr << "DONE" << std::endl;

	capd::ddeshelper::splot_many(plotpath, gplots, false, "all");
	// force plot, as xplot is precompiled without the flag...
	std::ostringstream cmd; cmd << "cd '" << plotpath << "' && gnuplot 'all.gp'";
	capd::ddeshelper::runSystemCommand(cmd.str());
}

/** TODO: Docs */
void computeCoordsByC1(system3d& roessler525, HSet2D& grid3, int CUTS_Y, int CUTS_Z){

	// TODO: coś mi tu nie poszły te obliczenia... przedyskutowac i przemyśleć w przyszłości...

	// roessler525.makeHistory(c3[0]);

	GridSet gridSet(2);
	grid3.gridSet(gridSet,CUTS_Y,CUTS_Z);
	auto box = gridSet.box();
	auto coords = gridSet.coordinateSystem();

	DVector zerov(roessler525.M());
	std::vector<std::vector<DVector>> finite_coords_all(gridSet.size(), {zerov, zerov, zerov});
	std::vector<DVector> x0_all(gridSet.size(), zerov);

	auto chunk = gridSet.begin();
	for(int i=0; chunk != gridSet.end(); ++chunk, ++i)
	{
		DMatrix C = capd::vectalg::midObject<DMatrix>(expand(coords));
		DVector v0 = capd::vectalg::midObject<DVector>(expand(*chunk));
		auto history = roessler525.makeHistoryD(v0, C, x0_all[i], finite_coords_all[i]);
//			printHistory(history);
//			printBase(x0_all[i], finite_coords_all[i], i);
	}
	// plot all the solutions in the Ambient space
	splotMany(roessler525, x0_all, "history");

	// Compute average coordinates and average middle
	// TODO: rethink summation and dividing, this might not be the most numerically accurate method...
	auto mean = [zerov](std::vector<DVector> items){
		auto vv = std::accumulate(items.begin(), items.end(), zerov);
		vv /= double(items.size());
		return vv;
	};
	DVector x0 = mean(x0_all);
	std::vector<DVector> coords_finite;
	for (int j = 0; j < 3; j++){
		std::vector<DVector> jcoords;
		for (auto const& coordv: finite_coords_all)
			jcoords.push_back(coordv[j]);
		coords_finite.push_back(mean(jcoords));
	}

	// Compute differences between reps and all vectors, plot scatterplot
	// First, i need 3 collections of many vectors, not many collections of 3 vectors.. (transpose the result)
	// I miss python :(
	std::vector<std::vector<DVector>> jcoords;
	for (int j = 0; j < 3; j++){
		std::vector<DVector> tmp;
		for (auto const& item: finite_coords_all)
			tmp.push_back(item[j]);
		jcoords.push_back(tmp);
	}
	// then prepare computations
	capd::vectalg::MaxLNorm<DVector, DMatrix> maxnorm;
	capd::vectalg::SumLNorm<DVector, DMatrix> sumnorm;
	capd::vectalg::EuclLNorm<DVector, DMatrix> euclnorm;
	std::vector<capd::vectalg::Norm<DVector, DMatrix>* > norms({&maxnorm, &sumnorm, &euclnorm});
	auto scatterPlot = [zerov, norms](std::ostream &out, DVector const& ref, std::vector<DVector> items){
		for (auto const& item: items){
			auto diff = item - ref;
			for (auto norm: norms)
				out << (*norm)(diff) << " ";
			out << diff << endl;
		}
	};
	// now compute the statistics
	cout << "x0: " << endl;
	scatterPlot(cout, x0, x0_all);
	cout << endl;
	for (int j =0; j < 3; ++j){
		cout << "C_{*" << j << "}: " << endl;
		scatterPlot(cout, coords_finite[j], jcoords[j]);
		cout << endl;
	}

	cout << "DIM: " << roessler525.M() << endl;

	// TODO: save coordinates for future use
	// x0 is the middle and coords_finite[j] is the j-th coordinate
	// we use only j=1 and 2, not 0 (section propagated)

	// Draw the grid set but not using history but coords:
	std::vector<DVector> the_set_chunks;
	the_set_chunks.push_back(x0);

	DVector ref = grid3.get_x();
//		auto selected_coords = coords_finite;
//		auto selected_x0 = x0;

	int sel = 0;
	auto selected_coords = finite_coords_all[0];
	auto selected_x0 = x0_all[0];
	int cnt = 0;
	for(auto chunk = gridSet.begin(); chunk != gridSet.end(); ++chunk)
		if (cnt == sel){
			ref = capd::vectalg::midObject<DVector>(*chunk);
			break;
		}
	DMatrix gridInvCoords = capd::vectalg::midObject<DMatrix>(grid3.get_I_invB());
	for(auto chunk = gridSet.begin(); chunk != gridSet.end(); ++chunk)
	{
		DVector mid_chunk = capd::vectalg::midObject<DVector>(*chunk);
		DVector diff = gridInvCoords*(mid_chunk - ref);
		cout << diff << endl;
		DVector v = x0 + diff[0] * coords_finite[1] + diff[1] * coords_finite[2];
		the_set_chunks.push_back(v);
	}
	splotMany(roessler525, the_set_chunks, "the_set_chunks");



	DMatrix infinite_coords;
}

void computeCoordsSimple(system3d& roessler525, HSet2D& grid3, int CUTS_Y, int CUTS_Z){

//	GridSet gridSet(2);
//	DMatrix grid_coords = capd::vectalg::midObject<DMatrix>(expand(gridSet.coordinateSystem()));
//	DVector grid_v0 = capd::vectalg::midObject<DVector>(expand(grid3.get_x()));
//	IVector grid_r0 = expand(grid3.box());
//	auto l_dy = grid_r0[1].leftBound();
//	auto r_dy = grid_r0[1].rightBound();
//	auto l_dz = grid_r0[2].leftBound();
//	auto r_dz = grid_r0[2].rightBound();

	GridSet gridSet(2);
	grid3.gridSet(gridSet, CUTS_Y, CUTS_Z);
	IVector box = expand(gridSet.box());
	DMatrix coords = capd::vectalg::midObject<DMatrix>(expand(gridSet.coordinateSystem()));
	auto l_dy = box[1].leftBound();
	auto r_dy = box[1].rightBound();
	auto l_dz = box[2].leftBound();
	auto r_dz = box[2].rightBound();
	//DVector set_mid = capd::vectalg::midObject<DVector>(expand(grid3.get_x()));
	DVector set_mid = expand(grid3.get_x());

	DVector zerov(roessler525.M());
	std::vector<DVector> x0_all(gridSet.size(), zerov);
	std::vector<DVector> zdirs, ydirs;
	auto chunk = gridSet.begin();
	for(int i=0; chunk != gridSet.end(); ++chunk, ++i)
	{
		DMatrix C = coords;
		DVector v0 = capd::vectalg::midObject<DVector>(expand(*chunk));
		roessler525.makeHistoryC0(v0, x0_all[i]);
		auto approx_direction = [&zerov, &v0, &roessler525](double l_d, double r_d, DVector coord){
			DVector l_coord(zerov), r_coord(zerov);
			roessler525.makeHistoryC0(v0 + l_d * coord, l_coord);
			roessler525.makeHistoryC0(v0 + r_d * coord, r_coord);
			return (r_coord - l_coord) / (r_d - l_d);
		};
		DVector ydir = approx_direction(l_dy, r_dy, coords.column(1));
		DVector zdir = approx_direction(l_dz, r_dz, coords.column(2));
		zdirs.push_back(zdir);
		ydirs.push_back(ydir);
	}
	// plot all the solutions in the Ambient space
	splotMany(roessler525, x0_all, "history-C0");

	// Compute average coordinates
	// The middle point is just just the image of the mid
	DVector x0(zerov);
	roessler525.makeHistoryC0(set_mid, x0);

	// TODO: rethink more stable numerical method to compute mean vector
	auto mean = [zerov](std::vector<DVector> items){
		auto vv = std::accumulate(items.begin(), items.end(), zerov);
		vv /= double(items.size());
		return vv;
	};
	// now we compute mean directions!
	DVector ydir = mean(ydirs);
	DVector zdir = mean(zdirs);
	// compute diffs from the mean (they should be small, relatively?
	capd::vectalg::MaxLNorm<DVector, DMatrix> maxnorm;
	capd::vectalg::SumLNorm<DVector, DMatrix> sumnorm;
	capd::vectalg::EuclLNorm<DVector, DMatrix> euclnorm;
	std::vector<capd::vectalg::Norm<DVector, DMatrix>* > norms({&maxnorm, &sumnorm, &euclnorm});
	auto scatterPlot = [zerov, norms](std::ostream &out, DVector const& ref, std::vector<DVector> items){
		for (auto const& item: items){
			auto diff = item - ref;
			for (auto norm: norms)
				out << (*norm)(diff) << " ";
			out << diff << endl;
		}
	};
	cout << "ydirs: " << endl; scatterPlot(cout, ydir, ydirs);
	cout << "zdirs: " << endl; scatterPlot(cout, zdir, zdirs);

	std::vector<DVector> the_set_chunks;
	DMatrix gridInvCoords = capd::vectalg::midObject<DMatrix>(grid3.get_I_invB());
	for(auto chunk = gridSet.begin(); chunk != gridSet.end(); ++chunk)
	{
		DVector mid_chunk = capd::vectalg::midObject<DVector>(*chunk);
		DVector diff = gridInvCoords*(mid_chunk - grid3.get_x());
		cout << diff << endl;
		DVector v = x0 + diff[0] * ydir + diff[1] * zdir;
		the_set_chunks.push_back(v);
	}
	splotMany(roessler525, the_set_chunks, "the_set_chunks_C0");

	// save coordinates to be used later
	IVector I_x0(x0);
	IMatrix I_C(roessler525.M(), roessler525.M());
	I_C.setToIdentity();
	I_C.column(1) = IVector(ydir);
	I_C.column(2) = IVector(zdir);
	IVector I_r0(roessler525.M()); I_r0 *= 0.;

	std::ofstream human_readable_coords("human_readable_coords.txt");
	human_readable_coords << "x0: " << endl;
	human_readable_coords << x0 << endl;
	human_readable_coords << "C: " << endl;
	human_readable_coords << capd::vectalg::midObject<DMatrix>(I_C) << endl;
	human_readable_coords << "r0: " << endl;
	human_readable_coords << I_r0 << endl;

	capd::ddeshelper::saveBinary("grid3_x0.ivector.bin", I_x0);
	capd::ddeshelper::saveBinary("grid3_C.imatrix.bin", I_C);
	capd::ddeshelper::saveBinary("grid3_r0.ivector.bin", I_r0);
	std::cout << "Please run from 'bin' folder: " << endl;
	std::cout << "../../../../bin/convmatrix ";
	std::cout << roessler525.M() << " bin grid3_C.imatrix.bin ";
	std::cout << "inv bin " << "grid3_invC.imatrix.bin" << endl;


	//DVector x0 = mean(x0_all);
//	std::vector<DVector> coords_finite;
//	for (int j = 0; j < 3; j++){
//		std::vector<DVector> jcoords;
//		for (auto const& coordv: finite_coords_all)
//			jcoords.push_back(coordv[j]);
//		coords_finite.push_back(mean(jcoords));
//	}
//
//	// Compute differences between reps and all vectors, plot scatterplot
//	// First, i need 3 collections of many vectors, not many collections of 3 vectors.. (transpose the result)
//	// I miss python :(
//	std::vector<std::vector<DVector>> jcoords;
//	for (int j = 0; j < 3; j++){
//		std::vector<DVector> tmp;
//		for (auto const& item: finite_coords_all)
//			tmp.push_back(item[j]);
//		jcoords.push_back(tmp);
//	}
//	// then prepare computations
//	capd::vectalg::MaxLNorm<DVector, DMatrix> maxnorm;
//	capd::vectalg::SumLNorm<DVector, DMatrix> sumnorm;
//	capd::vectalg::EuclLNorm<DVector, DMatrix> euclnorm;
//	std::vector<capd::vectalg::Norm<DVector, DMatrix>* > norms({&maxnorm, &sumnorm, &euclnorm});
//	auto scatterPlot = [zerov, norms](std::ostream &out, DVector const& ref, std::vector<DVector> items){
//		for (auto const& item: items){
//			auto diff = item - ref;
//			for (auto norm: norms)
//				out << (*norm)(diff) << " ";
//			out << diff << endl;
//		}
//	};
//	// now compute the statistics
//	cout << "x0: " << endl;
//	scatterPlot(cout, x0, x0_all);
//	cout << endl;
//	for (int j =0; j < 3; ++j){
//		cout << "C_{*" << j << "}: " << endl;
//		scatterPlot(cout, coords_finite[j], jcoords[j]);
//		cout << endl;
//	}
//
//	cout << "DIM: " << roessler525.M() << endl;
//
//	// TODO: save coordinates for future use
//	// x0 is the middle and coords_finite[j] is the j-th coordinate
//	// we use only j=1 and 2, not 0 (section propagated)
//
//	// Draw the grid set but not using history but coords:
//	std::vector<DVector> the_set_chunks;
//	the_set_chunks.push_back(x0);
//
//	DVector ref = grid3.get_x();
////		auto selected_coords = coords_finite;
////		auto selected_x0 = x0;
//
//	int sel = 0;
//	auto selected_coords = finite_coords_all[0];
//	auto selected_x0 = x0_all[0];
//	int cnt = 0;
//	for(auto chunk = gridSet.begin(); chunk != gridSet.end(); ++chunk)
//		if (cnt == sel){
//			ref = capd::vectalg::midObject<DVector>(*chunk);
//			break;
//		}
//	DMatrix gridInvCoords = capd::vectalg::midObject<DMatrix>(grid3.get_I_invB());
//	for(auto chunk = gridSet.begin(); chunk != gridSet.end(); ++chunk)
//	{
//		DVector mid_chunk = capd::vectalg::midObject<DVector>(*chunk);
//		DVector diff = gridInvCoords*(mid_chunk - ref);
//		cout << diff << endl;
//		DVector v = x0 + diff[0] * coords_finite[1] + diff[1] * coords_finite[2];
//		the_set_chunks.push_back(v);
//	}
//	splotMany(roessler525, the_set_chunks, "the_set_chunks");
//
//
//
//	DMatrix infinite_coords;
}

void computeCoordsForward(
		system3d& system, DVector const& ivp,
		HSet2D& hset, int CUTS_Y, int CUTS_Z,
		std::vector<DVector> fin_mid_others){

	int& p = system.p;
	int& order = system.order;
	// const int& d = system.d;
	int d = 3; // .... argh....
	DGrid &grid = system.dgrid;

	DEq rhs(system.a.leftBound(), 0.2, 0.);	// eps = 0.
	DDDEq vf(rhs, grid(p));
	DSection section(3, 0, 0.);
	DSolver solver(vf, order * 3);
	DPoincare P(solver, section, poincare::MinusPlus);

	DSolution constIVP(grid, -grid(p), grid(0), order, ivp);
	DSolution Px(grid, -grid(p), grid(0), order, ivp);
	DSolution solution = constIVP;

	std::vector<DVector> attractor;
	double tp;
	int NUM_ITERS = 100;
	for (int i = 0; i < NUM_ITERS; ++i){
		P(solution, Px, tp);
		DVector pv = Px;
		pv[0] = 0.0;
		attractor.push_back(pv);
	}
	splotMany(system, attractor, "forward-history");

	std::string plotpath = "plots/";
	std::string filename = "forward-history";
	capd::ddeshelper::plot_value(plotpath + filename, solution.pastTime(), solution.currentTime(), system.h.leftBound() / 10, solution, false);
	std::ostringstream plot;
	plot << "'" << filename << "ddes-plot.dat' u 3:5:7 with lines notitle";
	capd::ddeshelper::splot_many(plotpath, { plot.str() }, false, "trajectory");
	std::ostringstream cmd; cmd << "cd '" << plotpath << "' && gnuplot 'trajectory.gp'";
	capd::ddeshelper::runSystemCommand(cmd.str());

	DVector zerov(system.M());
	auto mean = [zerov](std::vector<DVector> items){
		auto vv = std::accumulate(items.begin(), items.end(), zerov);
		vv /= double(items.size());
		return vv;
	};

	DMatrix fin_C = hset.get_B();
	DVector fin_mid_main = expand(hset.get_x());

	DVector fin_ydir = expand(fin_C.column(0)); // fin_C is 2 dim...
	DVector fin_zdir = expand(fin_C.column(1));
	cout << fin_ydir << endl;
	cout << fin_zdir << endl;


	std::sort(attractor.begin(), attractor.end(), [](DVector const& a, DVector const& b){
		return a[1] < b[1];
	});

	// this seems like bad idea to compute average direction
	// istead, i just connect two farthest points and this seems ok-ish
	cout << "Computing ydir" << endl;
	DVector ydir = attractor.front() - attractor.back();
	double norm = ydir[1]; ydir /= norm;
	ydir = -ydir; // żeby bylo jak u ani

	// mean vector seems ok-ish
	cout << "Computing mid point" << endl;
	std::vector<DVector> mid_candidates;
	for (int i = 0; i < attractor.size() / 3; ++i)
		for (int j = attractor.size() -1; j >= 2 * attractor.size() / 3; --j)
			mid_candidates.push_back(mean({attractor[i], attractor[j]}));

	capd::vectalg::SumNorm<DVector, DMatrix> selnorm;
	std::sort(mid_candidates.begin(), mid_candidates.end(), [&fin_mid_main, &selnorm](DVector const& a, DVector const& b){
		DVector aa{a[1], a[2]};
		DVector bb{b[1], b[2]};
		DVector cc{fin_mid_main[1], fin_mid_main[2]};

		return selnorm(aa-cc) < selnorm(bb-cc);
	});
	DVector mid_point = mid_candidates[0];

	auto translate = [&](DVector const& fin_v, DVector const& full_v){
		DSolution const_fin_v(grid, -grid(p), grid(0), order, fin_v);
		DSolution const_proj_v(grid, -grid(p), grid(0), order, {full_v[0], full_v[1], full_v[2]});
		return full_v - DVector(const_proj_v) + DVector(const_fin_v);
	};

	auto compheads = [&](DVector const& fin_v, DVector const& full_v){
		DVector proj_v {full_v[0], full_v[1], full_v[2]};
		cout << "old: " << fin_v << endl;
		cout << "new: " << proj_v << endl;
		cout << "dif: " << proj_v - fin_v << endl;
	};

	cout << "midpoint:" << endl; compheads(fin_mid_main, mid_point);
	cout << "ydir:    " << endl; compheads(fin_ydir, ydir);

	// I DO a small correction to the coordinates (they are very close the new ones and the old ones)

	// this seems like bad idea:
	// (it moves past points away from their position too much)
//	ydir = translate(fin_ydir, ydir);
//	mid_point = translate(fin_mid_main, mid_point);

	cout << "translated midpoint:" << endl; compheads(fin_mid_main, mid_point);
	cout << "translated ydir:    " << endl; compheads(fin_ydir, ydir);

	DVector zdir(system.M());
	// zdir[0] = fin_zdir[0]; // this should be 0
	zdir[1] = fin_zdir[1];
	zdir[2] = fin_zdir[2];

	// const solutions is not great, as it does not scale in the range over z
	// and it does not contain the derivatives
//	// TODO: rethink if keep using this
//	DSolution const_fin_zdir(grid, -grid(p), grid(0), order, fin_zdir);
//	zdir = DVector(const_fin_zdir);

	// an idea: put zdir * e^lambda * -t, for tin range[-tau, 0] as
	// the vector. The lambda should be choosen to scale well with the attractor thickness

	// we need to do better...
//	double LAMBDA = 3; // ?? TODO: choose
//	auto h = system.h.leftBound();
//	for (int i = 0; i < p; ++i){
//		double elambdat = exp(LAMBDA * i * h);
//		cout << elambdat << endl;
//		for (int k = 0; k <= order; ++k){
//			if (i == 0 && k > 0) continue;
//			int bi = (i == 0 ? 0 : d * (1 + (i-1) * (order+1) + k)); // baseindex
//			DVector v = fin_zdir * elambdat;
//			for (int j = 0; j < d; ++j)
//				zdir[bi + j] = v[j];
//			elambdat *= -LAMBDA / (k+1); // (e^LAMBDA*-t)^(k) / k! (dla k+1)
//		}
//	}
	// just leave it 2 dim and check if the deviation is big or not....

	cout << "Saving coordinates" << endl;
	// save coordinates to be used later
	IVector I_x0(mid_point);
	IMatrix I_C(system.M(), system.M());
	I_C.setToIdentity();
	I_C.column(1) = IVector(ydir);
	I_C.column(2) = IVector(zdir);
	IVector I_r0(system.M()); I_r0 *= 0.;

	std::ofstream human_readable_coords("human_readable_coords.txt");
	human_readable_coords << "x0: " << endl;
	human_readable_coords << mid_point << endl;
	human_readable_coords << "C: " << endl;
	human_readable_coords << capd::vectalg::midObject<DMatrix>(I_C) << endl;
	human_readable_coords << "r0: " << endl;
	human_readable_coords << I_r0 << endl;

	capd::ddeshelper::saveBinary("grid3_x0.ivector.bin", I_x0);
	capd::ddeshelper::saveBinary("grid3_C.imatrix.bin", I_C);
	capd::ddeshelper::saveBinary("grid3_r0.ivector.bin", I_r0);
	std::cout << "Please run from 'bin' folder: " << endl;
	std::cout << "../../../../bin/convmatrix ";
	std::cout << system.M() << " bin grid3_C.imatrix.bin ";
	std::cout << "inv bin " << "grid3_invC.imatrix.bin" << endl;

	GridSet gridSet(2);
	hset.gridSet(gridSet, CUTS_Y, CUTS_Z);
	IVector box = expand(gridSet.box());
	std::vector<DVector> the_set_chunks;
	DMatrix gridInvCoords = capd::vectalg::midObject<DMatrix>(hset.get_I_invB());
	for(auto chunk = gridSet.begin(); chunk != gridSet.end(); ++chunk)
	{
		DVector mid_chunk = capd::vectalg::midObject<DVector>(*chunk);
		DVector diff = gridInvCoords*(mid_chunk - hset.get_x());
		//DVector diff = 0.5 * gridInvCoords*(mid_chunk - hset.get_x()); // TODO: why i need 0.5 to get the right picture??? - rethink!
		//cout << diff << endl;
		DVector v = mid_point + diff[0] * ydir + diff[1] * zdir;
		//DVector v = mid_point + diff[0] * ydir;
		the_set_chunks.push_back(v);
	}
	splotMany(system, the_set_chunks, "the_set_chunks_forward");

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

		std::string plotpath = "plots/scatter/";
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


//		ofstream gp(plotpath + gpprefix.str() + ".gp");
//		gp <<


	}
}

