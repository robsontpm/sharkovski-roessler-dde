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
		
		//const int CUTS_Y = 1000;
		//const int CUTS_Y = 100;
		const int CUTS_Y = 6; // TODO: for now smaller to faster test
		const int CUTS_Z = 3;
		//computeCoordsByC1(roessler525, grid3, CUTS_Y, CUTS_Z);
		computeCoordsSimple(roessler525, grid3, CUTS_Y, CUTS_Z);

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
	for (int i = 0; i < x0_all.size(); ++i){
		DVector& v = x0_all[i];
		auto seg = roessler525.makeDSegment(v);
		std::ostringstream filename; filename << "seg-" << i;
		std::string fileprefix = plotpath + filename.str();
		// cout << v << endl;
		cout << "plotting " << i << " " << fileprefix << endl;
		// capd::ddeshelper::plot_value(fileprefix.str(), seg, false);
		capd::ddeshelper::plot_value(fileprefix, seg.pastTime(), seg.currentTime(), roessler525.h.leftBound() / 10, seg, false);
		std::ostringstream plot;
		plot << "'" << filename.str() << "ddes-plot.dat' u 3:5:7 with lines notitle";
		gplots.push_back(plot.str());
	}

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

