//////////////////////////////////////////////////////////////////////////////
///
///  @file utils.cpp
///  
///  @author agrs  @date   Apr 18, 2024
//////////////////////////////////////////////////////////////////////////////

// Standard libraries and the CAPD library
#include <iostream>
#include "utils.h"
#include "capd/covrel/HSet2D.h"
#include "capd/dynsys/DynSysMap.h"
#include "capd/covrel/HSetMD.h"
#include "capd/intervals/lib.h"
#include "capd/capdlib.h"

using namespace capd;
using namespace capd::alglib;
using namespace capd::matrixAlgorithms;
using namespace std;

typedef capd::covrel::HSet2D<DMatrix, IMatrix> HSet2D;
typedef capd::dynsys::DynSysMap<IMap> DynSysMap;
typedef capd::covrel::HSetMD<DMatrix,IMatrix> HSet;
typedef capd::covrel::GridSet<IMatrix> GridSet;

// this forces compiler to generate code for all non-template members of the classes
// it might make compilation faster in some cases.
template class Rossler<interval>;
template class capd::ddeshelper::RigorousHelper<Eq, 1, IMatrix, IVector>;

//////////////////////////////////////////////////////////////////////////////
///	Some auxiliary functions
//////////////////////////////////////////////////////////////////////////////

/*
 *  Just plots some nice pictures from the non-rigorous computations.
 */
void splotMany(
		system3d& sys,
		std::vector<DVector>& segments,
		std::string const& name,
		std::string const& wd){
	// plot all the solutions in the Ambient space
	const std::string plotpath = wd + name + "/";
	if (wd != "./" && wd != "") capd::ddeshelper::mkdir_p(plotpath);
	std::vector<std::string> gplots;
	std::cerr << "Plotting " << segments.size() << " segments... " << std::flush;
	for (int i = 0; i < segments.size(); ++i){
		DVector& v = segments[i];
		auto seg = sys.makeDSegment(v);
		std::ostringstream filename; filename << "seg-" << i;
		std::string fileprefix = plotpath + filename.str();
		capd::ddeshelper::plot_value(fileprefix, seg.pastTime(), seg.currentTime(), sys.h.leftBound() / 10, seg, false);
		std::ostringstream plot;
		plot << "'" << filename.str() << "ddes-plot.dat' u 3:5:7 with lines notitle";
		gplots.push_back(plot.str());
	}
	std::cerr << "Plotting DONE" << std::endl;

	capd::ddeshelper::splot_many(plotpath, gplots, false, "all");
	#ifdef DDES_ALLOW_SYSTEM
	// TODO: we need to force plot, here as xplot is usually precompiled without the flag, unfortunatelly
	std::ostringstream cmd; cmd << "cd '" << plotpath << "' && gnuplot 'all.gp'";
	capd::ddeshelper::runSystemCommand(cmd.str());
	#endif
}

/** define this function for the fixed types in template arguments. See utils.h for the template */
IMatrix cut(IMatrix const& M){
	if(M.numberOfColumns()<2 or M.numberOfRows() < 2) throw std::logic_error("cut(matrix): too few dimensions!");
	if(M.numberOfColumns()==2 and M.numberOfRows()==2) return M;
	IMatrix y(2,2);
	for (int i = 0; i < 2; ++i)
		for (int j = 0; j < 2; ++j)
			y[i][j]=M[i+1][j+1];
	return y;
}

/** define this function for the fixed types in template arguments. See utils.h for the template */
IMatrix expand(IMatrix const& M)
{		
	if(M.numberOfColumns()==3 and M.numberOfRows()==3) return M;
	if(M.numberOfColumns()!=2 or M.numberOfRows()!=2) throw std::logic_error("expand(matrix): bad dimensions");
	IMatrix y(3,3);
	y[1][1]=M[0][0];
	y[1][2]=M[0][1];
	y[2][1]=M[1][0];
	y[2][2]=M[1][1];
	y[0][1]=y[1][0]=y[2][0]=y[0][2]=0.;
	y[0][0]=1.;
	return y;
}

//////////////////////////////////////////////////////////////////////////////
///	system3d structure methods
//////////////////////////////////////////////////////////////////////////////

/*
 * checks if the image of hset1 lies inside hset2
 * TODO: better docs.
 */
bool system3d::estimate_piece(const HSet2D &hset1, const HSet2D &hset2, int howManyPiecesH, int howManyPiecesV, int iy, int iz)
{
	// TODO: spraw, aby brał pod uwagę hset1 i hset2 !!!! WAZNE!
	// TODO: 1. wczytaj opis seta ogolnego
	IMatrix C(M(), M());
	IMatrix invC(M(), M());
	IVector x0(M());
	IVector r0(M());
	// TODO: get data from outside! as inside_piece
	capd::ddeshelper::readBinary("grid3_x0.ivector.bin", x0);
	capd::ddeshelper::readBinary("grid3_C.imatrix.bin", C);
	capd::ddeshelper::readBinary("grid3_r0.ivector.bin", r0);
	capd::ddeshelper::readBinary("grid3_invC.imatrix.bin", invC);
	// we assume x0 is generated from the mid point of hset2...
	// if not, we need to adjust... TODO: dorobić...

	IVector chunk_box = r0;

	double ry = hset2.get_r()[0];
	double rz = hset2.get_r()[1];
	IVector corner = x0 - interval(ry) * C.column(1) - interval(rz) * C.column(2);
	auto dy = ry / howManyPiecesH;
	auto dz = rz / howManyPiecesV;

	r0[1] = interval(-ry, ry);
	r0[2] = interval(-rz, rz);
	chunk_box[1] = interval(-dy, dy);
	chunk_box[2] = interval(-dz, dz);

	interval rt(0.);

	// compute images
	IVector new_r0 = r0;
	IVector new_Xi(p * d);

	auto chunk_x0 = corner + C.column(1) * interval(dy + dy * 2 * iy) + C.column(1) * interval(dz + dz * 2 * iz);
	auto segment = makeISegment(chunk_x0);
	segment.set_Cr0(C, chunk_box);
	auto Psegment = P(segment, rt);
	IVector Phull = invC * (Psegment.get_x() - x0) +
				(invC * Psegment.get_C()) * Psegment.get_r0() +
				(invC * Psegment.get_B()) * Psegment.get_r();

	IVector resthull{0.,0.,0.};
	for (int j = d; j < Phull.dimension(); j+=3){
		IVector u {Phull[j+0], Phull[j+1], Phull[j+2]};
		resthull = capd::vectalg::intervalHull(resthull, u);
	}
	std::ostringstream oss;
	oss << boolalpha;
	oss << iy << " " << iz;
	for (int j = 1; j < d; j++){
		 oss << " " << Phull[j].subset(r0[j]); // << " " << resthull[j];
	}
	cout << oss.str() << endl;

	new_Xi = Psegment.get_Xi();;
	new_r0 = capd::vectalg::intervalHull(new_r0, Phull);
	new_r0[0] = 0.;

	string dirpath = "estimate_piece/";
	capd::ddeshelper::mkdir_p(dirpath);

	ostringstream prefix;
	prefix << "_" << iy << "_" << iz;

	std::ofstream human_readable_new_r0(dirpath + "human_readable_new_radius" + prefix.str() + ".txt");
	human_readable_new_r0 << "r0:" << endl << new_r0 << endl;
	human_readable_new_r0 << "Xi:" << endl << new_Xi << endl;
	human_readable_new_r0.close();

	capd::ddeshelper::saveBinary(dirpath + "grid3_new_r0" + prefix.str() + ".ivector.bin", new_r0);
	capd::ddeshelper::saveBinary(dirpath + "grid3_new_Xi0" + prefix.str() + ".ivector.bin", new_Xi);

	return false;
}

/*
 * Checks if the image of hset1 lies inside hset2
 *
 */
bool system3d::refine_box(const HSet2D &hset1, const HSet2D &hset2, int howManyPiecesH, int howManyPiecesV, int iteration)
{
	// TODO: spraw, aby brał pod uwagę hset1 i hset2 !!!! WAZNE!
	// TODO: 1. wczytaj opis seta ogolnego
	IMatrix C(M(), M());
	IMatrix invC(M(), M());
	IVector x0(M());
	IVector r0(M());
	capd::ddeshelper::readBinary("grid3_x0.ivector.bin", x0);
	capd::ddeshelper::readBinary("grid3_C.imatrix.bin", C);
	capd::ddeshelper::readBinary("grid3_r0.ivector.bin", r0);
	capd::ddeshelper::readBinary("grid3_invC.imatrix.bin", invC);
	// we assume x0 is generated from the mid point of hset2...
	// if not, we need to adjust... TODO: dorobić...

	IVector chunk_box = r0;

	double ry = hset2.get_r()[0];
	double rz = hset2.get_r()[1];
	IVector corner = x0 - interval(ry) * C.column(1) - interval(rz) * C.column(2);
	auto dy = ry / howManyPiecesH;
	auto dz = rz / howManyPiecesV;

	r0[1] = interval(-ry, ry);
	r0[2] = interval(-rz, rz);
	chunk_box[1] = interval(-dy, dy);
	chunk_box[2] = interval(-dz, dz);

	auto segment = makeISegment(x0);
	interval rt(0.);
	auto Psegment = P(segment, rt);
	cout << "Midpoint: return time: " << rt << "(sanity check)" << endl;

	// compute images
	IVector new_r0 = r0;
	IVector new_Xi(p * d);
	int count = 0;
//	for (int iy = 0; iy < howManyPiecesH; ++iy){
//		for (int iz = 0; iz < howManyPiecesV; ++iz){
	for (int iy = howManyPiecesH-1; iy >= 0; --iy){
		for (int iz = howManyPiecesV-1; iz >= 0; --iz){
			auto chunk_x0 = corner + C.column(1) * interval(dy + dy * 2 * iy) + C.column(1) * interval(dz + dz * 2 * iz);
			auto segment = makeISegment(chunk_x0);
			segment.set_Cr0(C, chunk_box);
			auto Psegment = P(segment, rt);
			cout << "Piece: " << count << " return time: " << rt << endl;
			// TODO: check condition:
			IVector Phull = invC * (Psegment.get_x() - x0) +
						(invC * Psegment.get_C()) * Psegment.get_r0() +
						(invC * Psegment.get_B()) * Psegment.get_r();

			IVector resthull{0.,0.,0.};
			for (int j = d; j < Phull.dimension(); j+=3){
				IVector u {Phull[j+0], Phull[j+1], Phull[j+2]};
				resthull = capd::vectalg::intervalHull(resthull, u);
			}
			for (int j = 1; j < d; j++){
				cout << "    main:   " << Phull[j] << " vs " << r0[j] << " ? " << Phull[j].subset(r0[j]) << endl;
				cout << "    others: " << resthull[j] << endl;
			}

			IVector Xi = Psegment.get_Xi();
//			cout << Xi << endl;
//			cout << Xi.dimension() << " " << new_Xi.dimension() << endl;
			if (count == 0) new_Xi = Xi;
			else new_Xi = capd::vectalg::intervalHull(new_Xi, Xi);

			new_r0 = capd::vectalg::intervalHull(new_r0, Phull);

			++count;
			//if (++count > 0) break; // TODO: just for debug...
		}
	}

	new_r0[0] = 0.;
	std::ofstream human_readable_new_r0("human_readable_new_radius.txt");
	human_readable_new_r0 << "r0:" << endl << new_r0 << endl;
	human_readable_new_r0 << "Xi:" << endl << new_Xi << endl;
	human_readable_new_r0.close();

	new_r0 *= interval(-1.1, 1.1); // make it a little bigger (v(0) coordinates will be overwritten)
	new_Xi *= interval(-1.1, 1.1);
	capd::ddeshelper::saveBinary("grid3_new_r0.ivector.bin", new_r0);
	capd::ddeshelper::saveBinary("grid3_new_Xi0.ivector.bin", new_Xi);

	return false;


//	GridSet gridSet(2);																		// stores grid
//	hset1.gridSet(gridSet,howManyPiecesH,howManyPiecesV);
//	HSet2D hset2Base(IVector({0.,0.}),IMatrix::Identity(2), hset2.get_r());
//
//	IVector Pset2d(2);
//	interval rt(0.);
//
//	IMatrix expandedGridSetCoordinateSystem = expand(gridSet.coordinateSystem());			// expanded to 3d, where P is defined
//	IVector expandedGridSetBox = expand(gridSet.box());
//	IMatrix expandedHSet2InvCoordinateSystem = expand(hset2.invCoordinateSystem());
//	IVector expandedHSet2Center = expand(hset2.center());
//
//	bool liesInside = 1;
//	P.setRequiredSteps(0);
//		auto zero = grid(0);
//		auto tau = grid(p);
//	for(auto i = gridSet.begin();i!=gridSet.end();++i)
//	{
//		DDESolution segment(grid, -tau, zero, order, AmbientSet(expand(*i), expandedGridSetCoordinateSystem, expandedGridSetBox));				// (grid, from, to, order, constant function value)
//		auto Psegment = P(segment, rt);
//
//		IVector Phull = invC * (segment.get_x() - x0) +
//					(invC * segment.get_C()) * segment.get_r0() +
//					(invC * segment.get_B()) * segment.get_r();
//
//		// rt = P.getLastReachTime();
//		cout << "rt = " << rt << endl;
//
//		//Pset2d = cut(PV);		// the image of the small box
//
//		//liesInside=liesInside && hset2Base.inside(Pset2d);									// do all lie inside?
//
//		if(!liesInside)
//			{
//				// cout << "Does not lie inside, for this tiny box: " << Pset2d << " does not lie inside the box " << hset2Base << endl;
//				cout << "Failed the check..." << endl;
//				return 0;
//			}
//
//		break; // TODO: JUST FOR NOW to check only one piece.
//	}
//    return liesInside;
}

// checks if the image of hset1 lies inside hset2
bool system3d::inside(const HSet2D &hset1, const HSet2D &hset2, int howManyPiecesH, int howManyPiecesV, int iteration)
{

	// TODO: spraw, aby brał pod uwagę hset1 i hset2 !!!! WAZNE!

	// TODO:

	IMatrix C(M(), M());
	IMatrix invC(M(), M());
	IVector x0(M());
	IVector r0(M());
	IVector Xi0(d * p);
	capd::ddeshelper::readBinary("grid3_x0.ivector.bin", x0);
	capd::ddeshelper::readBinary("grid3_C.imatrix.bin", C);
	capd::ddeshelper::readBinary("grid3_new_r0.ivector.bin", r0);
	capd::ddeshelper::readBinary("grid3_new_Xi.ivector.bin", Xi0);
	capd::ddeshelper::readBinary("grid3_invC.imatrix.bin", invC);
	// we assume x0 is generated from the mid point of hset2...
	// if not, we need to adjust... TODO: dorobić...

	IVector chunk_box = r0;

	double ry = hset2.get_r()[0];
	double rz = hset2.get_r()[1];
	IVector corner = x0 - interval(ry) * C.column(1) - interval(rz) * C.column(2);
	auto dy = ry / howManyPiecesH;
	auto dz = rz / howManyPiecesV;

	r0[1] = interval(-ry, ry);
	r0[2] = interval(-rz, rz);
	chunk_box[1] = interval(-dy, dy);
	chunk_box[2] = interval(-dz, dz);

	auto segment = makeISegment(x0);
	interval rt(0.);
	auto Psegment = P(segment, rt);
	cout << "Midpoint: return time: " << rt << "(sanity check)" << endl;

	// compute images
	IVector new_r0 = r0;
	IVector new_Xi = r0;
	int count = 0;
	bool liesInside = 1;
	for (int iy = 0; iy < howManyPiecesH; ++iy){
		for (int iz = 0; iz < howManyPiecesV; ++iz){
			auto chunk_x0 = corner + C.column(1) * interval(dy + dy * 2 * iy) + C.column(1) * interval(dz + dz * 2 * iz);
			auto segment = makeISegment(chunk_x0);
			segment.set_Cr0(C, chunk_box);
			segment.set_Xi(Xi0);
			auto Psegment = P(segment, rt);
			cout << "Piece: " << count << " return time: " << rt << endl;
			// TODO: check condition:
			IVector Phull = invC * (Psegment.get_x() - x0) +
						(invC * Psegment.get_C()) * Psegment.get_r0() +
						(invC * Psegment.get_B()) * Psegment.get_r();

			for (int j = 1; j < Phull.dimension(); j++){
				bool testj = Phull[j].subset(r0[j]);
				liesInside = liesInside && testj;
				if (!testj) cerr << "    BAD r at " << iy << " " << iz << " " << j << ": " << Phull[j] << " vs " << r0[j] << " ? " << testj << endl;
			}
			IVector PXi = Psegment.get_Xi();
			for (int j = 0; j < PXi.dimension(); j++){
				bool testXi = PXi[j].subset(Xi0[j]);
				liesInside = liesInside && testXi;
				if (!testXi) cerr << "    BAD Xi at " << iy << " " << iz << " " << j << ": " << PXi[j] << " vs " << Xi0[j] << " ? " << testXi << endl;
			}
			++count;
		}
	}

	return liesInside;
}

/*
 * Checks if the (piece of) image of hset1 lies inside hset2:
 *
 * i.e: return P(hset1) \subset hset2
 *
 * This is the simplest form of covering relation, with all stable directions.
 */
bool system3d::inside_piece(
		const HSet2D &hset1, IVector const& mid1,
		const HSet2D &hset2, IVector const& mid2,
		const IVector& r0, const IVector& Xi0,
		const IMatrix& C, const IMatrix& invC,
		int howManyPiecesH, int howManyPiecesV,
		int pieceH, int pieceV,
		IVector &outPhull, IVector &outPXi)
{
	// TODO: spraw, aby brał pod uwagę hset1 i hset2 !!!! WAZNE!
	// TODO: check if the above is done now

	// TODO: Before we had control on M, invM, etc.
	// TODO: Now, we get all from outside. Add checks for the dimensions of arrays and vectors?
	// TODO: It is not so important, as the exceptions will be thrown otherwise.

	IVector chunk_box = r0;

	double ry1 = hset1.get_r()[0];
	double rz1 = hset1.get_r()[1];
	IVector corner = mid1 - interval(ry1) * C.column(1) - interval(rz1) * C.column(2);
	auto dy = ry1 / howManyPiecesH;
	auto dz = rz1 / howManyPiecesV;

	IVector hset1_r0 = r0, hset2_r0 = r0;

	hset1_r0[1] = interval(-ry1, ry1);
	hset1_r0[2] = interval(-rz1, rz1);
	chunk_box[1] = interval(-dy, dy);
	chunk_box[2] = interval(-dz, dz);

	double ry2 = hset2.get_r()[0];
	double rz2 = hset2.get_r()[1];
	hset2_r0[1] = interval(-ry2, ry2);
	hset2_r0[2] = interval(-rz2, rz2);

	auto& iy = pieceH;
	auto& iz = pieceV;

	bool liesInside = 1;
	interval rt(0.);
	auto chunk_x0 = corner + C.column(1) * interval(dy + dy * 2 * iy) + C.column(1) * interval(dz + dz * 2 * iz);
	auto segment = makeISegment(chunk_x0);
	segment.set_Cr0(C, chunk_box);
	segment.set_Xi(Xi0);
	auto Psegment = P(segment, rt);
	cerr << "Piece: " << iy << " " << iz << " return time: " << rt << endl;
	IVector Phull = invC * (Psegment.get_x() - mid2) +   // here we use the other mid point (we check condition in hset2)
				(invC * Psegment.get_C()) * Psegment.get_r0() +
				(invC * Psegment.get_B()) * Psegment.get_r();

	for (int j = 1; j < Phull.dimension(); j++){
		bool testj = Phull[j].subset(hset2_r0[j]);
		liesInside = liesInside && testj;
		if (!testj) cerr << "    BAD r at " << iy << " " << iz << " " << j << ": " << Phull[j] << " vs " << hset2_r0[j] << " ? " << testj << endl;
	}
	IVector PXi = Psegment.get_Xi();
	for (int j = 0; j < PXi.dimension(); j++){
		bool testXi = PXi[j].subset(Xi0[j]);
		liesInside = liesInside && testXi;
		if (!testXi) cerr << "    BAD Xi at " << iy << " " << iz << " " << j << ": " << PXi[j] << " vs " << Xi0[j] << " ? " << testXi << endl;
	}

	outPhull = Phull;
	outPXi = PXi;
	return liesInside;
}
