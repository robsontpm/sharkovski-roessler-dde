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
