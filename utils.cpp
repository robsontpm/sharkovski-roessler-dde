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
	// For Devs: here we need to force plot, here as xplot is usually precompiled without the flag, unfortunately
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
 *
 * NOTE:
 * this is specifically tailored to handle data from the accompanying paper.
 * The mid_i points are M dimensional vectors
 * (M big, dimension of the representation of functions in C^n_p, see paper)
 * The hset_i are 2 dimensional h-sets from work [GZ2022], that will define
 * the head \pi_{y,z} z(C_i), C_i \in C^n_p.
 * r0, Xi0 are to define thickness of the sets in C^n_p,
 * but such that x(r0) = [0,0] (this is the head, it will be taken from the hset_i)
 * 0 is x = 0 (we are on section).
 * C and invC = C^{-1} are the coordinate change (c_N in definition of hset)
 * C \in Matrix(M, M) (big matrix). C^{-1} is rigorous inverse, i.e. Id \in C * C^{-1}.
 * The sets C_i are defined as:
 * 		C_i = { 0 } x mid(hset_i) x \pi_{j >= 3} mid_i + C \cdot (({radius(hset_i)}, 0 ... ) + r0)
 * (this is ugly, see the accompanying papaer for more explanation)
 *
 * Note: radius(hset_i) := hset_i - mid(hset_i)
 *
 * Note:
 * the sets hset_1, hset_2 will be mathematically speaking hset_{i_1}, hset_{i_2}
 * and we will choose i_1, i_2 as needed to check the chain of covering relations
 * between many sets.
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
	// sanity check on input
	auto check_dim_vector = [&](IVector const& v, std::string const& name, int dim=-1) {
		if (dim < 0) dim = M();
		if (v.dimension() != dim){
			std::ostringstream info;
			info << "system3d::inside_piece(): vector '" << name << "' has bad dimension, ";
			info << "is: " << v.dimension() << ", should be: " << dim;
			throw std::logic_error(info.str());
		}
	};
	auto check_dim_matrix = [&](IMatrix const& C, std::string const& name) {
		if (C.numberOfColumns() != M() || C.numberOfRows() != M()){
			std::ostringstream info;
			info << "system3d::inside_piece(): vector '" << name << "' has bad dimension, ";
			info << "is: " << C.numberOfRows() << "x" << C.numberOfColumns() << ", should be: " << M() << "x" << M();
			throw std::logic_error(info.str());
		}
	};
	check_dim_vector(r0, "r0");
	check_dim_vector(mid1, "mid1");
	check_dim_vector(mid2, "mid2");
	check_dim_vector(Xi0, "Xi", p * d);
	check_dim_matrix(C, "M");
	check_dim_matrix(invC, "M^{-1}");
	// end of sanity check

	// we divide the head z(\cdot) of the h-set with tails on y and z coordinates
	// into uniform grid,  with the given number of howManyPieces and we select one chunk
	// of the whole thing to do the calculations
	IVector chunk_box = r0; // this will hold the diameter of the chunk of the set
	double ry1 = hset1.get_r()[0];	// hset is 2 dim on section x = 0,
	double rz1 = hset1.get_r()[1];  // so y is on index 0, z is on index 1
	IVector corner = mid1 - interval(ry1) * C.column(1) - interval(rz1) * C.column(2);
	auto dy = ry1 / howManyPiecesH;
	auto dz = rz1 / howManyPiecesV;

	// we fill all directions except the head
	IVector hset1_r0 = r0, hset2_r0 = r0;

	// we fill head with a chunk
	hset1_r0[1] = interval(-ry1, ry1);
	hset1_r0[2] = interval(-rz1, rz1);
	chunk_box[1] = interval(-dy, dy);
	chunk_box[2] = interval(-dz, dz);

	// this is the whole set on the head
	double ry2 = hset2.get_r()[0];
	double rz2 = hset2.get_r()[1];
	hset2_r0[1] = interval(-ry2, ry2);
	hset2_r0[2] = interval(-rz2, rz2);

	// renaming
	auto& iy = pieceH;
	auto& iz = pieceV;

	// this will stay true when the condition "P(C_{i_1}) lies inside C_{i_2}" is satisfied
	bool liesInside = 1;
	interval rt(0.);
	// compute the middle point of the chunk of hset_1
	auto chunk_x0 = corner + C.column(1) * interval(dy + dy * 2 * iy) + C.column(1) * interval(dz + dz * 2 * iz);
	// we configure the chunk to be a representation of a solution segment
	auto segment = makeISegment(chunk_x0); // first, we have point set (no radius)
	segment.set_Cr0(C, chunk_box); 	       // then we set the radius
	segment.set_Xi(Xi0);				   // and the estimates on the higher derivatives
	auto Psegment = P(segment, rt);		   // this computes the image and make sure the image is on the section.
	cerr << "Piece: " << iy << " " << iz << " return time: " << rt << endl;  // debug information
	IVector Phull = invC * (Psegment.get_x() - mid2) +   			// we transform the set into good coordinate frame
				(invC * Psegment.get_C()) * Psegment.get_r0() +		// of the set hset_2, both have the same c_N = C,
				(invC * Psegment.get_B()) * Psegment.get_r();		// but the mid point is different, so  we use mid_2

	// we check the lies inside on all coordinates except x (it is 0 on section)
	for (int j = 1; j < Phull.dimension(); j++){
		bool testj = Phull[j].subsetInterior(hset2_r0[j]);
		liesInside = liesInside && testj;
		if (!testj) cerr << "    BAD r at " << iy << " " << iz << " " << j << ": " << Phull[j] << " vs " << hset2_r0[j] << " ? " << testj << endl;
	}
	// we check the lies inside the higher derivatives
	// we can use here just a subset (not subsetInterior) as this is tail part of the hset with tail.
	// nevertheless, the computations are such that equality never happens.
	// see paper for details.
	IVector PXi = Psegment.get_Xi();
	for (int j = 0; j < PXi.dimension(); j++){
		bool testXi = PXi[j].subset(Xi0[j]);
		liesInside = liesInside && testXi;
		if (!testXi) cerr << "    BAD Xi at " << iy << " " << iz << " " << j << ": " << PXi[j] << " vs " << Xi0[j] << " ? " << testXi << endl;
	}

	// we return the images (in good coordinates) to use them in other computations
	// we reuse this code i.e. when finding good candidate for r0 and Xi.
	outPhull = Phull;
	outPXi = PXi;
	return liesInside;
}
