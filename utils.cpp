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

template class Rossler<interval>;
template class capd::ddeshelper::RigorousHelper<Eq, 1, IMatrix, IVector>;

//////////////////////////////////////////////////////////////////////////////
///	Some auxiliary functions
//////////////////////////////////////////////////////////////////////////////

// cut
// cuts the first coordinate of an interval 3-vector
IVector cut(IVector x)						
{
	if(x.dimension()==2) return x;
	IVector y(2);
	y[0]=x[1];
	y[1]=x[2];
	return y;
}

DVector cut(DVector x)						
{
	if(x.dimension()==2) return x;
	DVector y(2);
	y[0]=x[1];
	y[1]=x[2];
	return y;
}

// cuts an interval 3x3-matrix to 2x2
IMatrix cut(IMatrix M)						
{
	if(M.numberOfColumns()==2 and M.numberOfRows()==2) return M;
	IMatrix y(2,2);
	y[0][0]=M[1][1];
	y[0][1]=M[1][2];
	y[1][0]=M[2][1];
	y[1][1]=M[2][2];
	return y;
}

//expand
//prepends an interval 2-vector to 3-vector with zero 
IVector expand(IVector x)					
{	
	if(x.dimension()==3) return x;										
	IVector y({0.,0.,0.});
	y[1]=x[0];
	y[2]=x[1];
	return y;
}
DVector expand(DVector x)					
{	
	if(x.dimension()==3) return x;										
	DVector y({0.,0.,0.});
	y[1]=x[0];
	y[2]=x[1];
	return y;
}

//expands an interval 2x2 matrix to 3x3 with identity
IMatrix expand(IMatrix M)					
{		
	if(M.numberOfColumns()==3 and M.numberOfRows()==3) return M;									
	IMatrix y(3,3);
	y[1][1]=M[0][0];
	y[1][2]=M[0][1];
	y[2][1]=M[1][0];
	y[2][2]=M[1][1];
	y[0][1]=y[1][0]=y[2][0]=y[0][2]=interval(0.);
	y[0][0]=interval(1.);
	return y;
}

/////////////////////////////////////////////////////////////////////////
/// SecMap class
/////////////////////////////////////////////////////////////////////////

// image of the iteration iter only, no derivative calculation (faster)
IVector SecMap::image(const IVector &x, const int iter) const		
{
	IVector X({0.,0.,0.});
	X[1]=x[0];
	X[2]=x[1];
	
	//C0Rect2Set Q(X);
	DDESolution Q(*grid, (*grid)(0), X); 
	
	IVector Y=(*P)(Q,iter);
	
	IVector y(2);
	y[0]=Y[1];
	y[1]=Y[2];
	
	return y;
}


//////////////////////////////////////////////////////////////////////////////
///	system3d structure methods
//////////////////////////////////////////////////////////////////////////////

// checks if the image of hset1 lies inside hset2
bool system3d::estimate_piece(const HSet2D &hset1, const HSet2D &hset2, int howManyPiecesH, int howManyPiecesV, int iy, int iz)
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

// checks if the image of hset1 lies inside hset2
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

// checks if the (piece of) image of hset1 lies inside hset2
bool system3d::inside_piece(
		const HSet2D &hset1, IVector const& mid1,
		const HSet2D &hset2, IVector const& mid2,
		int howManyPiecesH, int howManyPiecesV,
		int pieceH, int pieceV,
		IVector &outPhull, IVector &outPXi)
{

	// TODO: spraw, aby brał pod uwagę hset1 i hset2 !!!! WAZNE!

	IMatrix C(M(), M());
	IMatrix invC(M(), M());
	//IVector x0(M());
	IVector r0(M());
	IVector Xi0(d * p);
	//capd::ddeshelper::readBinary("grid3_x0.ivector.bin", x0);
	//x0 = mid1;
	capd::ddeshelper::readBinary("grid3_C.imatrix.bin", C);
	capd::ddeshelper::readBinary("grid3_new_r0.ivector.bin", r0);
	capd::ddeshelper::readBinary("grid3_new_Xi.ivector.bin", Xi0);
	capd::ddeshelper::readBinary("grid3_invC.imatrix.bin", invC);
	// we assume x0 is generated from the mid point of hset2...
	// if not, we need to adjust... TODO: dorobić...

	IVector chunk_box = r0;

	double ry1 = hset1.get_r()[0];
	double rz1 = hset1.get_r()[1];
	IVector corner = /*x0*/ mid1 - interval(ry1) * C.column(1) - interval(rz1) * C.column(2);
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

void system3d::makeHistory(const HSet2D& hset){
	IMap backvf("par:a,b; var:x,y,z; fun: y + z, -x - b*y, -b - z*(x - a);");
	interval bb = interval(0.2);
	backvf.setParameter("a",a);
	backvf.setParameter("b",bb);
	IOdeSolver backsolver(backvf, 40);
	//odesolver.turnOffStepControl();
	//odesolver.setStep(-h);
	ITimeMap backmap(backsolver);
	C0Rect2Set set3D(expand(hset.center()), expand(hset.coordinateSystem()), expand(hset.box()));
	ofstream points("history.txt");
	std::vector<C0Rect2Set> history;
	for (int i = 0; i <= p; ++i){
//		odesolver.setStep(-h);
//		odesolver.turnOffStepControl();
		backmap(h * i, set3D);
		history.push_back(set3D);
		IVector hull = set3D, mid, r0;
		split(hull, mid, r0);

		cout << hull << endl;
		points << mid[0].leftBound() << " " << mid[1].leftBound() << " " << mid[2].leftBound() << " ";
		points << r0[0].rightBound() << " " << r0[1].rightBound() << " " << r0[2].rightBound() << endl;
	}
	points.close();
	ofstream gp("plot.gp");
	gp << "set terminal png size 800,600" << endl;
	gp << "set output 'xy.png'" << endl;
	gp << "plot 'history.txt' using 1:2:4:5 with boxxyerrorbars" << endl;
	gp << "set output 'xz.png'" << endl;
	gp << "plot 'history.txt' using 1:3:4:6 with boxxyerrorbars" << endl;
	gp << "set output 'yz.png'" << endl;
	gp << "plot 'history.txt' using 2:3:5:6 with boxxyerrorbars" << endl;
	gp.close();

	{ auto dump = system("gnuplot plot.gp"); }

	// let's try now to integrate last set forward in time and check how it
	// is compared to other sets on the backward trajectory.

	IMap forwardvf("par:a,b; var:x,y,z; fun: -y - z, x + b*y, b + z*(x - a);");
	forwardvf.setParameter("a",a);
	forwardvf.setParameter("b",bb);
	IOdeSolver forwardsolver(forwardvf, 40);
	ITimeMap forwardmap(forwardsolver);
	IMatrix Cinv = capd::matrixAlgorithms::inverseMatrix(set3D.get_C());
	IVector error_correction({0, 0, 0});
	// error_correction = (Cinv * set3D.get_B()) * set3D.get_r();
	C0Rect2Set testSet(set3D.get_x(), set3D.get_C(), set3D.get_r0() + error_correction); // to reset time
	for (int i = 1; i < p; ++i){

		forwardmap(h * i, testSet);
		auto histSet = history[p-i]; // we need to compare set with this one.
		IMatrix M = histSet.get_C();
		IMatrix Minv = capd::matrixAlgorithms::inverseMatrix(histSet.get_C());

		IVector hullInC = testSet.affineTransformation(Minv, histSet.get_x());
		cout << "forward " << i << endl;
		cout << histSet.get_r0() << endl;
		cout << hullInC << endl << endl;
	}

	// z historii wartości wylicz (p,n)-reprezentację. Na razie liczymy dla n = 0, 1, i 2

	// Df(v) = {
	// 0          -1        -1
	// 1           b        0
	// v[2]        0        v[0]
	// }
	//
	// niech będzie:
	//   v \in v0 + C * r_0   (to znamy z historii, tzn v(-h*i) =(jako zbiór)= history[i], i \in 0, ..., p;
	// to
	//   v' \in f(v0 + C * r_0) \approx f(v0) + (Df(v0 + C*r_0) * C) * r_0 \approx (usuwamy przedział wewnatrz DF)
	//                                f(v0) + (Df(v0) * C) * r_0 =: v0' + C' * r_0
	// teraz v'' spełnia
	//   v'' = Df(v) \circ Dv = Df(v) \cdot v' \approx Df(v) \cdot (v0' + C' * r_0) \approx
	//         Df(v0) v0'   +   (Df(v0) C') r_0   =:   v0'' + C'' r_0
	//
	// więc v^[2] = 1/2 (powyżej)
	// i tak można to pociągnąc rekurencyjnie wyżej, jakbysmy chcieli.

}

void system3d::makeHistoryC0(const DVector& v0, DVector& out_x0){
	DMap backvf("par:a,b; var:x,y,z; fun: y + z, -x - b*y, -b - z*(x - a);");
	interval bb = interval(0.2);
	double b = bb.leftBound();
	backvf.setParameter("a", a.leftBound());
	backvf.setParameter("b", bb.leftBound());
	DOdeSolver backsolver(backvf, 40);
	DTimeMap backmap(backsolver);
	HistoryType history;
	auto V = v0;
	double steph = h.leftBound();

	// value at t = 0
	for (int j = 0; j < d; ++j) out_x0[j] = v0[j];
	// propagate back
	for (int i = 0; i < p; ++i){
		auto newV = backmap(steph, V);
		V = newV;

		// signs are done for the backward iteration
		DMatrix Dfv({
			{0,      1,   1   },
			{-1,    -b,   0   },
			{-V[2],  0., -V[3]}
		});
		DVector VV[order+1]; // VV[k] = V^(k)
		VV[0] = V;
		VV[1] = -backvf(VV[0]); // -, bo potrzebuję do przodu! a vf jest do tyłu

		VV[2] = -Dfv * VV[1];

		VV[3] = -Dfv * VV[2];
		VV[3][2] += 2.*(VV[1][2]*VV[1][0]); // -2. * z'x'
		// TODO: compute regardless of order (now up to 3 (included))

		// make them coeffs (i.e. v^(k)/k!
		int factorial_k = 1;
		for (int k = 0; k <= order; k++){
			VV[k] *= 1./factorial_k;
			factorial_k *= (k + 1);
		}

		// I can simply fill the values of the vectors and the middle
		int base = d*(1 + (order+1)*i);
		for (int k = 0; k <= order; ++k){ // order
			for (int j = 0; j < d; ++j){ // ambient dimension
				auto row = base + k * d + j;
				out_x0[row] = VV[k][j];
			}
		}
	}
}

typename system3d::HistoryType
system3d::makeHistoryD(const DVector& v0, const DMatrix& C, DVector& out_x0, std::vector<DVector>& out_coords){
	DMap backvf("par:a,b; var:x,y,z; fun: y + z, -x - b*y, -b - z*(x - a);");
	// DMap forwardvf("par:a,b; var:x,y,z; fun: -y - z, x + b*y, b + z*(x - a);"); // TODO może niepotrzebne
	interval bb = interval(0.2);
	double b = bb.leftBound();
	backvf.setParameter("a", a.leftBound());
	backvf.setParameter("b", bb.leftBound());
	DOdeSolver backsolver(backvf, 40);
	DTimeMap backmap(backsolver);
	HistoryType history;
	auto V = v0;
	auto DV = C;
	double steph = h.leftBound();
	//bool NORMALIZE = true;
	bool NORMALIZE = false;
	SimpleJet jet;
	jet.push_back(std::make_pair(V, C)); // only value at t=0
	history.push_back(jet);

	// initialize coords at t=0
	for (int j = 0; j < d; ++j){ // ambient dimension
		out_x0[j] = V[j];
		for (int col = 0; col < d; ++col)
			out_coords[col][j] = C[j][col];
	}

	for (int i = 0; i < p; ++i){
		DMatrix newDV(d, d);
		auto newV = backmap(steph, V, DV, newDV);
		if (NORMALIZE){
			for (int j = 0; j < d; ++j)
				newDV.column(j).normalize();
		}
		DV = newDV;
		V = newV;

		// signs are done for the backward iteration
		DMatrix Dfv({
			{0,      1,   1   },
			{-1,    -b,   0   },
			{-V[2],  0., -V[3]}
		});
		DVector VV[order+1]; // VV[k] = V^(k)
		DMatrix CC[order+1];
		VV[0] = V;
		CC[0] = DV;
		VV[1] = -backvf(VV[0]); // -, bo potrzebuję do przodu! a vf jest do tyłu
		CC[1] = -Dfv * CC[0];

		VV[2] = -Dfv * VV[1];
		CC[2] = -Dfv * CC[1];
		// correct v''' representation with nonlinear terms (on z coordinate)
		// -(x' * C_z + z' * C_x)
		CC[2].row(2) += -(VV[1][0] * CC[0].row(2) + VV[1][2] * CC[0].row(0));

		VV[3] = -Dfv * VV[2];
		VV[3][2] += -2.*(VV[1][2]*VV[1][0]); // -2. * z'x'
		CC[3] = -Dfv * CC[2];      // 'linear' part Df(v) * v''
		// correct v''' representation with nonlinear terms (on z coordinate)
		// -2*(x' * C'_z + z' * C'_x)
		CC[3].row(2) += -2. * (VV[1][2] * CC[1].row(0) + VV[1][0] * CC[1].row(2));
		// TODO: compute regardless of order (now up to 3 (included))

//		// TODO: just check...
//		CC[1] *= 0.;
//		CC[2] *= 0.;
//		CC[3] *= 0.;

		// make them coeffs (i.e. v^(k)/k!
//		int factorial_k = 1;
//		for (int k = 0; k <= order; k++){
//			VV[k] *= 1./factorial_k;
//			CC[k] *= 1./factorial_k;
//			factorial_k *= (k + 1);
//		}

		// save full jet at the given time
		SimpleJet jet;
		for (int k = 0; k <= order; ++k)
			jet.push_back(std::make_pair(VV[k], CC[k]));
		history.push_back(jet);

		// I can simply fill the values of the vectors and the middle
		int base = d*(1 + (order+1)*i);
		for (int k = 0; k <= order; ++k){ // order
			for (int j = 0; j < d; ++j){ // ambient dimension
				auto row = base + k * d + j;
				out_x0[row] = VV[k][j];
				for (int col = 0; col < d; ++col)
					out_coords[col][row] = CC[k][j][col];
			}
		}
	}

	return history;

	// z historii wartości wylicz (p,n)-reprezentację. Na razie liczymy dla n = 0, 1, i 2

	// Df(v) = {
	// 0          -1        -1
	// 1           b        0
	// v[2]        0        v[0]
	// }
	//
	// niech będzie:
	//   v \in v0 + C * r_0   (to znamy z historii, tzn v(-h*i) =(jako zbiór)= history[i], i \in 0, ..., p;
	// to
	//   v' \in f(v0 + C * r_0) \approx f(v0) + (Df(v0 + C*r_0) * C) * r_0 \approx (usuwamy przedział wewnatrz DF)
	//                                f(v0) + (Df(v0) * C) * r_0 =: v0' + C' * r_0
	// teraz v'' spełnia
	//
	//   v'' = Df(v) \circ Dv = Df(v) \cdot v' = Df(v0 + C * r0) \cdot (v0' + C' * r_0) = ...
	//
	// i obliczamy bardzo delikatnie, wyrazenie z macierzy:
	//                            ( 0         0         0  )
	// Df(v0 + C * r0) = Df(v0) + ( 0         0         0  )  =    Df(v0) + A
	//                            ( Cz *r_0   0  Cx * r_0 )
	// Cz oznacza trzeci wiersz macierzy C (odpowiada z)
	//
	// wydaje się, że tutaj nie za bardzo da się wyjąć r_0 na zewnątrz, bez angażowania tensorów, ale...
	// popatrzymy na A*v0': (na A * C' * r_0 nie patrzymy to to wyraz drugiego rzędu względem r0!)
    //
	// A*v0' = (0, 0, x0'*(Cz*r_0) + z0'*(Cx* r0))^T i teraz możemy to zamienić
	//       = (0      0       0     )
	//         (0      0       0     )  * r0    = A * r0
	//         ( x0' * Cz + z0' * Cx )
	// i ostatnie wiersz to rzeczywiście wiersz! Nadużywam notacji i nazywam A te nową macierz.
	//
	// Ostatecznie, C'' = (Df(v0)C' + A) * r0
	//
	//
	// więc v^[2] = 1/2 (powyżej)
	// i tak można to pociągnąc rekurencyjnie wyżej, jakbysmy chcieli (nie jest tak prosto bo wchodzi D^2f...
	// v^(3) = Df(v) \cdot v'(t) =
	//       = ([-y'(t) - z'(t), x'(t) + b*y'(t), z(t)*x'(t) + x(t)*z'(t)]^T)' =
	//       = [-y'' - z'' + 0, x'' - by'' + 0, (zx'' + xz'') + (2*z'x')]^T  (***)
	// LUB
	//       = (D^2f(v) \circ v') \circ v' + Df(v) \cdot v''
	//       (pierwszy człon znika, gdy f jest liniowe (Df(v) niezależne od v)
	//       (drugi to po prostu liniowa część)
	//       we wzorze (***) dla x i y wystepuje tylko liniowa część,
	//       dla z mamy liniowa (pierwszy nawias ()) i nieliniową (drugi nawias)


	// UWAGA: dla odwzorowania w tył zmień znaki odpowiednie!

	// TODO: zastanowić się jak to się ma do Automatycnego Różniczkowania??

}
