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
bool system3d::inside(const HSet2D &hset1, const HSet2D &hset2, int howManyPiecesH, int howManyPiecesV, int iteration) 
{
	GridSet gridSet(2);																		// stores grid
	hset1.gridSet(gridSet,howManyPiecesH,howManyPiecesV);		
	HSet2D hset2Base(IVector({0.,0.}),IMatrix::Identity(2),hset2.get_r());
  
	//C0Rect2Set Set3d(expand(hset1.center()));

	IVector Pset2d(2);  
	interval rt(0.);
	IMatrix expandedGridSetCoordinateSystem = expand(gridSet.coordinateSystem());			// expanded to 3d, where P is defined
	IVector expandedGridSetBox = expand(gridSet.box());
	IMatrix expandedHSet2InvCoordinateSystem = expand(hset2.invCoordinateSystem());
	IVector expandedHSet2Center = expand(hset2.center());
	
	bool liesInside = 1;
	P.setRequiredSteps(0);
		auto zero = grid(0);
		auto tau = grid(p);
	for(auto i = gridSet.begin();i!=gridSet.end();++i)
	{
		// Set3d = C0Rect2Set(expand(*i),expandedGridSetCoordinateSystem,expandedGridSetBox);
		//Set3d = DDESolution(grid, -tau, zero, 3, expand(*i));				// (grid, from, to, order, constant function value)
		DDESolution segment(grid, -tau, zero, 3, AmbientSet(expand(*i),expandedGridSetCoordinateSystem,expandedGridSetBox));				// (grid, from, to, order, constant function value)
		//Set3d = DDESolution(grid, grid(0), expand(*i));
		//Set3d.getValueAtCurrent().set_Cr0(expandedGridSetCoordinateSystem, expandedGridSetBox);


		cout << "IVector " << IVector(segment.getValueAtCurrent()) << endl;
		auto Psegment = P(segment, rt);
		// PSet3d.affineTransform(expandedHSet2InvCoordinateSystem, expandedHSet2Center);


		AmbientSet PSet3d = Psegment.getValueAtCurrent();
		IVector PV = 
			expandedHSet2InvCoordinateSystem * (PSet3d.get_x() - expandedHSet2Center) + 
			(expandedHSet2InvCoordinateSystem * PSet3d.get_C()) * PSet3d.get_r0() +  
			(expandedHSet2InvCoordinateSystem * PSet3d.get_B()) * PSet3d.get_r();

		// rt = P.getLastReachTime();
		cout << "rt = " << rt << endl;

		Pset2d = cut(PV);		// the image of the small box

		liesInside=liesInside && hset2Base.inside(Pset2d);									// do all lie inside?
		
		if(!liesInside) 
			{
				cout << "Does not lie inside, for this tiny box: " << Pset2d << " does not lie inside the box " << hset2Base << endl;
				return 0;
			}
	}
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
