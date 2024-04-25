//////////////////////////////////////////////////////////////////////////////
///
///  @file utils.h
///  
///  @author ag  @date   Mar 5, 2020
//////////////////////////////////////////////////////////////////////////////


#ifndef _EXAMPLES_PROJECTSTARTER_UTILS_H_
#define _EXAMPLES_PROJECTSTARTER_UTILS_H_

#include "capd/intervals/lib.h"
#include "capd/capdlib.h"
#include <iostream>
#include "capd/covrel/HSet2D.h"
#include "capd/dynsys/DynSysMap.h"
#include "capd/covrel/HSetMD.h"

#include <capd/ddes/ddeslib.h>
#include <capd/ddeshelper/ddeshelperlib.h>
#include <capd/ddeshelper/DDEHelperRigorous.hpp>
#include "equation.h"

using namespace std;
using namespace capd;			

typedef capd::covrel::HSet2D<DMatrix, IMatrix> HSet2D;
typedef capd::dynsys::DynSysMap<IMap> DynSysMap; ///TODO: czy potrzebne ???
typedef capd::covrel::HSetMD<DMatrix,IMatrix> HSet;
typedef capd::covrel::GridSet<IMatrix> GridSet;

/**
 * we will use Rossler DDE with interval parameters.
 * This is the definition of f from x'(t) = f(x(t), x(t-tau)).
 */
typedef Rossler<interval> Eq;

/** this contains all necessary ingredients to do rigorous numerics in DDEs */
typedef capd::ddeshelper::RigorousHelper<Eq, 1, IMatrix, IVector> Setup;

// below are just renaming for shorter class names
typedef Setup::Grid Grid;					// grid for the solutions in the C^n_p space.
typedef Setup::DDEq DDEq;					// this is F from the abstract formulation x'(t) = F(x_t), that contains the information on delays. In our case F(x_t) := f(x_t(0), x_t(-tau)).
typedef Setup::Solution DDESolution;		// basic set type for DDEs
typedef Setup::SetType AmbientSet;		    // ??? TODO: docs
typedef Setup::Solver DDESolver;			// semidynamical system
typedef Setup::PoincareMap DDEPoincareMap;	// Poincare map for the semidynamical system
typedef Setup::Section DDESection;			// basic section for DDEs


//////////////////////////////////////////////////////////////////////////////
///	Some auxiliary functions
//////////////////////////////////////////////////////////////////////////////

// cuts the first coordinate of an interval vector
IVector cut(IVector x);	
DVector cut(DVector x);	
IMatrix cut(IMatrix M);	
			
//prepends an interval 2-vector (or 2x2 matrix) to 3-vector with zero
IVector expand(IVector x);	
DVector expand(DVector x);
IMatrix expand(IMatrix M);	

/////////////////////////////////////////////////////////////////////////
/// SecMap class
/////////////////////////////////////////////////////////////////////////

// 2D Poincare Map image of an interval vector on the x=0 section 
class SecMap	
{
public:
	DDEPoincareMap *P;
	Grid *grid;
	SecMap(DDEPoincareMap &Q, Grid &g) {P=&Q; grid = &g;}	
	SecMap(){};						
	
//	IVector image(const IVector &x,IMatrix &DP, const int iter = 1) const;
	IVector image(const IVector &x, const int iter = 1) const;
//	IVector operator()(const IVector &x,IMatrix &DP, const int iter = 1) const {return image(x,DP,iter);}
	IVector operator()(const IVector &x, const int iter = 1) const {return image(x,iter);}
	
};

//////////////////////////////////////////////////////////////////////////////
///	system3d structure
//////////////////////////////////////////////////////////////////////////////

// stores the Roessler system
struct system3d
{
int p;
interval tau, h;
interval a;
Grid grid;
int order;
Eq rhs;
DDEq vf;
DDESection section;
DDESolver solver;
DDEPoincareMap P;
SecMap P2d;
int iteration;

system3d(const interval &aa = interval(5.7)) :									// constructor, default parameter a=5.7
	p(64), tau(1.0), h(tau/p), grid(h), a(aa),
	order(3),
	rhs(aa, interval(0.2),0.00001),														// f(a, b,eps)
	vf(rhs, grid(p)), 															// f(x_t), grid(p) = p*h = tau
	section(3, 0),																// x=0 section
	solver(vf, order * 3),														//max order = order*3	
	P(solver, section, poincare::MinusPlus),
	P2d(P, grid),
	iteration(1)
	{
		P.setRequiredSteps(order * p);	
		P.setMaxSteps(1000);	
	}

bool inside(const HSet2D &hset1, const HSet2D &hset2, int howManyPiecesH=1, int howManyPiecesV=1, int iteration = 1);
void makeHistory(const HSet2D &hset);
};


#endif //_EXAMPLES_PROJECTSTARTER_UTILS_H_


