#ifndef _SHA_ROS_DDE_UTILS_
#define _SHA_ROS_DDE_UTILS_

/**
 * We are here building on top of the results from work
 * [GZ2022] A~Gierzkiewicz and P.~Zgliczy\'nski. ,,From the {S}harkovskii theorem to periodic orbits for the {R}\"ossler system.'', Journal of Differential Equations, 314:733--751, (2022).
 * So we inherited the structure from the computer assisted part there.
 *
 * You should check the 'utils.h' file there.
 */

// this tells the compilation system to have
#define DDES_ALLOW_SYSTEM

#include <iostream>
#include <vector>
#include <numeric>

// those are taken from the work [GZ2022]
// we will be building on top of data from this work
#include "capd/capdlib.h"
#include "capd/intervals/lib.h"
#include "capd/covrel/HSet2D.h"
#include "capd/dynsys/DynSysMap.h"
#include "capd/covrel/HSetMD.h"

// those add access to DDE codes.
#include <capd/ddes/ddeslib.h>
#include <capd/ddeshelper/ddeshelperlib.h>
#include <capd/ddeshelper/DDEHelperRigorous.hpp>
#include <capd/ddeshelper/DDEHelperNonrigorous.hpp>
#include <capd/ddeshelper/DDEHelperDrawing.h>
#include "equation.h"

using namespace std;
using namespace capd;			

/** this is standard def. from CAPD used by the [GZ2022] code */
typedef capd::covrel::HSet2D<DMatrix, IMatrix> HSet2D;
typedef capd::dynsys::DynSysMap<IMap> DynSysMap;
typedef capd::covrel::HSetMD<DMatrix,IMatrix> HSet;
typedef capd::covrel::GridSet<IMatrix> GridSet;

/**
 * we will use Rossler DDE with interval parameters.
 * This is the definition of f from x'(t) = f(x(t), x(t-tau)).
 */
typedef Rossler<interval> Eq;
/**
 * This is the equation for non-rigorous computations
 * done when generating good coordinate frame.
 */
typedef Rossler<double> DEq;

/** this contains all necessary ingredients to do rigorous numerics in DDEs */
typedef capd::ddeshelper::NonrigorousHelper<DEq, 1, DMatrix, DVector> Numerics;

/** this contains all necessary ingredients to do rigorous numerics in DDEs */
typedef capd::ddeshelper::RigorousHelper<Eq, 1, IMatrix, IVector> Setup;

// below just renaming for use in the code
// We use prefix 'D' to mark non-rigorous version
typedef Numerics::Grid DGrid;				// grid for basic curve
typedef Numerics::DDEq DDDEq;
typedef Numerics::Solution DSolution;		// basic description for a curve
typedef Numerics::Solver DSolver;			// semidynamical system
typedef Numerics::JetSection DSection;
typedef Numerics::PoincareMap DPoincare;

// below are just renaming for shorter class names
typedef Setup::Grid Grid;					// grid for the solutions in the C^n_p space.
typedef Setup::DDEq DDEq;					// this is F from the abstract formulation x'(t) = F(x_t), that contains the information on delays. In our case F(x_t) := f(x_t(0), x_t(-tau)).
typedef Setup::Solution DDESolution;		// basic set type for DDEs
typedef Setup::SetType AmbientSet;		    // I call 'Ambient' the space of \R^3, in which the head of the solution lives: z(x_t)
typedef Setup::Solver DDESolver;			// semidynamical system
typedef Setup::PoincareMap DDEPoincareMap;	// Poincare map for the semidynamical system
typedef Setup::Section DDESection;			// basic section for DDEs

/** just a forward declaration, so that helper functions see it */
struct system3d;


//////////////////////////////////////////////////////////////////////////////
///	Some auxiliary functions
//////////////////////////////////////////////////////////////////////////////

/**
 * nonrigorous. Makes a nice gnuplot of many segments of solutions (stored as a vector of DVector)
 * and saves it under appropriate file in the current working directory (wd).
 *
 * It interprets each DVector as a segment of solutions in the system sys.
 * NOTE: wd must have trailing '/'! Default is "./"
 */
void splotMany(system3d& sys, std::vector<DVector>& segs, std::string const& name, std::string const& wd="./");

// the projection of a vector (x, y, z, ...) to 2 dimensions (y, z). Unless the vector is already 2 dim, when it just returns it.
template<typename VT> VT cut(VT const& x){
	if(x.dimension()<2) throw std::logic_error("cut(vector): too few dimensions!");
	if(x.dimension()==2) return x;
	VT y(2);
	for (int i = 0; i < 2; ++i)	y[i]=x[i+1];
	return y;
}
// cuts a matrix to (2x2), removing leading row and column. Unless the matrix is already 2x2, when it just returns it.
IMatrix cut(IMatrix const& M);
			
//prepends a vector with an extra dimension at the front, zero on that extra dim. If is already 3 dim, then return it. S
template<typename VT> VT expand(VT const& x){
	if(x.dimension()==3) return x;
	if(x.dimension()!=2) throw std::logic_error("expand(vector): bad dimensions");
	VT y({0.,0.,0.});
	for (int i = 0; i < 2; ++i)	y[i+1]=x[i];
	return y;
}
//prepends a matrix with an extra dimension at the front, with identity on this block. If is already 3 dim, then return it. S
IMatrix expand(IMatrix const& M);

//////////////////////////////////////////////////////////////////////////////
///	system3d structure
//////////////////////////////////////////////////////////////////////////////

/*
 * stores the Roessler system with the DDE term.
 * it also stores all data necessary to integrate segments forward in time
 * and gives helper methods for various activities, like crating segments,
 * computing Poincare maps to section \pi_x z(X) = 0, etc.
 */
struct system3d
{
	static const int d = 3;
	int p;
	interval tau, h;
	interval a, eps;
	Grid grid;
	DGrid dgrid;
	int order;
	Eq rhs;
	DDEq vf;
	DDESection section;
	DDESolver solver;
	DDEPoincareMap P;
	int iteration;

	system3d(
			interval tau = interval(1.0),
			interval eps = interval(0.00001),
			interval aa = interval(5.7)
	):	// constructor, default parameter a=5.7
		p(64), tau(tau), h(tau/p), grid(h), dgrid(h.leftBound()), a(aa),
		order(3),
		rhs(aa, interval(0.2), eps),												// f(a, b,eps)
		vf(rhs, grid(p)), 															// f(x_t), grid(p) = p*h = tau
		section(3, 0),																// x=0 section
		solver(vf, order * 3),														//max order = order*3
		P(solver, section, poincare::MinusPlus),
		iteration(1)
		{
			cerr << "p=" << p << endl;
			cerr << "n=" << order << endl;
			cerr << "tau=" << tau << endl;
			cerr << "h=" << h << endl;
			cerr << "epsi=" << eps << endl;
			cerr << "a=" << a << endl;
			P.setRequiredSteps(order * p);
			P.setMaxSteps(1000);
		}

	/**
	 * checks if the image of one sets is inside the second set.
	 * \pi_{y,z} z(X_i) \subset hset1 \subset \R^2,
	 * CloseTail(X_i) \subset (Id - z(.))  (mid_i + C*r0)
	 * FarTail(X_i) \subset Xi_0
	 * (so that Close and Far tails are the same, only the heads given by hset1 changes)
	 */
	bool inside_piece(
		const HSet2D &hset1, IVector const& mid1,
		const HSet2D &hset2, IVector const& mid2,
		const IVector& r0, const IVector& Xi0,
		const IMatrix& C, const IMatrix& invC,
		int howManyPiecesH, int howManyPiecesV, int pieceH, int pieceV,
		IVector &outPimage, IVector &outPXi
	);

	// returns the M - size of the finite representation of the segments.
	int M() const { return (1 + (1+order) * p) * d; }

	// makes a rig segment, not a very DRY code, but the capdDDEs library forces that...
	DSolution makeDSegment(DVector vector){ DSolution sol(dgrid, -dgrid(p), -dgrid(0), order, {0.,0.,0.}); sol.set_x(vector); return sol; }
	// makes a rig segment, not a very DRY code, but the capdDDEs library forces that...
	DDESolution makeISegment(IVector vector){ DDESolution sol(grid, -grid(p), -grid(0), order, {0.,0.,0.}); sol.set_x(vector); return sol; }
};

// This will be used to turn on/off the debug in programs.
#define SHA_DEBUG(s) if(verbose) std::cerr << s << std::flush

#endif // _SHA_ROS_DDE_UTILS_


