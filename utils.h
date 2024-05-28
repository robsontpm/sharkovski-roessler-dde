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
#include <capd/ddeshelper/DDEHelperNonrigorous.hpp>
#include <capd/ddeshelper/DDEHelperDrawing.h>
#include "equation.h"

#include <vector>
#include <array>
#include <numeric>

using namespace std;
using namespace capd;			

typedef capd::covrel::HSet2D<DMatrix, IMatrix> HSet2D;
typedef capd::dynsys::DynSysMap<IMap> DynSysMap; ///TODO: czy potrzebne ???
typedef capd::covrel::HSetMD<DMatrix,IMatrix> HSet;
typedef capd::covrel::GridSet<IMatrix> GridSet;

//namespace capd{
//namespace covrel{
//
///**
// This class provides a h-set in an arbitrary dimension
// */
//template<typename MatrixT, typename IMatrixT>
//class MyHSetMD: capd::covrel::HSetMD<MatrixT, IMatrixT> {
//public:
//  typedef MatrixT MatrixType;
//  typedef IMatrixT IMatrixType;
//  typedef typename MatrixType::RowVectorType VectorType;
//  typedef typename IMatrixType::RowVectorType IVectorType;
//  typedef typename IVectorType::ScalarType ScalarType;
//  typedef typename VectorType::ScalarType BoundType;
//  typedef typename VectorType::size_type size_type;
//
//  // we assume that the set is represented as center + Base * B(0,r)
//  // uDim denotes a number of unstable directions
//  // sDim denotes a number of stable directions
//  MyHSetMD() {}
//  MyHSetMD(const VectorType& center, const MatrixType& Base, const IMatrixType& InvBase, size_type uDim, size_type sDim, const VectorType& r);
//  MyHSetMD(const IVectorType & center, const IMatrixType& Base, const IMatrixType& InvBase, size_type uDim, size_type sDim, const VectorType& r);
//  virtual ~MyHSetMD() {}
//
//  /**
//   * in the \c grid it returns uniform grid of the given face
//   */
//  template<typename IMatrix>
//  GridSet<IMatrix>& gridFace(
//      GridSet<IMatrix>& grid,
//      const std::vector<size_type>& gridSizes,
//      const std::vector<size_type>& dimensions,
//      size_type totalDimension,
//      size_type coordinateToFix,
//      Side side = bothSides
//  ) const;
//
//  /// this procedure creates a grid of the whole h-set
//  /// in the following form: G.center[i] + G.C * G.r
//  /// d is a vector of indices of coordinates if the set is embeded in higher dimension
//  template<typename IMatrix>
//  GridSet<IMatrix>& gridSet(
//      GridSet<IMatrix>& G,
//      const std::vector<size_type>& grid,
//      const std::vector<size_type>& d,
//      size_type totalDimension
//  ) const;
//
//
//  const IVectorType & center() const {return m_Ix;} ///< center of h-set
//  const IMatrixType & coordinateSystem() const {return m_IB;} ///< matrix of base vectors
//  const IMatrixType & invCoordinateSystem() const {return m_invIB;}
//  IVectorType box() const { ///< h-set in 'proper' coordinate system (it is a product of intervals)
//    IVectorType b(m_r.dimension());
//    for(size_type i=0; i<m_r.dimension(); ++i)
//    b[i] = ScalarType(-m_r[i],m_r[i]);
//    return b;
//  }
//  const VectorType & radius() {return m_r;} ///< returns vector of radiuses (in each direction)
//  size_type unstableDimension() const {return m_uDim;} ///< number of unstable dimensions
//  size_type stableDimension() const {return m_sDim;} ///< number of stable dimensions
//  const VectorType & get_x() const { return m_x;}
//  const MatrixType & get_B() const { return m_B;}
//  const VectorType & get_r() const { return m_r;}
//  const IVectorType & get_I_x() const { return m_Ix;}
//  const IMatrixType & get_I_B() const { return m_IB;}
//  const IMatrixType & get_I_invB() const { return m_invIB;}
//
//  virtual std::string show() const;
//  std::string showInfo() const;
//  friend std::ostream & operator << (std::ostream & stream, const HSetMD & set ){
//  		stream << set.show();
//  		return stream;
//  }
//protected:
//  VectorType m_x; ///< center point
//  MatrixType m_B; ///< coordinate system
//  IVectorType m_Ix; ///< rigorous center point
//  IMatrixType m_IB, m_invIB; ///< rigorous coordinate system and its inverse
//
//  size_type m_uDim, ///< number of unstable dimensions
//            m_sDim; ///< number of stable dimensions
//  VectorType m_r; ///< radiuses of balls (used both for non-rigorous and rigorous version)
//}; // class HSetMD
//
//}// namespace covrel
//}// namespace capd

/**
 * we will use Rossler DDE with interval parameters.
 * This is the definition of f from x'(t) = f(x(t), x(t-tau)).
 */
typedef Rossler<interval> Eq;

/** this contains all necessary ingredients to do rigorous numerics in DDEs */
typedef capd::ddeshelper::NonrigorousHelper<Eq, 1, DMatrix, DVector> Numerics;

/** this contains all necessary ingredients to do rigorous numerics in DDEs */
typedef capd::ddeshelper::RigorousHelper<Eq, 1, IMatrix, IVector> Setup;

// below just renaming
typedef Numerics::Grid DGrid;				// grid for basic curve
typedef Numerics::Solution DSolution;		// basic description for a curve
typedef Numerics::JetSection DSection;
typedef Numerics::PoincareMap DPoincareMap;

// below are just renaming for shorter class names
typedef Setup::Grid Grid;					// grid for the solutions in the C^n_p space.
typedef Setup::DDEq DDEq;					// this is F from the abstract formulation x'(t) = F(x_t), that contains the information on delays. In our case F(x_t) := f(x_t(0), x_t(-tau)).
typedef Setup::Solution DDESolution;		// basic set type for DDEs
typedef Setup::SetType AmbientSet;		    // ??? TODO: docs
typedef Setup::Solver DDESolver;			// semidynamical system
typedef Setup::PoincareMap DDEPoincareMap;	// Poincare map for the semidynamical system
typedef Setup::Section DDESection;			// basic section for DDEs

typedef std::pair<DVector, DMatrix> AffineCoordinates;
typedef std::vector<AffineCoordinates> SimpleHistory;
typedef std::vector<AffineCoordinates> SimpleJet;
typedef std::vector<SimpleJet> HighOrderHistory;


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
SecMap P2d;
int iteration;

typedef HighOrderHistory HistoryType;

system3d(const interval &eps = interval(0.00001), const interval &aa = interval(5.7)):	// constructor, default parameter a=5.7
	p(64), tau(0.5), h(tau/p), grid(h), dgrid(h.leftBound()), a(aa),
	order(3),
	rhs(aa, interval(0.2), eps),												// f(a, b,eps)
	vf(rhs, grid(p)), 															// f(x_t), grid(p) = p*h = tau
	section(3, 0),																// x=0 section
	solver(vf, order * 3),														//max order = order*3	
	P(solver, section, poincare::MinusPlus),
	P2d(P, grid),
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

bool refine_box(const HSet2D &hset1, const HSet2D &hset2, int howManyPiecesH=1, int howManyPiecesV=1, int iteration = 1);
bool inside(const HSet2D &hset1, const HSet2D &hset2, int howManyPiecesH=1, int howManyPiecesV=1, int iteration = 1);
bool inside_piece(const HSet2D &hset1, const HSet2D &hset2, int howManyPiecesH, int howManyPiecesV, int pieceH, int pieceV);
void makeHistory(const HSet2D &hset);
HistoryType makeHistoryD(const DVector& v0, const DMatrix& C, DVector& out_x0, std::vector<DVector>& out_coords);
void makeHistoryC0(const DVector& v0, DVector& out_x0);

int M() const { return (1 + (1+order) * p) * d; }

DSolution makeDSegment(DVector vector){ DSolution sol(dgrid, -dgrid(p), -dgrid(0), order, {0.,0.,0.}); sol.set_x(vector); return sol; }
DDESolution makeISegment(IVector vector){ DDESolution sol(grid, -grid(p), -grid(0), order, {0.,0.,0.}); sol.set_x(vector); return sol; } // TODO: check if this works...
};


#endif //_EXAMPLES_PROJECTSTARTER_UTILS_H_


