#ifndef EXAMPLES_ROSSLER_ODE_VS_DDE_CODE_EQUATION_H_
#define EXAMPLES_ROSSLER_ODE_VS_DDE_CODE_EQUATION_H_

/**
 * This class defines f in the right hand side of the DDE
 * (1) x'(t) = f(x(t), x(t-tau), ...)
 *
 * Here, we define DDE perturbation of R\"ossler ODE of the following form:
 *
 * (2) x'(t) = f(x(t)) + \epsilon * f(x(t - \tau))
 *
 * Please note that the \tau is not the part of the parameters, as this is the parameter
 * of the phasespace C^0([-\tau, 0], \R^3).
 *
 * It is a little harder than in the CAPD itself, as CAPD allows
 * for specifying RHS as a simple string. Here, the approach is similar
 * to specifying  RHS as a function Node f(Node t, Node in[], int inDim, Node out[], int outDim, Node params[], int parDim)
 * used for example in ./external/capd/capdDynSys4/examples/poincare/IPoincareMapExample.cpp
 *
 * This does not need to be template class, but I like to do it that way so that i can use it for
 * rigorous and non-rigorous computations!
 */
template<typename ParamSpec>
class Rossler{
public:
	/** this is not required */
	typedef unsigned int size_type;
	/** RigorousHelper requires to define explicitly what is the type of parameters! It must be named ParamType! */
	typedef ParamSpec ParamType;
	/** RigorousHelper requires to define explicitly what is the type of a vector of parameters! It must be named ParamsVectorType! */
	typedef capd::vectalg::Vector<ParamType, 0> ParamsVectorType;

	/** default constructor for default values of parameters, its good to have one. */
	Rossler(): a(ParamType(57) / ParamType(10)), b(ParamType(2) / ParamType(10)), eps(0) {}
	/** construct eq. for given parameter values */
	Rossler(ParamType a, ParamType b, ParamType eps): a(a), b(b), eps(eps) {}
	/** construct eq. taking parameters from a vector. It must be of this form for RigorousHelper to work! */
	Rossler(ParamsVectorType const& params): a(params[0]), b(params[1]), eps(params[2]) {}

	/** output dimension of the map f in DDE (1), R\"ossler is 3D  */
	static size_type imageDimension() { return 3; }
	/**
	 * input dimension of the map f in DDE (1), should be (1+m) * imageDimension() for m delays,
	 * since we do ODE (2), we set m = 0
	 * */
	static size_type dimension() { return 6; }
	/** required by the RigorousHelper */
	static size_type getParamsCount() { return 3; }

	/**
	 * computation of RHS of the equation x'(t) = f(x(t), x(t-tau), ...),
	 *
	 * It must be template, as we use autodiff to generate derivatives of this.
	 *
	 * Let d = imageDimension(). We follow the convention that
	 * InVectorSpec x = (x_1(t), x_2(t), ..., x_d(t), x_1(t-tau), x_2(t-tau), ..., x_d(t-tau), ...)
	 * Therefore in code below, e.g. x[1] is y coordinate in the original Rossler ODE.
	 * (Please remember, that indices start from 0 in InVectorSpec x)
	 *
	 * The t is current time for non-autonomous equations.
	 * The fx parameter (output), fx means 'f(x)', is the output of the operator.
	 * It is for convenience of calling operator
	 * (as the return type cannot be auto-deduced in C++!),
	 * and also for avoiding excessive copy-constructing, move-copying, etc.
	 * It just alter the value of already existing vector given by reference.
	 */
	template<typename RealSpec, typename InVectorSpec, typename OutVectorSpec>
	void operator()(const RealSpec& t, const InVectorSpec x, OutVectorSpec& fx) const {
		fx[0] = -(x[1]+x[2])    - eps*(x[4]+x[5]);
		fx[1] =  x[0]+b*x[1]    + eps*(x[3]+b*x[4]);
		fx[2] = b+x[2]*(x[0]-a) + eps*(b+x[5]*(x[3]-a));
	}

protected:
	ParamType a;
	ParamType b;
	ParamType eps;
};

#endif
