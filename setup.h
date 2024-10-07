#ifndef _SHA_ROS_DDE_SETUP_
#define _SHA_ROS_DDE_SETUP_

/**
 * We are here building on top of the results from work
 * [GZ2022] A~Gierzkiewicz and P.~Zgliczy\'nski. ,,From the {S}harkovskii theorem to periodic orbits for the {R}\"ossler system.'', Journal of Differential Equations, 314:733--751, (2022).
 *
 * In this file we gather common definitions, that will be used across all programs.
 * In particular, we define here the sets from the work [GZ2022].
 *
 * You can change this setup and recompile all files, if you want to prove something by yourself.
 */

#include "utils.h"

struct Config {
	interval a = 5.25; 									// this is representable number, so we can left it that way
	interval epsi = interval(1.) / interval(10000.); 	// [0.0001], as this is not representable
	interval tau = 0.5;									// this is representable number, so we can left it that way
	system3d roessler525;
	HSet2D grid3;
	vector<HSet2D> c3;


	Config(interval par_epsi): roessler525(a, par_epsi, tau) {
		// the data for the sets and their names are taken from [GZ2022]
		grid3 = HSet2D(
					IVector ({-6.38401, 0.0327544}),
					IMatrix ({{-1., 0.000656767}, {-0.000656767, -1.}}),
					DVector({3.63687,0.0004})
				);
		c3.resize(3);
		c3[0] = HSet2D(IVector ({-3.46642, 0.0346316}) , IMatrix ({{-1., 0.000656767}, {-0.000656767, -1.}}) , DVector({0.072,0.00048}));			// cube 1
		c3[1] = HSet2D(IVector ({-6.26401, 0.0326544}) , IMatrix ({{-1., 0.000656767}, {-0.000656767, -1.}}) , DVector({0.162,0.00066}));			// cube 2
		c3[2] = HSet2D(IVector ({-9.74889, 0.0307529}) , IMatrix ({{-1., 0.000656767}, {-0.000656767, -1.}}) , DVector({0.036,0.00072}));			// cube 3
	}

	Config() { Config(epsi); }
} config;


#endif // _SHA_ROS_DDE_SETUP_


