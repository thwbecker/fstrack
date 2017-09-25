//
// real constants (some of those are also still set in input_para.c
// $Id: constants.h,v 1.11 2010/12/31 18:59:17 becker Exp $
//
#ifndef SINGLE_PRECISION
// Runge Kutta eps default precision factor
#define RK_EPS_DEF 1.0e-5
#else
#define RK_EPS_DEF 1.0e-5
#endif
// max time bailout default for backward strain [Ma]
#define MTIME_BAILOUT_DEF 43.0
// max strain default bailout value (log strain as in mstrain.c)
#define MSTRAIN_BAILOUT_DEF 0.25
//
// given velocities in cm/yr, a nice timescale might be 1e6 yr
// give this constant in units of [yrs]
//
#define TIMESCALE 1e6
//
// timescale in [s]
#define TIMESCALE_S (TIMESCALE * 365.25 * 24. * 3600.)
//
// max stepper increment for time in characterstic times
//
#define HMAX_DEF 1.0
// Earth's radius in [km]
#define R_E 6371.0087714
//
// radius of CMB in non-dim r units (0...1: surface), [non. dim]
#define CMB_R 0.5448
// limit for lowermost layer detection
#define CMB_RL 0.5463
// 410 in [non. dim]
#define P410_RL 0.93564589546382
//
// additional bailout criteria for max strain bailout mode
//
// bailout depth for backtracing strain if tracers are deeper than this level
// [non. dim]
// applies for -bs and -btd 
//
#define STRAINMAX_DEPTH_BAILOUT P410_RL
// initial tracer depth in km if not read from tdepth.dat file [km]
#define INIT_TRACER_DEPTH_DEF 200.0

// limit for max initial timestep
#define MAX_DT_INIT 1.0
/*

ISA related

parameter values determined empirically
*/
#define ISA_PI_CRIT_DEF 1e-2
/*

default number of tracers

*/
#define DEF_NR_TRACERS 500
/*

Lyapunov related

*/
#define LYA_RENORM_DEF 10.0// renorm if L(i) > renorm
#define LYA_RLBAILOUT_DEF 1.0e-6// bailout at this maximum fractional change
#define LYA_TMAX_DEF 10000.0// this needn't be a realistic time necessarily
#define LYA_TMIN_DEF 50.0 	/* minimum time before bailout is allowed */

#ifdef USE_DREX
/* 

DREX related

*/

/* 
   angles to rotate enstatite to align with olivine in flow
*/
#include "drex.h"

#define ENSTATITE_EULER {90.0,90.0,0.0}


#define DREX_YERRFAC_DEF -0.005;	/* 
					   by default, the maximum timestep for KR computations
					   is 0.01 the inverse of the max strain rate. this was found to 
					   yield good misfits (almost indistinguishable from finer stepping) at
					   about eight times the execution times than no control
					*/
#define DREX_FAVG 0;		/*  0: voigt 0.5: VHR 1: reuss, default is Voigt */

#endif
