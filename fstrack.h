//
// header files with constants for fstrack
// $Id: fstrack.h,v 1.37 2016/09/05 04:44:58 becker Exp $
//
//#define FSTRACK_DEBUG
#ifndef __FSTRACK_HEADER_SOURCED__

//
//
// for calculation the strain bailout using the absolute size of the eigenvalue
//#define MAX_EVAL_STRAIN
// for using the fraction between e1/e2 and e2/e3, whichever is larger
#define MAX_FRAC_EVAL_STRAIN
// 
//

#include <stdio.h>
#include <stdlib.h>
//#include <malloc.h> // is this maybe not needed anymore?
#include <string.h>
#include <limits.h>
#include <time.h>
#include <math.h>

/*  
    filenames are defined here
*/
#include "filenames.h"

/* 
   
precision settings

*/
#include "precision.h"

//
// length of string variables
#define STRLEN 300
//
// FORTRAN routines underscore calling convention
//
#ifndef NO_UNDERSCORE
/* regular ones */
#define ellsphere_cart ellsphere_cart_
#define cart_vel_int cart_vel_int_
#define invert3x3 invert3x3_
#define dgpadm dgpadm_
#define odeintell odeintell_
#define gramschmidt gramschmidt_
#define inter_vel inter_vel_
#define svdcmp svdcmp_

#endif
#ifdef USE_DREX
/* DREX stuff */
#include "drex.h"
#endif



#define MY_USHRT_MAX 65535

/* 

   PREM STUFF

 */
#include "prem.h"
/* 

splitting analysis stuff

*/
#include "vera_util.h"
#define MAX_NR_HARM 10		/* max number of harmonics for
				   analysis */

/* 

for the model structure, which holds a few states

*/
#define NR_LOC_STATES 4
#define INITIAL 0
#define ORIG 1 
#define FINAL 2
#define TEMP 3

/* 

surface wave sensitivities, list of specified periods in [s]

 */

//#define N_SW_SENS 3		/* three periods */
//#define SW_SENS_PERIODS {50,100,150} /* of those second */
#define N_SW_SENS 12		
#define SW_SENS_PERIODS {25,29,33,40,50,67,83,100,111,125,143,150}
//
// definition of model, tracer, and state structures
//
#include "structures.h"
//
// index and type conventions
//
#
#define FSTRACK_R 0       /* for the spherical location vectors */
#define FSTRACK_THETA 1
#define FSTRACK_PHI 2
#define TPPROD 3  /* those two are only used for the n[] vector */
#define NRNTNP 4
#define FSTRACK_X 0       /* for the cartesian location vectors */
#define FSTRACK_Y 1
#define FSTRACK_Z 2
#define FSTRACK_E1 2      /* the eigenvalues returned by 
		     calc_eigensystem_sym are sorted ascendingly,
		     hence if e1 > e2 > e3, those are the indices */
#define FSTRACK_E2 1
#define FSTRACK_E3 0 


// macro for C style refence to velocity field at depth level i
// theta row j and phi column k 
#define VOFF(i,j,k) ((i)*MDP->n[TPPROD]+(j)*MDP->n[FSTRACK_PHI]+(k))
// 3x3 matrix elements
// RR RT RP   0 1 2
// TR TT TP   3 4 5
// PR PT PP   6 7 8
#define RR 0 // first six are used to initialize symmetric matrices
#define RT 1
#define RP 2
#define TT 4
#define TP 5
#define PP 8
#define TR 3 // these are the (possibly) symmetric parts
#define PR 6
#define PT 7
/* those are the same, really */
#define XX 0 // first six are used to initialize symmetric matrices
#define XY 1
#define XZ 2
#define YY 4
#define YZ 5
#define ZZ 8
#define YX 3 // these are the (possibly) symmetric parts
#define ZX 6
#define ZY 7

//
// operational modes for the program
//
// init tracer modes
#define DIST_EVEN_AREA 0   
#define DIST_EVEN_DX 1
#define SPOTTED_LAT_FROM_FILE 2
#define SPOTTED_3D_FROM_FILE 3
#define SPOTTED_3D_WITH_ATTR 4
// all tracer output modes, see output.c
#define ALL_TRACER_RTPT 0	/* location in spherical coordinates 

				*/
#define ALL_TRACER_LLZT 1	/* lon lat depth time */
#define ALL_TRACER_LLZAT 2	// lon lat depth attribute(s) time
#define ALL_TRACER_STRAIN_COMP 3 /* location and left stretch strain tensor L */
#define ALL_TRACER_STRAIN_EIGEN 4 /* location and eigen axis of LS strain (FSE) */
#define ALL_TRACER_STRAIN_EIGEN_ANGLE 5	/* angles of those axis */
#define ALL_TRACER_DEFORMATION 6 /* location and deformation F  */
#define ALL_TRACER_LYA 7	/* location and lyapunov  */
#define ALL_TRACER_ISA 8	/* ISA axis */
#define ALL_TRACER_TI 9	/* 
			   tranverse anisotropy from textured derived
			   stiffness tensor
			*/
#define ALL_TRACER_RPHI 10 	/* 2phi/4phi orientations for surface
				   waves  */
#define ALL_TRACER_SAV 11 	/* stiffness tensor */


// individual tracers, see output.c for description of formats
#define LOCATION_TRTP 0		/* location in  t,r,theta,phi formay */
#define LOCATION_LLZ 1		/* location in lon lat z format */
#define LOCATION_XYZ 2		/* location in cartesian format  */
#define LOCATION_TLLZ 3		/* time lon lat z */
#define DEFORMATION_TENSOR 4	/* time deformation tensor F  */
#define STRAIN_COMP 5		/* time L */
#define STRAIN_EVAL 6		/* location and EV(L) */
#define CAUCHY_STRAIN_EVAL 7/* Cauchy tensor strain: time 1/sqrt(eigenvalue(B^1)) eigenvectors format, */
#define STRAIN_EVAL_CART 8	/* cartesian system location EV(L)  */
#define LOCATION_VGM_RTPTVGM 9/* r theta phi t velocity_gradient_matrix */
#define LOCATION_VGM_RTPTVGM_CART 10/* r theta phi t velocity_gradient_matrix in cartesian */
#define LOCATION_VGM_XYZTVGM_CART 11 /* x y z t VGM_cart */
#define XYZ_STRAIN_CART_COMP 12	/* cartesian system location and L */
#define TEXTURE_SAV 13		/* stiffness matrix from texture */
#define TEXTURE_TI 14  	/* best fit transverse anisotropy  */
#define TEXTURE_ODF 15		/* orientation density function  */
#define TEXTURE_RAYLEIGH_RPHI 16		/* 2phi and 4phi terms from Sav for Rayleigh waves  */
#define ISA_AXES 17		/* local ISA axis */
#define TEXTURE_STATS 18		/* statistical values for texture */
#define STRAINRATE_VORTICITY 20	/* strain rate and vorticity */
//
// advect or other general operational modes
// (don't change numbering, this will screw up advect_and_search_for_state)
//
#define FORWARD_TIME 0 // advect forward in time
#define BACKWARD_DEPTH 1 // advect backward in time until tracer arrives from beneath critical depth
#define BACKWARD_STRAIN 2 // advect backward in time until tracer has critical strain
#define BACKWARD_TIME 3 // advect backward in time
#define FORWARD_STRAIN 4 	/* forward compute up to limit strain */
#define CALC_LYAPUNOV 5 // calculate lyapunov exponents
#define ONLY_VEL_STATS 6 // calculate velocity field statistics 
#define ISA 7 // infinite strain axis
#define BACKWARD_TIME_DEPTH 8 // advect backward in time, but only down to certain depth

/* texture development codes */
#define NO_TEXTURE 0
#define KR_TEXTURE 1


/* elastic modul codes */
#define EMOD_CONST 0
#define EMOD_PT 1
#define EMOD_FIXPT 2
#define EMOD_PT_NEW 3
//
// shortcuts
//
#define MEMERROR(x) {fprintf(stderr,"memory allocation error: %s\n",x);exit(-1);}
#define DIFFERENT(x, y) ((fabs((x)-(y))>1.0e-04)?(1):(0))
#define PE(x) {fprintf(stderr,"%s\n",x);}
#define PEE(x) {PE(x);exit(-1);}
#ifndef MAX
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#endif


#define RAND_GEN ran2		/* which random number generator? */

// constant and default values for parameters
#include "constants.h"
// trigonometric constants
#include "trig.h"
// non-dim radius of a level of z km depth
#define ND_RADIUS(x) ((R_E - (x))/R_E)
// depth in km of a non-dimensional radius, z counted positive
#define ZDEPTH(x) (R_E * (1.0-(x)))
//
// function declarations
//
#include "misc.h"
#include "nr_defines.h"
// automatically generated C prototypes, type 'make proto' to generate
#ifndef SINGLE_PRECISION
#include "auto_proto.h"
#else
#include "auto_proto.sgl.h"
#endif
#include "fstrack_ftrn.h"


extern void ellsphere_cart(VPREC *,VPREC *,VPREC *,
			   int *, int *, int *,
			   EXT_FTRN_PREC *,EXT_FTRN_PREC *, 
			   EXT_FTRN_PREC *, EXT_FTRN_PREC *, int *,
			   EXT_FTRN_PREC *, 
			   EXT_FTRN_PREC *,EXT_FTRN_PREC *, EXT_FTRN_PREC *,
			   int *,EXT_FTRN_PREC *,int *,EXT_FTRN_PREC *,
			   EXT_FTRN_PREC *,EXT_FTRN_PREC *,
			   EXT_FTRN_PREC *,int *,EXT_FTRN_PREC *);

extern void ellderivs_cart(EXT_FTRN_PREC *,EXT_FTRN_PREC *,EXT_FTRN_PREC *,
			   int *,
			   VPREC *,VPREC *,VPREC *,
			   int *,int *,int *,EXT_FTRN_PREC *,EXT_FTRN_PREC *,
			   EXT_FTRN_PREC *,EXT_FTRN_PREC *,int *,
			   int *,EXT_FTRN_PREC *);


extern void cart_vel_int(EXT_FTRN_PREC *,EXT_FTRN_PREC *,EXT_FTRN_PREC *,int *,
			 VPREC *,VPREC *,VPREC *,
			 int *,int *,int *,
			 EXT_FTRN_PREC *,EXT_FTRN_PREC *,
			 EXT_FTRN_PREC *,EXT_FTRN_PREC *,EXT_FTRN_PREC *,
			 int *,int *,EXT_FTRN_PREC *);

extern void odeintell_cart(EXT_FTRN_PREC *,int *,EXT_FTRN_PREC *,EXT_FTRN_PREC *,
			   EXT_FTRN_PREC *,EXT_FTRN_PREC *,EXT_FTRN_PREC *,
			   VPREC *,VPREC *,VPREC *,
			   int *,int *,int *,EXT_FTRN_PREC *,EXT_FTRN_PREC *,
			   EXT_FTRN_PREC *,EXT_FTRN_PREC *,int *,int *,
			   EXT_FTRN_PREC *,EXT_FTRN_PREC *,EXT_FTRN_PREC *,
			   EXT_FTRN_PREC *,int *,EXT_FTRN_PREC *) ;

extern void inter_vel(EXT_FTRN_PREC *, EXT_FTRN_PREC *,
		      int *, int *, int *,EXT_FTRN_PREC *,
		      EXT_FTRN_PREC *);

extern void invert3x3(EXT_FTRN_PREC *,EXT_FTRN_PREC *);

extern void svdcmp(EXT_FTRN_PREC *a,int *,int *,int *,int *,
		   EXT_FTRN_PREC *,EXT_FTRN_PREC *,EXT_FTRN_PREC *);

extern void gramschmidt(EXT_FTRN_PREC *,EXT_FTRN_PREC *);

//
// Eispack eigensystem stuff
extern void rs_(int *, int *, double *,double *, int *, 
		double *, double *,double *, int *);
extern void s_rs_(int *, int *, float *,float *, int *, 
		 float *, float *,float *, int *);

extern void rg_(int *, int *, double *,double *, double *,
		int *, double *, int *,double *, int *);
extern void s_rg_(int *, int *, float *,float *, float *,
		  int *, float *, int *,float *, int *);
//
// EXPOKIT
extern void dgpadm(int *,int *,EXT_FTRN_PREC *,EXT_FTRN_PREC *,
		   int *,EXT_FTRN_PREC *,int *,int *,int *,
		   int *,int *);


#define __FSTRACK_HEADER_SOURCED__
#endif
