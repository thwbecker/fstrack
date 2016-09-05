/*


 STRUCTURES for fstrack program 

 $Id: structures.h,v 1.18 2006/06/15 20:23:06 becker Exp becker $

*/
#include "precision.h"

/* 
   
surface wave sensitivity kernel

*/
struct swss_layer{
  COMP_PRECISION r,RA,RF,RL;	/* depth and parameters for SW sens layer, Rayleigh parmaters */
  COMP_PRECISION LL, LN;	/* Love */
};
struct swss{
  COMP_PRECISION p;		/* period */
  int n;			/* number of depth layers */
  struct swss_layer *l;		/* one layer */

};
/*

  (deformation) state (of-a-tracer) structure
  
  holds the strain matrix s, which is (normally) the squared
  left-stretch tensor, L, L^2,
  the position x[0..2] of the tracer at time t
  and the deformation matrix F, which is in x[3....8]

  IF YOU ADD VARIABLES HERE, MAKE SURE TO ADJUST HANDLE_STATE ROUTINES SUCH 
  AS COPY_STATE

*/
struct stt{
  COMP_PRECISION left_stretch[9];/* 
				    L^2, the squared left stretch strain matrix or 
				    Lyapunov exponents
				 */
  COMP_PRECISION x[12];/* 
			  position vector plus deformation matrix, f[i]=x[i+3]
			  x: 0 1 2 3 4 5 6 7 8 9 10 11
			  r: 0 1 2                      r vector
			  f:       0 1 2 3 4 5 6  7  8  deformation
			           rr      tt        pp
		       */

  COMP_PRECISION t; // time
  //
  // for ISA calculation: grain orientation lag parameter
  COMP_PRECISION pipar;
  //
  // ISA vector
  COMP_PRECISION isa[3];
  // does ISA exist?
  my_boolean isa_exists;
  /* 
     texture related
  */
  COMP_PRECISION Sav[36];	/* average stiffness tensor from
				   texture computation, stored in C
				   style (doesn't matter, since
				   symmetric)
				*/
  int npole;			/* tracer pole density size counter */
  COMP_PRECISION *pdens;		/* pole figure values  */
#define MAX_OL_AXES_ENTRIES 9
  COMP_PRECISION olivine_axes_stats[MAX_OL_AXES_ENTRIES]; /* olivine axes stats,
							     unweighted and weighted */
};

/*

  TRACER structure, holds nstate number of states

*/
#define NR_T_ATTR 1
struct trc{
  struct stt *state; // states, ie. strain and position at certain time
  unsigned short int nstate;// number of different states
  my_boolean discarded; // some tracers can be discarded for various reasons
  COMP_PRECISION attr[NR_T_ATTR];// scalar tracer attribute(s)
};


/*

  depth weight structure, used by average_tracers

*/
struct dw{
  int n; // nr of weights
  COMP_PRECISION *z,*w;/* depth where weight is specified and weight at 
			  that depth */
};
/*

strain structure, used by average_tracers

*/
// definitions for strain averaging arrays
#define STR_N_ESI 5 // nr of fields
#define STR_ESI_S_H 0 // the horizontal component of e1/e2 times the e1 vector
#define STR_ESI_AZI 1  // the azimuth of that componennt (degrees)
#define STR_ESI_LRR 2  // L_rr
#define STR_ESI_AX 3   // the x part of the horizontal projection 
#define STR_ESI_AY 4   // the y part 
struct str{
  COMP_PRECISION x[3],l[9],age,ws,wss;/* position: x 
					 L and age at
					 each depth, sum of weights,
					 and strain weighted sum of weights */
  COMP_PRECISION es[STR_N_ESI];/* log(e1/e2)_h azi[deg], 
				  lrr, azi_x, azi_y at each depth 
			   */
  COMP_PRECISION dazi,oldazi; // accumulated absolute rotation and old azimuth
  COMP_PRECISION pipar;/* grain orientation lag GOL factor */
};


/*
  
 MODEL STRUCTURE

 holds pointers to the velocities at certain depth levels, the tracers
 and other control variables

*/
#define MAX_NR_TINIT_LEVELS 30
struct mod{
  //
  // tracer related
  //
  int tinitmode;// tracer init mode
  unsigned int ntracer; // number of tracers
  struct trc *tracer;// tracer structures
  char tlatfilename[STRLEN];// filename with lateral tracer distribution
  // depth levels
  int ntlevels;// nr of levels
  COMP_PRECISION *tlevel,init_tracer_depth;// depth levels, single layer 
  char tdepthfilename[STRLEN];// depth level filename
  unsigned int *ntd;// nr of tracers at this depth
  my_boolean tdepth_from_file;
  int ntout;// nr of tracers that are output 
  //
  // advect mode
  int amode;
 
  // nsteps of advection for tracer history output
  int nsteps;
  //
  // bailout values for advect
  //
  COMP_PRECISION sdepth;// target depth 
  COMP_PRECISION maxstrain;// 
  COMP_PRECISION tf;// final time, forward or backward
  COMP_PRECISION itime;// intial time for forward advection, usually zero
  COMP_PRECISION bstime;// time for bailout when tracing strains backward

  my_boolean use_initial_strain; /* this is normally off, 
				    start with unity matrix */
  COMP_PRECISION initial_strain[9]; /* possible initial strain in spherical 
				       system */
  /*

  infinite strain axis related: 
  critical pi parameter, if larger, ISA doesn't exist
  */
  COMP_PRECISION pi_crit;
  /* 

  load the surface wave sensitivity functions?
     
  */
  my_boolean use_sw_sens;
  my_boolean sw_sens_init;
  struct swss *swsens;
  char sw_sens_file[STRLEN];	/* starting filename */
  // IO related
  my_boolean read_gmt;
  my_boolean verbose;
  char ostring[STRLEN];		/* suffix for output files */
  my_boolean sav_out;		/* output of stiffness tensor for texture */
  /* 
     derivative parameters  
  */
  struct der_par *dp;
};
/* shortcut to access the derivative related parameters */
#define MDP model->dp
/* shortcut to access the DREX parameters */
#define DREX model->dp->drex

