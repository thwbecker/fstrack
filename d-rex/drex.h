#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>

#ifndef __DREX_HRD_READ__

/* size of pole figure storage */
#define DREX_POLE_NX 120
#define DREX_POLE_NY 30

#include "precision.h"

/* 

   differerent conventions for Cartesian coordinate systems

*/
/* DREX_REG2SW_CONVENTION:   x=North, y=East, z=down, for surface wave (sens_handling) */
/* DREX_REG_CART_CONVENTION: x=South, y=East, z=up  */
/* DREX_REG2SW_CONVENTION:   go from regular convention to surface wave convention */
/* this is in separate file to allow sourcing from C and fortran */
#include "drex_coord_convention.h"



#include "drex_fconst.h"
#define DREX_PI_C 3.1415926535897932


#ifndef MEMERROR
#define MEMERROR(x) {fprintf(stderr,"memory allocation error: %s\n",x);exit(-1);}
#endif
#ifndef PE
#define PE(x) {fprintf(stderr,"%s\n",x);}
#endif
#ifndef PEE
#define PEE(x) {PE(x);exit(-1);}
#endif

struct drex_para{
  int size,size3,bsize;		/* size: total number of grains 
				   size3: number of grains in each direction
				   bsize: size of all time dependent arrays,
				   12 + 2 * size + (size * 9 ) * 2


				*/
  COMP_PRECISION lambda, Mob, chi;
  // lambda = nucleation parameter
  // Mob = grain mobility
  // chi = threshold volume fraction for activation of grain boundary sliding
  COMP_PRECISION Xol;		/* fraction of olivine */
  COMP_PRECISION tau[4],tau_ens; /*  */
  /* 

  NOTE THAT ALL ARRAYS WILL BE addressed FORTRAN STYLE!

  */
  COMP_PRECISION alt[27];	/* epsilon ijk tensor */
  COMP_PRECISION stressexp;	/*  */
  COMP_PRECISION del[9];	/* kronecker */
  int ijkl[9];			/* indices for eijkl */
  int l1[6],l2[6];			/* tensor for Sij converions */
  COMP_PRECISION S0[36],S0_ens[36]; /* elastic tensors */

  
  COMP_PRECISION *x;		/* this is the master array for odf, odf_ens,
				   acs, and acs_ens 

				*/
  /* integer pointers to locations with an x-type array, those might
     or might not point to the locations pointed to by odf, odf_ens
     etc */
  unsigned int podf,podf_ens,pacs,pacs_ens;
  COMP_PRECISION *odf,*odf_ens;	/* those will just point to locations within x */
  COMP_PRECISION *acs,*acs_ens;
  /* 
     average a axes
  */
  COMP_PRECISION avg_axes[6], 	/* avg a axis, and odf-weighted avg a
				   axes */
    avg_axes_ens[6],dev_axes[2],
    dev_axes_ens[2];		/* deviations, unweighted and weighted */
  COMP_PRECISION olivine_jindex; /* j index */
  /* random initial arrays */
  COMP_PRECISION *acs0,*acs0_ens;
  /* for pole figures */
  int np[2],npole;
  COMP_PRECISION *pdens;
  /*  */
  COMP_PRECISION *rt,*rt_ens;
  /* save the gamma strain factors? */
  int isave_gamma;
  float *gamma_save;

  my_boolean initialized;
};
/* 

drex fabric codes

 */
#define DREX_FABRIC_FREE -1
#define DREX_FABRIC_ATYPE 0
#define DREX_FABRIC_BTYPE 1 
#define DREX_FABRIC_CTYPE 2
#define DREX_FABRIC_DTYPE 3 
#define DREX_FABRIC_ETYPE 4
#define DREX_FABRIC_HIGHP 5
/* switch between A and HP slip system */
#define DREX_FABRIC_SWITCH_A_HIGHP 6
#define DREX_SS_PSWITCH 10 	/* default transition pressure, 
				   should be < 11 Gpa 
				*/
/* 
   Kaminski & Ribe texture parameters

*/
#define DREX_SIZE3_DEF 12	/* number of grains per direction */
#define DREX_XOL_DEF 70.0	/* olivine fraction in percent */
#define DREX_TAU_OLIVINE_DEF_ATYPE {1.0,2.0,3.0,1e60}	/* normal slip
							   system for
							   olivine,
							   type A */
#define DREX_TAU_OLIVINE_DEF_BTYPE {3.0,2.0,1.0,1e60}	/* B type slip
							   system for
							   olivine */
#define DREX_TAU_OLIVINE_DEF_CTYPE {3.0,1e60,2.0,1.0}	/* C type slip
							   system for
							   olivine */
#define DREX_TAU_OLIVINE_DEF_DTYPE {1.0,1.0,3.0,1e60}	/* D type slip
							   system for
							   olivine */
#define DREX_TAU_OLIVINE_DEF_ETYPE {2.0,1.0,1e60,3.0}	/* E type slip
							   system for
							   olivine */
#define DREX_TAU_OLIVINE_DEF_HIGHP {3.0,1e60,2.0,1.0}	/* high pressure type slip
							   system for olivine, experimental,
							*/
#define DREX_TAU_ENS_DEF 1.0	/* enstatite factor for slip system */
#define DREX_MOB_DEF 125.0	/* grain boundary mobility, 
				   ~125 */
#define DREX_CHI_DEF 0.30	/* 
				   chi (threshold volume fraction for
				   GBS) this should be ~0.3
				*/
#define DREX_LAMBDA_DEF 5.0	/* this shouhld be around 5, the grain nucleation factor */




/* individual elastic constant modes */
#define DREX_CTP_OL_ESTEY 0
#define DREX_CTP_EN_ESTEY 1
#define DREX_CTP_OL_NEW 6
#define DREX_CTP_EN_NEW 7
#define DREX_C_OL_KR 2
#define DREX_C_EN_KR 3
#define DREX_C_OL_JULES 4
#define DREX_C_EN_JULES 5
/* combined modes modes */
#define DREX_CTP_CONST_KR 0		/* constant, KR */
#define DREX_CTP_ESTEY 1	/* p,T from Estey */
#define DREX_CTP_NEW 2		/* newer sources */

/* 


routine to be called from C

*/
void drex_initialize(struct drex_para **,int ,int , 
		     COMP_PRECISION ,COMP_PRECISION *, 
		     COMP_PRECISION ,COMP_PRECISION , 
		     COMP_PRECISION,COMP_PRECISION, int *, my_boolean,
		     my_boolean);
void drex_isacalc_from_norm_vgmc(COMP_PRECISION *, COMP_PRECISION *,
				 COMP_PRECISION *, COMP_PRECISION *);
void drex_init_pathline(struct drex_para *);
void drex_deriv(COMP_PRECISION *,COMP_PRECISION *,COMP_PRECISION *, 
		struct drex_para *);
void drex_check_phys_limits(COMP_PRECISION *,int ,struct drex_para *);
void drex_vhr(struct drex_para *, COMP_PRECISION *,int,COMP_PRECISION,
	      COMP_PRECISION,COMP_PRECISION);
void drex_decsym(COMP_PRECISION *,COMP_PRECISION *,COMP_PRECISION *,COMP_PRECISION *,
		 COMP_PRECISION *,COMP_PRECISION *,COMP_PRECISION *,COMP_PRECISION *,
		 COMP_PRECISION *,COMP_PRECISION *,COMP_PRECISION *,COMP_PRECISION *,
		 COMP_PRECISION *,COMP_PRECISION *,int *);
void drex_compute_pdens(struct drex_para *,COMP_PRECISION *,
			COMP_PRECISION,my_boolean,my_boolean,int , my_boolean);
my_boolean drex_is_symmetric(COMP_PRECISION *, int );
void drex_test_if_initialized(struct drex_para *,char *);
void drex_compute_avg_axes(struct drex_para *,COMP_PRECISION *,my_boolean);
void drex_get_sav_constants(COMP_PRECISION *, int , COMP_PRECISION ,COMP_PRECISION );
void drex_compute_symm_tensors(COMP_PRECISION *, int , COMP_PRECISION *,COMP_PRECISION *,COMP_PRECISION *,
			       COMP_PRECISION *,COMP_PRECISION *,char *);
void drex_compute_symten_frac(COMP_PRECISION *,COMP_PRECISION *,COMP_PRECISION *,COMP_PRECISION *);
COMP_PRECISION drex_ctn_36(COMP_PRECISION *);
void drex_voigt_avg(COMP_PRECISION *,COMP_PRECISION *,COMP_PRECISION ,
		    COMP_PRECISION *);
void drex_calc_euler_rad_from_xp(COMP_PRECISION *,COMP_PRECISION *,
				 COMP_PRECISION *, COMP_PRECISION *,
				 int *);
void drex_calc_euler_rad_from_xp_(COMP_PRECISION *,COMP_PRECISION *,
				  COMP_PRECISION *, COMP_PRECISION *,
				  int *);
void drex_assign_tau_from_system(COMP_PRECISION *, int type, 
				 COMP_PRECISION *, COMP_PRECISION, COMP_PRECISION);


/* 
   those are in misc_c in this directory, or should function like those 
   routines if linked from somewhere else
*/
void my_vecalloc(COMP_PRECISION **,int ,char *);
void my_svecalloc(float **,int ,char *);
void my_vecrealloc(COMP_PRECISION **,int ,char *);
void a_equals_b_vector(COMP_PRECISION *,COMP_PRECISION *,int);
void a_is_b_transpose_3by3(COMP_PRECISION *,COMP_PRECISION *);
void transpose_3by3_inplace(COMP_PRECISION *);
void a_is_b_transpose_nbyn(COMP_PRECISION *,COMP_PRECISION *, int);
void transpose_nbyn_inplace(COMP_PRECISION *,int );
void sub_a_from_b_vector(COMP_PRECISION *, COMP_PRECISION *,int );
COMP_PRECISION mean(COMP_PRECISION *,int);
COMP_PRECISION norm(COMP_PRECISION *,int);
COMP_PRECISION vecdotp(COMP_PRECISION *,COMP_PRECISION *,int);
void print_vector(COMP_PRECISION *,int ,FILE *);

/* 

deal with FORTRAN calling conventions

*/
#ifndef NO_UNDERSCORE
#define drex_init_para drex_init_para_
#define drex_vhr_avg_ftrn drex_vhr_avg_ftrn_
#define drex_init_acs_random_ftrn drex_init_acs_random_ftrn_
#define drex_isacalc_ftrn drex_isacalc_ftrn_
#define drex_init_pathline_ftrn drex_init_pathline_ftrn_
#define drex_deriv_ftrn drex_deriv_ftrn_
#define drex_strain_rate_ftrn drex_strain_rate_ftrn_
#define drex_check_phys_limits_ftrn drex_check_phys_limits_ftrn_
#define drex_decsym_ftrn drex_decsym_ftrn_
#define drex_compute_pdens_ftrn drex_compute_pdens_ftrn_
#define drex_compute_swpar_ftrn drex_compute_swpar_ftrn_
#define drex_compute_avg_axes_ftrn drex_compute_avg_axes_ftrn_
#define drex_rotate_6x6_rad_ftrn drex_rotate_6x6_rad_ftrn_
#define drex_rotate_6x6_deg_ftrn drex_rotate_6x6_deg_ftrn_
#define drex_rotate_6x6_rot_ftrn drex_rotate_6x6_rot_ftrn_
#define drex_get_olivine_sav_ftrn drex_get_olivine_sav_ftrn_
#define drex_get_olivine_new_sav_ftrn drex_get_olivine_new_sav_ftrn_
#define drex_get_olivine0_sav_ftrn drex_get_olivine0_sav_ftrn_
#define drex_get_olivine1_sav_ftrn drex_get_olivine1_sav_ftrn_
#define drex_get_enstatite_sav_ftrn drex_get_enstatite_sav_ftrn_
#define drex_get_enstatite_new_sav_ftrn drex_get_enstatite_new_sav_ftrn_
#define drex_get_enstatite0_sav_ftrn drex_get_enstatite0_sav_ftrn_
#define drex_get_enstatite1_sav_ftrn drex_get_enstatite1_sav_ftrn_
#define drex_fullsym6_ftrn drex_fullsym6_ftrn_
#define drex_st_vhr_avg_ftrn drex_st_vhr_avg_ftrn_
#define drex_scca drex_scca_
#define drex_v21d_ftrn drex_v21d_ftrn_
#define drex_d21v_ftrn drex_d21v_ftrn_
#define drex_projec_ti drex_projec_ti_
#define drex_projec_mo drex_projec_mo_
#define drex_projec_or drex_projec_or_
#define drex_projec_te drex_projec_te_
#define drex_projec_is drex_projec_is_
#define drex_tens4_norm2_ftrn  drex_tens4_norm2_ftrn_
#define drex_tens4_ftrn drex_tens4_ftrn_
#define drex_rot4 drex_rot4_
#define drex_calc_rotmat_cart drex_calc_rotmat_cart_
#endif

/* 

declarations for F90 subroutines from D-REX 

*/
/* drex_util */
extern void drex_check_phys_limits_ftrn(double *,double *,double *,double * , 
					int *);
extern void drex_init_para(int *,int *,double *,double *,
			   double *,double *,double *,
			   double *,double *,double *,
			   double *,int *,int *,int *,
			   double *,double *,int *);
extern void drex_init_acs_random_ftrn(int *,int *,double *,int *);

extern void drex_vhr_avg_ftrn(double *,double *,double *,
			      double *,double *,double *,
			      double *,int *,int *,int *,int *,
			      double *,double *);
extern void drex_compute_pdens_ftrn(double *,double *,double *,double *,int *,
				    double *,double *,int *, int *,double *,
				    double *,int *,int *,int *,int *);
extern void drex_compute_avg_axes_ftrn(double *,double *,double *,double *,
				       int *,double *,int *,double *,double *,
				       double *,double *,double *);
extern void drex_fullsym6_ftrn(double *);
extern void drex_v21d_ftrn(double *,double *);
extern void drex_d21v_ftrn(double *,double *);
extern void drex_projec_ti(double *,double *, double *);
extern void drex_projec_mo(double *,double *, double *);
extern void drex_projec_or(double *,double *, double *);
extern void drex_projec_te(double *,double *, double *);
extern void drex_projec_is(double *,double *, double *,double *, double *);
extern void drex_st_vhr_avg_ftrn(double *,double *,double *,double *,double *);
/* drex_elast_const */
extern void drex_get_olivine_sav_ftrn(double *,double *, double *);
extern void drex_get_olivine_new_sav_ftrn(double *,double *, double *);
extern void drex_get_olivine0_sav_ftrn(double *);
extern void drex_get_olivine1_sav_ftrn(double *);
extern void drex_get_enstatite_sav_ftrn(double *,double *, double *);
extern void drex_get_enstatite_new_sav_ftrn(double *,double *, double *);
extern void drex_get_enstatite0_sav_ftrn(double *);
extern void drex_get_enstatite1_sav_ftrn(double *);

/* drex_deriv */
extern void drex_deriv_ftrn(double *,double *,
			    double *,int *,double *,double *, 
			    double *,double *,double *,double *,
			    double *,double *,double *,double *,
			    double *,double *,double *,
			    double *,double *,double *,double *,double *,
			    int *, float *);

extern void drex_strain_rate_ftrn(double *,double *,double *);

/* drec_decmod */
extern void drex_decsym_ftrn(double *,double *,double *,double *,double *,double *,
			     double *,double *,double *,double *,double *,double *,
			     double *,double *, int *);

extern void drex_isacalc_ftrn(double *, double *, double *,double *);
extern void drex_init_pathline_ftrn(double *, double *,double *, double *,int *);

extern void drex_compute_swpar_ftrn(double *,double *,double *,
				    double *,double *,double *,
				    double *, double *,double *,
				    double *, double *);
extern void drex_rotate_6x6_rad_ftrn(double *,double *,double *,double *,double *);
extern void drex_rotate_6x6_deg_ftrn(double *,double *,double *,double *,double *);
extern void drex_rotate_6x6_rot_ftrn(double *,double *,double *);

extern void drex_tens4_norm2_ftrn(double *, double *);
extern void drex_tens4_ftrn(double *,double *);
extern void drex_scca(double *,double *,double *,double *,double *,double *,double *,double *,double *);
extern void drex_rot4(double *,double *,double *);
extern void drex_calc_rotmat_cart(double *,double *,double *,double *);

#define __DREX_HRD_READ__
#endif
