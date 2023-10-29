#include "drex.h"
/* 
   check is Sav really get's passed as symmetric 
*/

//#define DREX_DEBUG_SYMMETRY


/* 

driver routines for Kaminski's DREX code

Thorsten Becker (thwbecker@post.harvard.edu)

$Id: drex_driver.c,v 1.16 2011/04/12 06:18:44 becker Exp becker $


*/

/* 
   
initialize the parameters

imode = 0: use input parameters
imode = 1: read parameters from file
imode = 2: use standard parameters

odf_np[2]: dimensions of pole save ODF array
*/
#define DREX (*drex)

void drex_initialize(struct drex_para **drex,int imode,
		     int size3, /* nr_grains^(1/3) */
		     COMP_PRECISION Xol, /* ~70 (in %) */
		     COMP_PRECISION *tau, /* 1,2,3,1e60 */
		     COMP_PRECISION tau_ens, /* 1 */
		     COMP_PRECISION Mob, /* 125 */
		     COMP_PRECISION chi, /* 0.3 */
		     COMP_PRECISION lambda, /* 5.0 */
		     int *odf_np,
		     my_boolean use_oriented_grains,
		     my_boolean save_gamma_factors)
{
  static my_boolean init = FALSE;
  int size9,i;
  int ran_irmode = 1;		/* random number generator mode */
  int orient_irmode = 2;	/* oriented grains mode */
  if(!init){
    *drex = (struct drex_para *)calloc(1,sizeof(struct drex_para));
    if(!*drex)
      MEMERROR("drex_initialize: 1");
    init = TRUE;
  }else
    PEE("drex_initialize: error: should only be called once");
  switch(imode){
  case 0:
    /* parameters are input */
    PE("drex_initialize: using values as input:");
    fprintf(stderr,"drex_initialize: size3: %i Xol: %g%%\n",
	    size3,Xol);
    fprintf(stderr,"drex_initialize: tau[4]: %g, %g, %g, %g tau_ens: %g\n",
	    tau[0],tau[1],tau[2],tau[3],tau_ens);
    fprintf(stderr,"drex_initialize: Mob: %g chi: %g lambda: %g\n",
	    Mob,chi,lambda);
    /* 
       
    make sure that all input parameters are assigned to the structure

    */

    DREX->size3 = size3;
    DREX->Xol = Xol;
    a_equals_b_vector(DREX->tau,tau,4);
    DREX->tau_ens = tau_ens;
    DREX->Mob = Mob;
    DREX->chi = chi;
    DREX->lambda = lambda;
    break;
  case 1:
    /* reading from file */
    PE("drex_initialize: reading parameters from file");
    break;
  case 2:
    PE("drex_initialize: using default parameters");
    break;
  default:
    fprintf(stderr,"drex_initialize: error: imode %i undefined\n",
	    imode);
    exit(-1);
    break;
  }
  /* 
     initialize the main parameters 
  */
  drex_init_para(&DREX->size,&DREX->size3,&DREX->Xol,
		 DREX->tau,&DREX->tau_ens,&DREX->Mob,&DREX->chi,
		 &DREX->lambda,&DREX->stressexp,
		 DREX->alt,DREX->del,DREX->ijkl,DREX->l1,DREX->l2,
		 DREX->S0,DREX->S0_ens,&imode);
  
  /* 
     allocate the grain parameters arrays
  */
  size9 = DREX->size * 9;
  /* 
     this is the size for 12 + odf, odf_ens, acs, and acs_ens 
  */
  DREX->bsize = 12 + DREX->size * 2 + size9 * 2;
  my_vecalloc(&DREX->x, DREX->bsize,"drex_init: x");
  /* 
     assign integer pointers to allow addressing within and x-type array
     (needn't be *the* x)
  */
  i = 12;         DREX->podf=i;		/* odf */
  i += DREX->size;DREX->podf_ens=i; /* odf_ens */
  i += DREX->size;DREX->pacs=i; /* acs */
  i += size9;     DREX->pacs_ens=i;	/* acs_ens */
  /* 
     assign locations in the real x to other pointers 
  */
  DREX->odf =     (DREX->x + DREX->podf);
  DREX->odf_ens = (DREX->x + DREX->podf_ens);
  DREX->acs     = (DREX->x + DREX->pacs);
  DREX->acs_ens = (DREX->x + DREX->pacs_ens);
  /* other arrays */
  my_vecalloc(&DREX->rt,DREX->size,"drex_init: 2");
  my_vecalloc(&DREX->rt_ens,DREX->size,"drex_init: 2");
  /* those are only need to initialize */
  my_vecalloc(&DREX->acs0,size9,"drex_init: 3");
  my_vecalloc(&DREX->acs0_ens,size9,"drex_init: 3");
  /* 

  initialize the original, random initial distribution for each path
  this will be reused each time a new pathline gets initialized
  (will be initialized FORTRAN style)

  */
  if(use_oriented_grains){
    fprintf(stderr,"drex_init: WARNING: initializing oriented grains\n");
    drex_init_acs_random_ftrn(&DREX->size,&DREX->size3,DREX->acs0,&orient_irmode);
    drex_init_acs_random_ftrn(&DREX->size,&DREX->size3,DREX->acs0_ens,&orient_irmode);
  }else{			/* regular random */
    drex_init_acs_random_ftrn(&DREX->size,&DREX->size3,DREX->acs0,&ran_irmode);
    drex_init_acs_random_ftrn(&DREX->size,&DREX->size3,DREX->acs0_ens,&ran_irmode);
  }
  if(save_gamma_factors){	/* make room for saving gamma strain rate factors */
    DREX->isave_gamma = 1;
    my_svecalloc(&DREX->gamma_save,DREX->size*3,"drex_init: 3b");
  }else{
    DREX->isave_gamma = 0;
  }
  /* 
     
  initialize the array for pole figures

  */
  for(i=0;i<2;i++)
    DREX->np[i] = odf_np[i];
  /* 
     density entries for three axes and two components 
  */
  DREX->npole = DREX->np[0] * DREX->np[1] * 3 * 2;
  my_vecalloc(&DREX->pdens,DREX->npole,"drex_init: 4");


  DREX->initialized = TRUE;
  fprintf(stderr,"drex_initialize: done, using %i grains\n",
	  DREX->size);
}

#undef DREX
/* 

intialize the grain arrays for each pathline by assigning the acs0 arrays
and setting an even orientation density function

*/
void drex_init_pathline(struct drex_para *drex)
{
  int i;
#ifdef DREX_DEBUG
  drex_test_if_initialized(drex,"drex_init_pathline");
#endif
  for(i=0;i<drex->npole;i++)
    drex->pdens[i] = 0.0;	/* this isn't really necessary */
  /* 
     olivine component 
  */
  drex_init_pathline_ftrn(drex->rt,drex->odf,drex->acs,
			  drex->acs0,&drex->size);
  /* 
     enstatite component 
  */
  drex_init_pathline_ftrn(drex->rt_ens,drex->odf_ens,
			  drex->acs_ens,drex->acs0_ens,
			  &drex->size);
}

/* 
   
driver for the FORTRAN routine to compute the ISA axis 
(GOL only gets modified if F==0)

vgm_c is input as the NORMALIZED velocity gradient matrix in C fashion

routine returns the ISA vector and gol, the latter is used a parameter
for the existence of the ISA

if the GOL is undefined because of F==0, then isa=(0,0,0) and GOL=-1

if the ISA is the velocity vector, then ISA=(-1,-1,-1)

omega is the rotation rate
*/
void drex_isacalc_from_norm_vgmc(COMP_PRECISION *isa, 
				 COMP_PRECISION *gol,
				 COMP_PRECISION *vgm_c,
				 COMP_PRECISION *omega)
{
  COMP_PRECISION vgm_f[9];
  /* 
     resort the matrix , G_FTRN = G_C^T
  */
  a_is_b_transpose_3by3(vgm_f,vgm_c);
  /* 
     call the FORTRAN based routine 
  */
  drex_isacalc_ftrn(isa,gol,vgm_f,omega);
}


/* 
   
driver for the derivative routines. this is pretty much identical to a call 
to drex_deriv_ftrn, only that the velocity gradient matrix is input C style 
as vgm_c . we have also rephrased all timevariing variables (odf,odf_ens,acs,acs_ens) 
and their derivatices as x[] and dx[]. they are both dimensionsed drex->bsize

*/
void drex_deriv(COMP_PRECISION *vgm_c, COMP_PRECISION *x, 
		COMP_PRECISION *dx, struct drex_para *drex)

{

  COMP_PRECISION vgm_f[9],e[9],epsnot;
#ifdef DREX_DEBUG
  drex_test_if_initialized(drex,"drex_deriv");
#endif
  /* flip the velocity gradient matrix since internally, we are using FORTRAN
     convention */
  a_is_b_transpose_3by3(vgm_f,vgm_c);
  /* 
     compute the strain rate matrix
  */
  drex_strain_rate_ftrn(vgm_f,e,&epsnot);
  /* 
     
  compute the derivatives for the arrays, use addressing within x (and
  not out internal pointers to ACS and such), since x might be some general 
  array

  */
  drex_deriv_ftrn(vgm_f,e,&epsnot,&drex->size,
		  (x+drex->pacs),(x+drex->pacs_ens),
		  (dx+drex->pacs),(dx+drex->pacs_ens),
		  (x+drex->podf),(x+drex->podf_ens),
		  (dx+drex->podf),(dx+drex->podf_ens),
		  drex->rt,drex->rt_ens,
		  drex->tau,
		  &drex->tau_ens,&drex->stressexp,drex->alt,
		  &drex->Xol,&drex->lambda,&drex->Mob,&drex->chi,
		  &drex->isave_gamma,drex->gamma_save);
}
/* 

limit acs and odf arrays to physical limits

*/
void drex_check_phys_limits(COMP_PRECISION *x,int n,
			    struct drex_para *drex)
{
#ifdef DREX_DEBUG
  int s9,sc;
  /* check if the DREX structure was initialized */
  drex_test_if_initialized(drex,"drex_check_phys_limits");
  /* check if the pointers work out OK */
  sc = 12;s9 = drex->size * 9;
  sc  = 12;    
  if(drex->podf != sc)
    PEE("drex_check_phys_limits: pointer error: podf");
  sc += drex->size;
  if(drex->podf_ens != sc)
    PEE("drex_check_phys_limits: pointer error: podf_ens");
  sc += drex->size;
  if(drex->pacs != sc)
    PEE("drex_check_phys_limits: pointer error: pacs");
  sc += s9;        
  if(drex->pacs_ens != sc)
    PEE("drex_check_phys_limits: pointer error: pacs_ens");
  sc += s9;
  if(n != sc)
    PEE("drex_check_phys_limits: pointer error: n");
#endif
  drex_check_phys_limits_ftrn((x+drex->pacs),(x+drex->pacs_ens),
			      (x+drex->podf),(x+drex->podf_ens),
			      &drex->size);
}
/* 

compute the average stiffness matrix Sav[36] in FORTRAN sorting, sav
is output aggregate using Voigt (f=0), Reuss (f=1), or VHR (f=0.5)
averaging

*/
void drex_vhr(struct drex_para *drex, COMP_PRECISION *sav,
	      int var_mod, 
	      COMP_PRECISION temp, COMP_PRECISION pressure,
	      COMP_PRECISION f)
{
  COMP_PRECISION sol[36],sen[36];
  switch(var_mod){
  case DREX_CTP_ESTEY:
    /* 
       use varying moduli from Estey
    */
    drex_get_olivine_sav_ftrn(sol,&temp,&pressure);
    drex_get_enstatite_sav_ftrn(sen,&temp,&pressure);
    drex_vhr_avg_ftrn(sol,sen,drex->acs,drex->acs_ens,
		      drex->odf,drex->odf_ens,
		      &drex->Xol,&drex->size,
		      drex->ijkl,drex->l1,drex->l2,sav,&f);
    break;
  case DREX_CTP_CONST_KR:
    /* use the tensors that were initialized once in the beginning */
    drex_vhr_avg_ftrn(drex->S0,drex->S0_ens,drex->acs,
		      drex->acs_ens,drex->odf,drex->odf_ens,
		      &drex->Xol,&drex->size,drex->ijkl,
		      drex->l1,drex->l2,sav,&f);
    break;
  case DREX_CTP_NEW:
    /* 
       use varying moduli from newer sources
    */
    drex_get_olivine_new_sav_ftrn(sol,&temp,&pressure);
    drex_get_enstatite_new_sav_ftrn(sen,&temp,&pressure);
    drex_vhr_avg_ftrn(sol,sen,drex->acs,drex->acs_ens,
		      drex->odf,drex->odf_ens,
		      &drex->Xol,&drex->size,
		      drex->ijkl,drex->l1,drex->l2,sav,&f);
    break;
  default:
    fprintf(stderr,"drex_vhr: var_mod mode %i undefined\n",
	    var_mod);exit(-1);
    break;
  }
  /* we don't have to back-rotate Sav, since it's symmetric */

#ifdef DREX_DEBUG_SYMMETRY
  if(!drex_is_symmetric(sav,6)){
    fprintf(stderr,"drex_vhr: error: Sav not symmetric\n");
    exit(-1);
  }
#endif
}
/* 
   
compute the decomposition into transverse isotropic symmetry from a 
stiffness matrix sav[36] which is input in C ctyle

input:
ced[6,6] tensor


for finding best axes:
scca_old_mode: 1: only find best-fit, and align (old mode)
               0: find best and worst fit, and use those as a coordinate system

output:

k,g: elastic moduli
vel[9]:  [0,1]vpmax/min*sqrt(dens), [2]: c_33^0 [3]: v_44^0
                 [4]: epsilon*c_33^0 [5]: gamma*c_44^0 [6]: delta*c_33^0
		 [7,8] vs1*sqrt(dens),vs2*sqrt(dens)
tiaxis[6]: heaxgonal symmetry axis, [0..2]: best [3..5]: worst
symm_frac[6]: ! fractional norms for:
               ! 1: isotropy
	       ! 2: anisotropy due to VTI
	       ! 3: anisotropy due to tet
	       ! 4: anisotropy due to ort
	       ! 5: anisotropy due to mon
	       ! 6: anisotropy due to tri

dc_...: decomposed tensors[6,6] ALL in SCC system

ced_scc: rotated original tensor

scc_rot: inverse rotation matrix [3,3] to go from SCC back to orig

*/
void drex_decsym(COMP_PRECISION *ced, COMP_PRECISION *k, COMP_PRECISION *g, 
		 COMP_PRECISION *vel,COMP_PRECISION *symm_frac, 
		 COMP_PRECISION *tiaxis,
		 COMP_PRECISION *dc_iso, COMP_PRECISION *dc_hex, 
		 COMP_PRECISION *dc_tet, COMP_PRECISION *dc_ort, 
		 COMP_PRECISION *dc_mon, COMP_PRECISION *dc_tri,
		 COMP_PRECISION *ced_scc, 
		 COMP_PRECISION *scc_irot,
		 int *scca_old_mode)

{
  COMP_PRECISION ccp[36],loc_rot[9];
  int i,j;
#ifdef DREX_DEBUG_SYMMETRY
  if(!drex_is_symmetric(ced,6)){
    fprintf(stderr,"drex_decsym: error: C not symmetric\n");
    exit(-1);
  }
#endif
  a_equals_b_vector(ccp,ced,36);
  /* 
     call the D-REX routines, doesn't matter that Sav is C style,
     since it's symmetric
  */
  drex_decsym_ftrn(ccp,k,g,symm_frac,tiaxis,
		   dc_iso,dc_hex,dc_tet,dc_ort,dc_mon,dc_tri,vel, 
		   ced_scc,loc_rot,scca_old_mode); 
  /* we still need to rotate the matrix to get the C style inverse
     (yes, i checked) */
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      scc_irot[j*3+i] = loc_rot[i*3+j];
}
/* 

compute the pole figure density, the summed projection of the cartesian
reference system grains at location xp[r,theta,phi]
time is only passed for reference, so is the deformation matrix F


if rotate is true, will rotate the ODF density function accordin to
xp into a local Cartesian system else, simply use the global Caretsian
system


icall is just an output fielname flag

*/
//#define SUPER_DEBUG
void drex_compute_pdens(struct drex_para *drex,
			COMP_PRECISION *xp,
			COMP_PRECISION time,
			my_boolean rotate,
			my_boolean o1xyz_out,int icall,
			my_boolean debug_out)
{
  char fname[200];
  FILE *out;
  int i,j;
#ifdef DREX_SUPER_DEBUG
  int nxny,os1,os2;
#endif
  int irotate,io1xyz_out,iuse_old_pdens;
  int offset;
  /* 
     
  use the old pdens routines?

  */
  iuse_old_pdens = 1;

  /* length of one whole array  */
  offset = drex->np[0] * drex->np[1] * 3;
  /* logic flags */
  irotate = (rotate)?(1):(0);
  io1xyz_out = (int)o1xyz_out;
  /* 
     call the fortran routine that deals with the computation 
  */
  drex_compute_pdens_ftrn(drex->acs,drex->acs_ens,
			  drex->odf,drex->odf_ens,
			  &drex->size,
			  drex->pdens,(drex->pdens+offset),
			  (drex->np+0),(drex->np+1),
			  xp,&time,&irotate,&io1xyz_out,
			  &icall,&iuse_old_pdens);
  if(debug_out){
    /* output of olivine component, first axis */
    sprintf(fname,"o1pdens.%0g.data",time);
    fprintf(stderr,"drex_compute_pdens: writing olivine [100] at time %g to %s\n",
	    time,fname);
    out=fopen(fname,"w");
    for(i=0;i  < drex->np[1];i++){/* lat loop */
      for(j=0;j < drex->np[0];j++){ /* lon loop */
	fprintf(out,"%.3e\n",drex->pdens[i*drex->np[0]+j]);
      }
    }
    fclose(out);
  }
}
/* 
   
for debugging purposes

*/
my_boolean drex_is_symmetric(COMP_PRECISION *x, int n)
{
  int i,j;
  COMP_PRECISION diff;
  for(i=0;i<n;i++)
    for(j=i;j<n;j++){
      diff = fabs(x[i*n+j]-x[j*n+i]);
      if(diff > 5e-7){
	fprintf(stderr,"i: %i j: %i x(i,j): %g x(j,i): %g dx: %g\n",
		i+1,j+1,x[i*n+j], x[j*n+i],x[i*n+j]-x[j*n+i]);
	return FALSE;
      }
      if(diff > 5e-15)		/* average */
	x[i*n+j] = x[j*n+i] = (x[i*n+j] + x[j*n+i])/2.0;
    }
  return TRUE;
}


void drex_test_if_initialized(struct drex_para *drex,char *program)
{
  if(!drex->initialized){
    fprintf(stderr,"%s: error: DREX structure was not initialized\n",
	    program);
    exit(-1);
  }
}
/* 

compute average axes orientations for olivine and enstatite [100] (a)
using unweighted and odf weighted averages

*/
void drex_compute_avg_axes(struct drex_para *drex,
			   COMP_PRECISION *xp,
			   my_boolean rotate)
{
  int rot,i,j;
  COMP_PRECISION co[6],ce[6];
  rot = (rotate)?(1):(0);
  drex_compute_avg_axes_ftrn(drex->acs,drex->acs_ens,
			     drex->odf,drex->odf_ens,
			     &drex->size,xp,&rot,
			     drex->avg_axes,drex->avg_axes_ens,
			     drex->dev_axes,drex->dev_axes_ens,
			     &drex->olivine_jindex);
  /* convert from fortran to C so that the fast index is the dimension,
     the slow one the unweighted/weighted axis */
  a_equals_b_vector(co,drex->avg_axes,6);
  a_equals_b_vector(ce,drex->avg_axes_ens,6);
  for(i=0;i<2;i++)
    for(j=0;j<3;j++){
      drex->avg_axes[i*3+j]     = co[i+j*2];
      drex->avg_axes_ens[i*3+j] = ce[i+j*2];
    }
}

/* 

get stiffness matrices for a few minerals

sav[6,6] stiffness matrix, returned C style

temp: temperature in deg 
depth: depth in km

modes:

ol and en from Estey
 DREX_CTP_OL_ESTEY 
 DREX_CTP_EN_ESTEY 
ol and en from newer sources
 DREX_CTP_OL_NEW 
 DREX_CTP_EN_NEW 
original KR
 DREX_C_OL_KR 
 DREX_C_EN_KR 
jules
 DREX_C_OL_JULES 
 DREX_C_EN_JULES 

*/
void drex_get_sav_constants(COMP_PRECISION *sav, int mode, 
			    COMP_PRECISION temp,
			    COMP_PRECISION pressure)
{
  switch(mode){
  case DREX_CTP_OL_ESTEY:		/* olivine at T, p */
    drex_get_olivine_sav_ftrn(sav, &temp, &pressure);
    break;
  case DREX_CTP_EN_ESTEY:			/* enstatite at T, p */
    drex_get_enstatite_sav_ftrn(sav, &temp, &pressure);
    break;
  case DREX_CTP_OL_NEW:			/* new olivine at T, p */
    drex_get_olivine_new_sav_ftrn(sav, &temp, &pressure);
    break;
  case DREX_CTP_EN_NEW:			/* new enstatite at T, p */
    drex_get_enstatite_new_sav_ftrn(sav, &temp, &pressure);
    break;
    /* 
       KR reference values
    */
  case DREX_C_OL_KR:			/* olivine reference from orig KR*/
    drex_get_olivine0_sav_ftrn(sav);
    break;
  case DREX_C_EN_KR:			/* enstatite reference from orig KR */
    drex_get_enstatite0_sav_ftrn(sav);
    break;
    /* 
       Browaeys & Chevrot (2004) reference
    */
  case DREX_C_OL_JULES:			/* olivine at 1500 K as in Browaeys */
    drex_get_olivine1_sav_ftrn(sav);
    break;
  case DREX_C_EN_JULES:			/* enstatite at room T as in Browaeys */
    drex_get_enstatite1_sav_ftrn(sav);
    break;
  default:
    fprintf(stderr,"drex_get_sav_constants: error: mode %i undefined\n",
	    mode);
    exit(-1);
    break;
  }
  /* this is unnecessary, as the tensors are symmetric, but we might
     have filled in the other part of the triangle */
  transpose_nbyn_inplace(sav,6); /* flip to C style */
}

/* 
   
compute the tensor norm from a [6,6] elasticity matrix

*/
COMP_PRECISION drex_ctn_36(COMP_PRECISION *sav)
{
  COMP_PRECISION sav21[21];
  drex_fullsym6_ftrn(sav);	/* make sure we have symmetric
				   elements (this way we don't 
				   care about C or F sorting */
  drex_v21d_ftrn(sav,sav21); /* convert [6,6] into [21] */
  return sqrt(vecdotp(sav21,sav21,21));
}

/* 
   
simple VHR average of tensor s1 (frac1) + s2 (1-frac1) 

*/
void drex_voigt_avg(COMP_PRECISION *s1,COMP_PRECISION *s2,COMP_PRECISION frac1,
		    COMP_PRECISION *so)
{
  int i;
  COMP_PRECISION frac2;
  frac2 = 1.0-frac1;
  for(i=0;i<36;i++)
    so[i] = frac1 * s1[i] + frac2 * s2[i];
}

/* 

given a location on the surface of a sphere given in spherical
coordinates r,theta,phi obtain three Euler angles alpha, beta, gamma
as in Dahlen and Tromp, p. 920 such that the corresponding rotation
matrix will rotate a globally Cartesian system vector or tensor into a
local Cartesian system where 

REG_CART_CONVENTION: x = South, y = East, z = up
SW_CART_CONVENTION:  x = North, y = east, z = down
REG2SW_CONVENTION: go from regular to SW convention

cart_convention is a pointer for compatibility with fortran

output is in radians

*/
void drex_calc_euler_rad_from_xp(COMP_PRECISION *xp,COMP_PRECISION *alpha,
				 COMP_PRECISION *beta, COMP_PRECISION *gamma,
				 int *cart_convention)
{
  switch(*cart_convention){
  case DREX_NO_ROTATION:	/* no rotation */
    *alpha = 0.0;
    *beta  = 0.0;
    *gamma = 0.0;
    break;
  case DREX_REG2SW_CONVENTION:	/* go from regular convention to
				   surface wave convention */
    *alpha = 0.0;
    *beta  = DREX_PI_C;
    *gamma = 0.0;
    break;
  case DREX_REG_CART_CONVENTION:	/* regular: S, E, U */
    *alpha = xp[DREX_PHI];
    *beta  = xp[DREX_THETA];
    *gamma = 0.0;
    break;
  case DREX_SW_CART_CONVENTION:	/* for surface wave stuff: N, E, D */
    *alpha = xp[DREX_PHI];
    *beta  = xp[DREX_THETA]+ DREX_PI_C;
    *gamma = 0.0;
    break;
  default:
    fprintf(stderr,"drex_calc_euler_rad_from_xp: error: Cartesian system convention %i undefined\n",
	    *cart_convention);
    break;
  }
}
void drex_calc_euler_rad_from_xp_(COMP_PRECISION *xp,COMP_PRECISION *alpha,
					   COMP_PRECISION *beta, COMP_PRECISION *gamma,
					   int *cart_convention)
{
  drex_calc_euler_rad_from_xp(xp,alpha,beta, gamma,cart_convention);
}


void drex_assign_tau_from_system(COMP_PRECISION *tau, /* output slip system */
				 int type, /* type switch */
				 COMP_PRECISION *tau_rss, /* input slip ssytem */
				 COMP_PRECISION pressure, /* for adjust mode GPa */
				 COMP_PRECISION ptrans /* transition pressure */
				 )
{
  static COMP_PRECISION 
    tau_def_a[4]=DREX_TAU_OLIVINE_DEF_ATYPE,
    tau_def_b[4]=DREX_TAU_OLIVINE_DEF_BTYPE,
    tau_def_c[4]=DREX_TAU_OLIVINE_DEF_CTYPE,
    tau_def_d[4]=DREX_TAU_OLIVINE_DEF_DTYPE,
    tau_def_e[4]=DREX_TAU_OLIVINE_DEF_ETYPE,
    tau_def_hp[4]=DREX_TAU_OLIVINE_DEF_HIGHP;

  switch(type){
  case DREX_FABRIC_FREE:
    a_equals_b_vector(tau,tau_rss,4);	/* free format*/
    break;
  case DREX_FABRIC_ATYPE:
    a_equals_b_vector(tau,tau_def_a,4);	/* regular, A type slip system */
    break;
  case DREX_FABRIC_BTYPE:
    a_equals_b_vector(tau,tau_def_b,4);/* type B  */
    break;
  case DREX_FABRIC_CTYPE:
    a_equals_b_vector(tau,tau_def_c,4);/* type C */
    break;
  case DREX_FABRIC_DTYPE:
    a_equals_b_vector(tau,tau_def_d,4);/* type D  */
    break;
  case DREX_FABRIC_ETYPE:
    a_equals_b_vector(tau,tau_def_e,4);/* type E  */
    break;
  case DREX_FABRIC_HIGHP:
    a_equals_b_vector(tau,tau_def_hp,4);/* high pressure  */
    break;
  case DREX_FABRIC_SWITCH_A_HIGHP:
    if(pressure < ptrans)
      a_equals_b_vector(tau,tau_def_a,4);	/* regular, A type slip system */
    else
      a_equals_b_vector(tau,tau_def_hp,4);/* high pressure  */
    break;
  default:
    fprintf(stderr,"drex_assign_tau_from_system: slip system %i undefined\n",
	    type);
    exit(-1);
    break;
  }
}
