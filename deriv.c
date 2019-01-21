#include "fstrack_flow.h"
/*

routines dealing with Runge Kutta derivative parameters 
and forming the derivative of all integrated properties

also holds the physical limits routine and error checkers

$Id: deriv.c,v 1.20 2010/12/31 18:59:17 becker Exp $

*/


/* 


DERIVATIVE STUFF


*/

/* 

   driver for all the derivative routines 

*/
void rk_derivs(COMP_PRECISION x, COMP_PRECISION *y,
	       COMP_PRECISION *dydx, int n, 
	       struct der_par *dp)
{
  /* use the Cartesian system internally to keep track of strain */
  static my_boolean vgm_cart = TRUE;
  my_boolean frozen;
  COMP_PRECISION xp[3];
#ifdef SUPER_FSTRACK_DEBUG
  lonlatz_from_xp(y,xp);
  fprintf(stderr,"rk_derivs: n: %6i nok: %6i nbad: %6i t: %11g llz: %.3e,%.3e,%.3e\n",
	  n,dp->nok,dp->nbad,x,xp[FSTRACK_X],xp[FSTRACK_Y],xp[FSTRACK_Z]);
#endif
  if(n > 12){
    /* 
       compute the derivatives with respect to time for dydx(0,1,2):
       coordinates (NOT velocities, those are only related, see
       fse_deriv) x(3...11): strain, ie. VGM . strain. this routine also
       computes the local velocities, dp->v_loc
       
    */
    fse_derivs_wrapper(x,y,dydx,12,dp,dp->v_loc,dp->vgm,vgm_cart,TRUE,&frozen,
		       dp->strain_fraction_from_gamma);   
    /* do we need to change the slip system? */
    if(dp->drex_type == DREX_FABRIC_SWITCH_A_HIGHP){
      xp[FSTRACK_Z] = ZDEPTH(y[FSTRACK_R]);
      drex_assign_tau_from_system(dp->drex->tau,dp->drex_type,
				  dp->drex_rss_tau,
				  pressure_model(xp[FSTRACK_Z]),
				  dp->drex_sstp);
      //fprintf(stderr,"z: %11g p: %11g pc: %11g --> {%11g, %11g, %11g, %11g}\n",xp[FSTRACK_Z],pressure_model(xp[FSTRACK_Z]), dp->drex_sstp,dp->drex->tau[0],dp->drex->tau[1],dp->drex->tau[2],dp->drex->tau[3]);
    }
    /* 
       compute the DREX derivatives from the velocity gradient matrix
       by calling the driver for drex_deriv_ftrn
    */
    drex_deriv(dp->vgm,y,dydx,dp->drex);
  }else{
    /* just velocities (and strain) */
    fse_derivs_wrapper(x,y,dydx,n,dp,dp->v_loc,dp->vgm,vgm_cart,TRUE,&frozen,
		       dp->strain_fraction_from_gamma);   
  }
}
/* 

wrapper for the FORTRAN based routine which calculates coordinate time
derivatives (NOT velocities), the velocities, and the velocity
gradient matrix. 

those are all normally in spherical coordinates, unless
strainmat_cartesian is set. in this case, the strain derivatives
(dydx+3), VGM, and the velocities will be in the Cartesian system


allow frozen strain will allow the remove_symm_strain option, else off

*/
void fse_derivs_wrapper(COMP_PRECISION time, COMP_PRECISION *y,
			COMP_PRECISION *dydx, int n,
			struct der_par *dp,
			COMP_PRECISION *vp,
			COMP_PRECISION *vgm,
			my_boolean strainmat_cartesian,
			my_boolean allow_frozen_strain,
			my_boolean *frozen,
			my_boolean alpha_strain_fraction)
{
  int irs;			/* remove strain flag */
  int ism_cart=1;			/* cartesian strain flag */
  int ivgm_const=0;		/* constant VGM flag */
  int ifrozen;			/* for the query of fse */
  int ialpha_strain;	/* only use fraction of strain_rate matrix? */
  static my_boolean warned[3]={FALSE,FALSE,FALSE};
  COMP_PRECISION alpha;
  if(allow_frozen_strain){	/* allow the option of frozen strain, but 
				   don't switch on yet */
    if(dp->remove_strain){
      irs = 1;			/* indeed remove symm strain */
      if(!warned[0]){
	fprintf(stderr,"\nWARNING\n\nfse_derivs_wrapper: removing strains in fse_deriv according to freezing rule\n\n");
	warned[0] = TRUE;
      }
    }else{
      irs = 0;
    }
  }else{
    irs = 0;
  }
  ism_cart = (strainmat_cartesian)?(1):(0);
  if(alpha_strain_fraction){
    if(!warned[1]){
      fprintf(stderr,"\nfse_derivs_wrapper: WARNING: using alpha fraction of the strain\n");
      warned[1] = TRUE;
    }
    ialpha_strain = 1;
#ifdef FSTRACK_USE_GGRD    
    /* 

    determine alpha from the r,theta,phi location

    */
    if(!ggrd_grdtrack_interpolate_rtp(y[0],y[1],y[2],dp->ggrd_alpha,&alpha,FALSE,FALSE,PREM_RE_KM)){
      PE("fse_derivs_wrapper: interpolation error with alpha grids");
      exit(-1);
    }
    /* adjust alpha */
    if(alpha < 0)
      alpha = 0.0;
    if(alpha > 1)
      alpha = 1.0;
    //fprintf(stderr,"lon: %11g lat: %11g z: %11g\talpha: %11g\n",PHI2LONGITUDE(y[2]),THETA2LATITUDE(y[1]),ZDEPTH(y[0]),alpha);
#else
    fprintf(stderr,"fse_derivs_wrapper: need to compile with FSTRACK_USE_GGRD for alpha\n");
    exit(-1);
#endif
  }else{
    ialpha_strain = 0;
  }
  if(dp->constant_vgm_for_debugging){
    /* 
       assign a constant velocity gradient matrix for testing purposes
    */
    if(!warned[2]){
      fprintf(stderr,"\nfse_derivs_wrapper: WARNING: using constant VGM (S: %g T: %g W: %g)\n\n",
	      dp->cvfd_stw[0],dp->cvfd_stw[1],dp->cvfd_stw[2]);
      
      warned[2] = TRUE;
    }
    ivgm_const = 1;
#ifndef FSTRACK_DEBUG
    fprintf(stderr,"fse_derivs_wrapper: error: only debugging mode can deal with constant VGM\n");
    exit(-1);
#endif
  }
  fse_derivs(&time,y,dydx,&n,dp->vr,dp->vt,dp->vp,&dp->n[FSTRACK_R],
	     &dp->n[FSTRACK_THETA],&dp->n[FSTRACK_PHI],&dp->dtheta,&dp->dphi,dp->rlevels,
	     dp->vtimes,&dp->nvtimes,&irs,dp->vhdivmax,
	     vp,vgm,&ism_cart,&ivgm_const,dp->cvfd_stw,&ifrozen,
	     &ialpha_strain,&alpha);
  if(ifrozen)
    *frozen = TRUE;
  else
    *frozen = FALSE;
}

/* 

   INITIALIZATION OF DERIVATIVE PARAMETER STRUCTURE

   rot_g_coord_sys: rotate the grain coordinate system to E-N-U
   or not
*/
void init_der_par(struct der_par **dp)
{
  *dp = (struct der_par *)calloc(1,sizeof(struct der_par));
  if(!(*dp))
    MEMERROR("init_der_par");
  (*dp)->rkeps = RK_EPS_DEF; // Runge Kutta precision
  (*dp)->remove_strain = FALSE;
  (*dp)->strain_fraction_from_gamma = FALSE;
  /* 
     other computational modes
   */
  (*dp)->calc_lyapunov = FALSE;
  (*dp)->history = FALSE;
  (*dp)->hmax = HMAX_DEF;
  /* 
     
  different ways to compute the error scale
  
  0: default
  1: make all F components vary on the same scale
  
  */
  (*dp)->error_mode = 0;
  /* 
     go from cm/yr to R/timescale[yr], 
     this is a velocity in cm/yr
  */
  (*dp)->velscale = (R_E/TIMESCALE)*1e5;
  //
  // lyapunov related
  //
  (*dp)->renorm = LYA_RENORM_DEF;// renorm if L(i) > renorm
  (*dp)->rlbailout = LYA_RLBAILOUT_DEF;
  (*dp)->tmax = LYA_TMAX_DEF;/* this needn't be a realistic 
				time necessarily */
  /* texture stuff */

   /* averaging: 0 hill 0.5 VHR 1: reuss*/
  (*dp)->vhr_favg = DREX_FAVG;



  (*dp)->texture_mode = NO_TEXTURE;		/* texture mode */
  (*dp)->drex_epsfac = DREX_YERRFAC_DEF; /* multiple of RK eps to use */
  (*dp)->drex_Mob = DREX_MOB_DEF; /* GBM parameter */
  (*dp)->drex_lambda = DREX_LAMBDA_DEF; /* GBM parameter */

  (*dp)->drex_sstp = DREX_SS_PSWITCH; /* transition pressure GPa */
  (*dp)->drex_Xol = DREX_XOL_DEF; /* olivine fraction */
  (*dp)->drex_chi = DREX_CHI_DEF; /* chi parameter */
  (*dp)->drex_size3 = DREX_SIZE3_DEF;
  (*dp)->drex_type = DREX_FABRIC_ATYPE;	/* type A deformation is default */
  (*dp)->drex_save_pole = FALSE;	/* don't save the ODFs by default */
  (*dp)->drex_module = EMOD_PT; /* by default, use the old varying
				   moduli p,T derivatives*/
  (*dp)->constant_vgm_for_debugging = FALSE; /* this should be FALSE,
						normally */
  (*dp)->rotate_grain_coord_sys = TRUE;	/* rotate grains to E-N-U
							   system */
  (*dp)->drex_start_grains_oriented = FALSE; /* start with random */
  (*dp)->lpo_tracer_count = -1;	/* count the number of times we
				   compute LPOs, this way the first
				   one will start at zero */
  (*dp)->print_o1xyz = 0;	/* don't print olivine [100] (for =1 ,
				   or all axes for = 2) at each
				   timestep
				*/
  /* 
     this has to be initialized here so that 
     the tracer routines know how much storage to allocate
     for poles, if needed
  */
  (*dp)->drex_np[0] = DREX_POLE_NX;
  (*dp)->drex_np[1] = DREX_POLE_NY; 
  (*dp)->drex_npole = (*dp)->drex_np[0] * 
    (*dp)->drex_np[1] * 3 * 2;
}
/* 




ERROR CHECKING STUFF


*/


/* 

routine to check the physical limits of all quantities
this gets called after each RK sub-step.

the called_from_rkck flag is true when called from within rkck,
else it should be FALSE. this allows one to distinguish different operating 
procedures

*/
void rk_check_phys_limit(COMP_PRECISION *x, int n, 
			 struct der_par *dp,
			 my_boolean called_from_rkck)
{
  if(!called_from_rkck){
    /* 
       check the tracer locations and adjust if needed
    */
    check_phys_lim_tracer(x,x);
  }
  /* 
     texture part
  */
  if(n > 12){
    drex_check_phys_limits(x,n,dp->drex);
  }
}

/*     
     check the physical bounds of a spherical system point
     in r,theta,phi space
     
     returns the corrected location

*/     
/* fortran wrapper */
void check_phys_lim_tracer__(COMP_PRECISION *xloc,
			     COMP_PRECISION *xloc_new)
{check_phys_lim_tracer(xloc,xloc_new);}
void check_phys_lim_tracer_(COMP_PRECISION *xloc,
			    COMP_PRECISION *xloc_new)
{check_phys_lim_tracer(xloc,xloc_new);}

void check_phys_lim_tracer(COMP_PRECISION *xloc,
			   COMP_PRECISION *xloc_new)
{
  //     a little smaller than actual CMB radius
  static COMP_PRECISION r_cmb = 0.5462;
  /* 
     work on new vector copy, don't modify xloc by defaults 
  */
  a_equals_b_vector3d(xloc_new, xloc);
  /*

  adjust position values, first radially

  */
  if(xloc_new[FSTRACK_R] > 1.0){
    
    fprintf(stderr,"check_phys_lim_tracer: r surface adjust: %g\n",
	    xloc_new[FSTRACK_R]);
    xloc_new[FSTRACK_R]=1.0;
  }
  if(xloc_new[FSTRACK_R] < r_cmb){
    fprintf(stderr,"check_phys_lim_tracer: r CMB adjust: %g\n",
	    xloc_new[FSTRACK_R]);
    
    xloc_new[FSTRACK_R] = r_cmb;
  }
  //
  // check the theta range
  //
  // 0 <= theta <= pi
  //
  if (xloc_new[FSTRACK_THETA] < 0.0){
    xloc_new[FSTRACK_THETA] = -xloc_new[FSTRACK_THETA];
    xloc_new[FSTRACK_PHI] += PI;
  }
  if (xloc_new[FSTRACK_THETA] > PI){
    xloc_new[FSTRACK_THETA] = TWOPI - xloc_new[FSTRACK_THETA];
    xloc_new[FSTRACK_PHI] += PI;
  }
  //
  // check the phi range
  //
  // 0<= phi <= 2*pi
  if (xloc_new[FSTRACK_PHI] > TWOPI){
    xloc_new[FSTRACK_PHI] -= TWOPI;
  }
  if (xloc_new[FSTRACK_PHI] < 0.0){
    xloc_new[FSTRACK_PHI] += TWOPI;
  }
}
//     
//     find a physics a priori knowledge based error condition
//     in our case, we want to avoid tracers exiting the surface of the 
//     Earth
//     
//     the results should be returned in units of eps precision
//     
//     
void check_physics_based_error(COMP_PRECISION *y,
			       COMP_PRECISION *dydx,
			       int *n,
			       COMP_PRECISION *eps,
			       COMP_PRECISION *errmax,
			       COMP_PRECISION *h)
{
  COMP_PRECISION tmp_err;
  //     
  //     check for radial direction coordinate
  //
  tmp_err = (y[FSTRACK_R] - 1.0)/ (*eps);
  if(tmp_err > 1.0){
    *errmax = MAX(tmp_err,*errmax);
#ifdef FSTRACK_DEBUG
    fprintf(stderr,"check_phy_based error: dydx/h/r: %g %g %g\n",
	    dydx[FSTRACK_R],*h,y[FSTRACK_R]);
#endif
  }
}
/* 
   
given n dimensional vectors for Runge Kutta, determine
up to which n the error checking should proceed

 */
int determine_n_error(int n,struct der_par *dp)
{
  if(dp->strain_rate_control)
    return MIN(12,n);		/* check only up to finite
				   strain, tops, the texture variables are 
				   controled by 1/strainrate max time stepping
				*/
  else				/* check all variables */
    return n;
}
/* 
   determine the scale for each component of y[nvar] which
   will be used to determine the relative error
*/
void determine_error_scale(COMP_PRECISION time,
			   COMP_PRECISION *y,
			   COMP_PRECISION *dydx,
			   COMP_PRECISION h,
			   int nvar,struct der_par *dp,
			   COMP_PRECISION *yscal)
{
  int ilim,i;
  COMP_PRECISION ys_f_min,ytmp;
  static my_boolean init = FALSE;
  static int n_error;
  if(!init){
    /* determine up to which variable dimension the error should be determined */
    n_error = determine_n_error(nvar,dp);
    if(dp->drex_epsfac < 0){
      fprintf(stderr,"determine_error_scale: error: drex_epsfac < 0\n");
      exit(-1);
    }
    init = TRUE;
  }
  /* 
     find the error scaling using different approaches
  */
  switch(dp->error_mode){	/* 
				   switch through 
				   different ways of determining
				   the error 
				*/
  case 0:
    /* 
       
    all error scales are relative

    */
    ilim = MIN(12,nvar);	/* limit to locations or strain */
    for (i=0;i < ilim;i++){
      /* location and deformation */
      yscal[i] = fabs(y[i]) + fabs(dydx[i]*h) + dp->rkeps; 
    }
    /* 
       if drex_epsfac is smaller than DREX_RK_ERR_OFF_FACTOR,
       use also the error bounds for texture, else only for flow
    */
    for (i=12;i < n_error;i++){
      /* texture  */
      yscal[i]  = fabs(y[i]) + fabs(dydx[i]*h) + dp->rkeps; 
      yscal[i] *= dp->drex_epsfac;
    }
    break;
  case 1:
    /* 
       
    error scale relative but strain has all components 
    normalized the same

    */
    for (i=0;i < 3;i++){	
      /* first, for coordinates */
      yscal[i] = fabs(y[i]) + fabs(dydx[i]*h) + dp->rkeps; 
    }
    /* 
       now for deformation, find the minimum and have all scale
       with this yscale
    */
    ys_f_min = 1e20;
    for (i=3;i < 12;i++){
      ytmp = fabs(y[i]) + fabs(dydx[i]*h) + dp->rkeps;
      ys_f_min = MIN(ys_f_min,ytmp);
    }
    for (i=3;i < 12;i++)
      yscal[i] = ys_f_min;
    /* now texture */
    for(i=12;i < n_error;i++){
      yscal[i] = fabs(y[i]) + fabs(dydx[i]*h) + dp->rkeps; 
      yscal[i] *= dp->drex_epsfac;
    }
    break;
  case 2:
    /* 
       all scales fixed to constant
    */
    ilim = MIN(12,nvar);
    for (i=0;i < ilim;i++)	/* location and deformation */
      yscal[i] = dp->rkeps; 
    ytmp =  dp->rkeps * dp->drex_epsfac;
    for(i=12;i < n_error;i++) {
      /* texture */
      yscal[i] = ytmp;
    }
    break;
  default:
    fprintf(stderr,"determine_error_scale: error: error mode %i undefined\n",
	    dp->error_mode);
    exit(-1);
    break;
  }
}

/* 

given  velocity fields in the internal format stored in dp, compute
velocities, vp, and the strain rate tensor, e, in intrinsic units at time `time'
and location x[12], given in r, theta, phi system

dimensions: x[12],vp[3],e[9]

if cart_system is set, will return v[3] and e[9] in Cartesian system,
else both are spherical system


WARNING: x has to be dimensioned [12]
*/

void calc_vel_and_strain_rate(COMP_PRECISION time, COMP_PRECISION *x,
			      struct der_par *dp,COMP_PRECISION *v,
			      COMP_PRECISION *e,my_boolean cart_system)
{
  COMP_PRECISION dx[12];
  my_boolean frozen;
  //
  // obtain velocities v and 
  // velocity gradient matrix at x[3] and time `time'
  //
  if(dp->strain_fraction_from_gamma)
    fprintf(stderr,"calc_vel_and_strain_rate: WARNING: fse_deriv using only alpha strain fraction\n");
  fse_derivs_wrapper(time,x,dx,12,dp,v,dp->vgm,cart_system,TRUE,&frozen,
		     dp->strain_fraction_from_gamma);
  if(frozen){
    fprintf(stderr,"calc_vel_and_strain_rate: WARNING: fse_derivs returned frozen\n");
  }
  //
  //
  //
  // calculate strain-rate matrix from velocity gradient
  calc_cd_symm_part(dp->vgm,e);
}
/* 
   compute the maximum timestep from scaling factor times the inverse
   of the characteristic strain rate. need to input time and position
   x[3], will also return the eigenvalues of the strain rate matrix
*/
COMP_PRECISION calc_max_dt_from_strain_rate(COMP_PRECISION time, COMP_PRECISION *x,
					    struct der_par *dp)
{
  COMP_PRECISION dx[12],e[9],v[3],eval[3];
  my_boolean frozen;
  // obtain velocities v and velocity gradient matrix at x[3] and time `time', don't allow 
  // frozen strain
  fse_derivs_wrapper(time,x,dx,12,dp,v,dp->vgm,FALSE,FALSE,&frozen,
		     dp->strain_fraction_from_gamma);
  // calculate strain-rate matrix from velocity gradient
  calc_cd_symm_part(dp->vgm,e);
  return dp->eps_strain / char_strain_rate_from_e(e,eval);
}

/* 

given velocity fields in the internal format stored in dp, compute
velocities, vp, and the strain rate tensor, e, in intrinsic units at
time `time' and location x[0,1,2 ...], given in r, theta, phi system.

WARNING: x has to be dimensioned [12]

returns velocities in input units (cm/yr) and strainrates in 1/s timescale


if cart_system is TRUE, will use cartesian system for e and v


*/
void calc_vel_and_strain_rate_units(COMP_PRECISION time, 
				    COMP_PRECISION *x,
				    struct der_par *dp,
				    COMP_PRECISION *v,
				    COMP_PRECISION *e,
				    my_boolean cart_system)
{
  calc_vel_and_strain_rate(time,x,dp,v,e,cart_system);
  scale_vector3d(v,dp->velscale);
  /* [e']=1/t' , t'=t/t_c */
  scale_vector(e,1.0/TIMESCALE_S,9);
}

