/*


calculate the infinite strain axis in the global Cartesian system at
any given point following Kaminski and Ribe (G^3, 2002GC000222, 2002)

uses the tracer's astate state. 

$Id: calc_isa.c,v 1.17 2006/04/18 17:44:51 twb Exp twb $

*/
#include "fstrack_flow.h"
void calc_isa(struct trc *tracer,int astate,struct mod *model)
{
  COMP_PRECISION omega[3],oamp;
  calc_pipar(tracer->state[astate].x,tracer->state[astate].t,
	     tracer->state[astate].isa,&tracer->state[astate].pipar,omega,
	     model->dp);

  if(norm3d(tracer->state[astate].isa) > EPS_COMP_PREC){
    tracer->state[astate].isa_exists = TRUE;
    //oamp=norm(omega,3);fprintf(stderr,"1 %11g %11g %11g: %11g\n",omega[0]/oamp,omega[1]/oamp,omega[2]/oamp,oamp);
  }else{			/* 
				   i checked the omega values
				   here. for non-existing ISA axes,
				   there were very large rotational
				   components. for some \Omega > 1
				   values, ISA axes could be found,
				   though not for the largest ones
				*/
    tracer->state[astate].isa_exists = FALSE;
    //oamp=norm(omega,3);fprintf(stderr,"0 %11g %11g %11g: %11g\n",omega[0]/oamp,omega[1]/oamp,omega[2]/oamp,oamp);
  }
}


/* 

calculate the ISA axis, grain orientation lag (PI parameter) using the
drex routines, and modelled after the drex example program

input is a position xp[r,theta,phi]

output:

       isac[3]: ISA axis at xp in the global Caretsian system
       gol:     Grain Orientation lag, the Pi parameter
       omega:   rotation vector

this routine is modeled after the pipar routine in Drex



*/

void calc_pipar(COMP_PRECISION *xp, COMP_PRECISION time,
		COMP_PRECISION *isac, COMP_PRECISION *golc,
		COMP_PRECISION *omega,
		struct der_par *dp)
{
  COMP_PRECISION xstate[12],xstate_orig[12],dxstate[12],dt;
  COMP_PRECISION vpc[3],vpb[3],vpf[3];
  COMP_PRECISION vgmc[9],vgmb[9],vgmf[9],isab[3],
    isaf[3],char_strainc,tmp;
  my_boolean frozen;
  COMP_PRECISION theta_isa_f,theta_isa_b,dtheta;
  static int nwork = 12;
  *golc = 10.0;			/* max GOL is 10 */
  /* 
     assign location 
  */
  a_equals_b_vector3d(xstate,xp);	
  unity_matrix((xstate+3),3);	/* unnecessary, really */
  /* save center point state */
  a_equals_b_vector(xstate_orig,xstate,12);
  /* 
     compute velocity and velocity gradient matrix at center point
     xp. velocities and vgm are in Cartesian
  */
  if(dp->strain_fraction_from_gamma)
    fprintf(stderr,"calc_pipar: WARNING: fse_deriv using only alpha strain fraction\n");
  fse_derivs_wrapper(time,xstate_orig,dxstate,nwork,dp,vpc,vgmc,TRUE,TRUE,
		     &frozen,dp->strain_fraction_from_gamma);
  if(frozen)
    fprintf(stderr,"calc_pipar: WARNING: fse_deriv returned frozen\n");
  /* compute dt as a small fraction of 1/|v| */
  dt = 0.05 / (norm3d(vpc) * dp->velscale);
  /* 
     compute the ISA at the center from the normalized vgm:
     
     normalize VGM by strain-rate
  */
  char_strainc = char_strain_rate_from_vg(vgmc);
  scale_vector(vgmc,1.0/char_strainc,9);
  /* 
     compute the ISA axis and GOL parameter at the center from the
     normalized velocity gradient matrix. this routine uses a C style
     velocity gradient matrix
  */
  drex_isacalc_from_norm_vgmc(isac,golc,vgmc,omega);
  /* 
     
  march back by dt 
  
  */
  rk_wrapper(xstate,time,time-dt,0,TRUE,FALSE,dp->rkeps,dp,FALSE);
  /* compute velocity and gradient there */
  fse_derivs_wrapper(time-dt,xstate,dxstate,nwork,dp,vpb,vgmb,TRUE,TRUE,
		     &frozen,dp->strain_fraction_from_gamma);
  /* normalize velocity */
  normalize3d(vpb);
  /* normalize the velocity gradient matrix */
  scale_vector(vgmb,1.0/char_strain_rate_from_vg(vgmb),9);
  /* 
     compute the ISA axis and GOL parameter from the normalized 
     velocity gradient matrix
  */
  drex_isacalc_from_norm_vgmc(isab,golc,vgmb,omega);
  if(fabs(vec_sum(isab,3) - 3.0) < EPS_COMP_PREC){
    /* 
       special case, the velocity is the ISA, assign velocity to ISA
       axis
    */
    a_equals_b_vector3d(isab,vpb);
  }
  /* angle between ISA and flow direction */
  theta_isa_b = acos(vec3ddotp(vpb,isab));
  /* 

     march forward by dt 

  */
  a_equals_b_vector(xstate,xstate_orig,12);
  rk_wrapper(xstate,time,time+dt,0,TRUE,FALSE,dp->rkeps,dp,FALSE);
  /* compute velocity and gradient there */
  fse_derivs_wrapper(time+dt,xstate,dxstate,nwork,dp,vpf,vgmf,TRUE,TRUE,&frozen,
		     dp->strain_fraction_from_gamma);
  /* normalize velocity */
  normalize3d(vpf);
  /* normalize the velocity gradient matrix */
  scale_vector(vgmf,1.0/char_strain_rate_from_vg(vgmf),9);
  /* 

     compute the ISA axis and GOL parameter 

  */
  drex_isacalc_from_norm_vgmc(isaf,golc,vgmf,omega);
  if(fabs(vec_sum(isaf,3) - 3.0) < EPS_COMP_PREC){
    /* special case, assign velocity to ISA axis */
    a_equals_b_vector3d(isaf,vpf);
  }
  theta_isa_f = acos(vec3ddotp(vpf,isaf));
  /* 
     Pi parameter, compute from central difference 
     change between -dt and +dr isa axes
  */

  dtheta = fabs(theta_isa_b - theta_isa_f);
  if(dtheta > PI)
    dtheta -= PI;
  /* 

  now assign to the GOL parameter for the central location
  
  \Pi = 1/\eps | D\theta/dt|

  */
  tmp = (dtheta/(2.0*dt))/char_strainc;
  *golc = MIN(*golc, tmp);

  //fprintf(stderr,"theta_1: %11g theta_2: %11g dtheta: %11g eps: %11g dt: %11g Pi: %11g\n",
  //theta_isa_b,theta_isa_f,dtheta,char_strainc,dt,*golc);

}



