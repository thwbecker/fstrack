#include "fstrack.h"
//
//
// $Id: ellsphere.c,v 1.1 2004/02/13 04:07:11 becker Exp becker $
//
//
//     runge kutta integration driver, calls odeintell
//
//     this routine assumes that velocities are given as 
//     phi = 0. ... (nphi-1)*dphi
//     theta = dtheta/2 ... pi-dtheta/2
//
//
void ellsphere(COMP_PRECISION *work,COMP_PRECISION ti,COMP_PRECISION tf,
	       COMP_PRECISION dt,COMP_PRECISION eps,
	       BOOLEAN calc_strain,BOOLEAN calc_lya,BOOLEAN remove_strain,
	       struct der_par *dp)
{
//_______________________________________________________________________
//
//     ti  -> initial time 
//     tf  -> final time   both of these have to be in the right units
//                         for the time interpolation scheme to work
//     dt -> suggested first time step
//     work -> array with r,theta,phi coordinates and F components, if 
//     ifellipse is set, ie. work is work(3) or work(12)
//     ifremovestrain: remove the symmtetric part of the velocity gradient
//                    matrix depending on the rules laid out in ellderiv
//
//_______________________________________________________________________

  int nwork,i,nok,nbad;
  static COMP_PRECISION smaxstep = 1.0;
  COMP_PRECISION stepmin,alpha[3];
  //
  //set parameters and call odeint
  //
  if (!calc_strain){
    //
    // no strain tracing and vel gradient matrix needed, just velocity advection
    //
    nwork = 3;
    if(calc_lya){
      fprintf(stderr,"ellsphere: You cannot calculate lyapunov exponents\n");
      fprintf(stderr,"ellsphere: without calculating strain (ifellipse)!\n");
      fprintf(stderr,"ellsphere: ifellipse,iflyap: %i, %i\n",
	      calc_strain,calc_lya);
      exit(-1);
    }
  }else{			/* need to compute strain */
    nwork = 12;
  }
  //
  //     max step is internally limited to smaxstep
  //
  if(fabs(dt) > smaxstep){
    dt = (fabs(dt)/dt) * smaxstep;
  }
  /* there is no minimum step */
  stepmin=0.0;
  //
  //     integrate to obtain position and deformation matrix at t2
  //
  nr_odeint(work,nwork,ti,tf,eps,dt,stepmin,&nok,&nbad,dp);
  if (calc_lya){
    /* 
       Lyapunov part 
    */
    // normalize
    gramschmidt(work,alpha);
    // and increment lyapunov exponents
    for(i=0;i<3;i++)
      dp->rlyap[i] += log(alpha[i]);
  }
  /* tracer should be within radial bounds since this is 
     corrected for in  rkqsell */
  if(work[0] > 1.0){
    fprintf(stderr,"ellsphere: tracer x_r > 1: %g\n",work[0]);
    exit(-1);
  }
  //     adjust the location vector, work(1...3)
  check_phys_lim_tracer(work);
}

