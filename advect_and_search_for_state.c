#include "fstrack_flow.h"
/*

  follow tracer to 
  - some backward time, 
  - depth sdepth, or 
  - max strain equals mstrain or depth larger than sdepth
  (see bailoutcrit.c)

  and return the final state

  final state will have the location at positive time t as if the
  velocity were inverted

  return total number of steps if succesful and negative error 
  code if failure

  this is all pretty terrible but works ok so far at some point, maybe
  come up with some nicer scheme to find the critical strain and stuff


  $Id: advect_and_search_for_state.c,v 1.7 2010/12/31 18:59:17 becker Exp $

*/
#ifdef FSTRACK_DEBUG
#define FSTRACK_DEBUG_ITER 
#endif

#define MAXITER1 8000
#define MAXITER2 50
#define MAXITER3 100
#define MAXITER4 50
#define INIT_OVERSHOOT_PREC 0.1// initial relative overshoot tolerance
#define FIN_OVERSHOOT_PREC 0.01/* 
				  final relative overshoot tolerance
				  this will depend on the R-K
				  precision, say if strains are to be
				  correct up to 5e-3, then velocities
				  should be correct up to 2.5e-5
			       */
#define DT_INITS {0.0, 30.0, 1, 0.0, 0.1}/* initial guesses of DOUBLE the 
					    negative timestepping for all
					    modes of bailout types:
					    first:  dummy (forward time, no bailout)
					    second: depth 
					    third:  strain
					    fourth: dummy (backward time, exact)
					    fifth: forward strain
					 */
#define NR_BAILOUT_CRITERIA 5
int advect_and_search_for_state(struct mod *model,int mode,
				struct stt *initial_state,
				struct stt *final_state,
				my_boolean verbose,
				my_boolean itrack_strain,
				COMP_PRECISION *last_dt)
{
  static COMP_PRECISION dt0[NR_BAILOUT_CRITERIA]=DT_INITS;
  COMP_PRECISION dt,overshoot=0.0,dtnew,mstrain,tlast=0.0,final_overshoot_prec;
  int i1,i2,iter1,iter2,tstep;
  static my_boolean bailout[3];// bailout flags, determined in bailout_crit_fulfilled
  i1=i2=iter1=iter2=0;
#ifdef FSTRACK_DEBUG
  /*
    
    do some sanity checks, probably not needed for optimized code
    
  */
  switch(mode){
  case FORWARD_TIME:
  case BACKWARD_DEPTH:
  case BACKWARD_STRAIN:
  case BACKWARD_TIME:
  case FORWARD_STRAIN:
    if(bailout_crit_fulfilled(initial_state,mode,model,&overshoot,bailout,&mstrain)){
      fprintf(stderr,"advect_and_search_for_state: bailout already at init: mode: %i z: %g(%g) e: -(%g) t: %g(%g)\n",
	      mode,ZDEPTH(initial_state->x[FSTRACK_R]),ZDEPTH(model->sdepth),
	      model->maxstrain,initial_state->t,model->tf);
      exit(-1);
    }
    if(((mode == BACKWARD_STRAIN)||(mode == FORWARD_STRAIN))&&(!itrack_strain)){
      fprintf(stderr,"advect_and_search_for_state: error in search_for_state: itrack should be TRUE for mode %i\n",
	      mode);
      exit(-1);
    }
    break;
  default:
    fprintf(stderr,"advect_and_search_for_state: mode %i undefined\n",mode);
    exit(-1);
    break;
  }
#endif
  /*
    
    determine initial (double) trial timesteps
    
  */
  if(mode == BACKWARD_TIME){// go backward in time, thus dt can be exact
    dt = -model->tf*2.0;    // choose double since later halfed
  }else{
    if(mode != FORWARD_STRAIN){
      /* backward strain, depth, or time */
      dt = -dt0[mode];
    }else{
      /* forward strain */
      dt = dt0[mode];
    }
  }
  if(mode == BACKWARD_DEPTH)	/* higher precision for backward depth
				   computation */
    final_overshoot_prec = FIN_OVERSHOOT_PREC/10.0;
  else
    final_overshoot_prec = FIN_OVERSHOOT_PREC;
  //
  // start main iteration
  //
  tstep=0;// total timesteps
  iter1=iter2=0;// iteration counters
  do{
    /* 

    start backward loop in strain or time

    this loop will decrease the initial timestep until overshoot is within bounds

    */
    if(iter1) {
      /* save the last, unsuccesful time */
      tlast = final_state->t - dt;
    }
    /*
      half the time step 
    */
    dt /= 2.0;
    if(fabs(dt) < EPS_PREC){
      fprintf(stderr,"advect_and_search_for_state: refine 1: dt too small: %g\n",dt);
      *last_dt = dt;
      return -2;
    }
    /*
      
       reset all times and states here 

    */
    copy_state(initial_state,final_state,FALSE,"local initial","local final");
    /* 
       if second iter, advance to previous time 
    */
    if(iter1){
      if(fabs(tlast) > EPS_COMP_PREC){
	if(rk_wrapper(final_state->x,final_state->t,(final_state->t+tlast),
		      0,itrack_strain,FALSE,MDP->rkeps,MDP,FALSE)){
	  *last_dt = tlast;
	  return -10;
	}
	final_state->t += tlast;
      }
    }
    /* 
    
    start  incremental approach

    */
    i1=0;
#ifdef FSTRACK_DEBUG_ITER
    iter_report_progress(final_state,mode,iter1,iter2,i1,i2,dt,overshoot,mstrain,
			 bailout);
#endif
    while((!bailout_crit_fulfilled(final_state,mode,model,&overshoot,bailout,&mstrain))
	  &&(i1 < MAXITER1)){
      /* 

      backward (or forward for FORWARD_STRAIN) advection in steps
      until criterion fullfilled

      */
      if(rk_wrapper(final_state->x,final_state->t,(final_state->t+dt),
		    0,itrack_strain,FALSE,MDP->rkeps,MDP,FALSE)){
	*last_dt = dt;
	return -10;
      }
      final_state->t += dt;
      i1++;
      tstep++;
#ifdef FSTRACK_DEBUG_ITER
      iter_report_progress(final_state,mode,iter1,iter2,i1,i2,dt,overshoot,mstrain,
			   bailout);
#endif
    } /* end backward loop */
#ifdef FSTRACK_DEBUG_ITER
    iter_report_progress(final_state,mode,iter1,iter2,i1,i2,dt,overshoot,mstrain,
			   bailout);
#endif
    if(i1 == MAXITER1){
      fprintf(stderr,"advect_and_search_for_state: too many advection steps (%i) at dt: %g t: %g\n",
	      i1,dt,final_state->t);
      *last_dt = dt;
      return -3;
    }
    iter1++;
    /* 
       
    end overshoot loop

    */
  }while((overshoot > INIT_OVERSHOOT_PREC)&&(iter1 < MAXITER2));
  if(iter1 == MAXITER2){
    fprintf(stderr,"advect_and_search_for_state: couldn't find an initial timestep dt: %g iter: %i\n",
	    dt,iter1);
    *last_dt = dt;
    return -4;
  }
  while((overshoot > final_overshoot_prec) && (iter2 < MAXITER4)){
    /* 

    start refinement loop if overshoot isn't small enough 

    */
    iter2++;
    // remove last timestep 
    dt = -dt;
    tstep++;
    if(rk_wrapper(final_state->x,final_state->t,(final_state->t+dt),
		  0,itrack_strain,FALSE,MDP->rkeps,MDP,FALSE)){
      *last_dt = dt;
      return -10;
    }
    final_state->t += dt;
    /*
      decrease dt and flip back to negative again
    */
    if(fabs((dtnew = -dt/2.0)) < EPS_PREC)
      dt= -EPS_PREC;
    else
      dt = dtnew;
    /*
      
      go backward until bailout reached again

    */
    i2=0;
    while((!bailout_crit_fulfilled(final_state,mode,model,&overshoot,bailout,
				   &mstrain)) && (i2 < MAXITER3)){
      if(rk_wrapper(final_state->x,final_state->t,(final_state->t+dt),
		    0,itrack_strain,FALSE,MDP->rkeps,MDP,FALSE)){
	*last_dt = dt;
	return -10;
      }
      final_state->t += dt;
      i2++;
      tstep++;
#ifdef FSTRACK_DEBUG_ITER
      iter_report_progress(final_state,mode,iter1,iter2,
			   i1,i2,dt,overshoot,mstrain,
			   bailout);
#endif
    } /* end backward loop */
    if(i2 == MAXITER4){
      fprintf(stderr,"advect_and_search_for_state: too many refinement advection steps (%i) at dt: %g\n",i2,dt);
      *last_dt = dt;
      return -6;
    }
  } /* end refinement loop */
  if(iter2 == MAXITER2){
    fprintf(stderr,"advect_and_search_for_state: couldn't refine stepping dt: %g iter: %i\n",
	    dt,iter2);
    *last_dt = dt;
    return -7;
  }
  if(verbose)
    fprintf(stderr,"advect_and_search_for_state: overshoot: %8.5f%% fi: %i si: %i ti: %i total: %4i steps at t: %11g bailout: s:%1i d:%1i t:%1i\n",
	    overshoot*100.0,iter1,iter2,i1,tstep,final_state->t,bailout[0],bailout[1],bailout[2]);
  *last_dt = dt;
  return tstep;
}
void iter_report_progress(struct stt *state, int mode,
			  int i1,int i2,int i3,int i4,
			  COMP_PRECISION dt, 
			  COMP_PRECISION overshoot,
			  COMP_PRECISION mstrain,
			  my_boolean *bailout)
{
  COMP_PRECISION xp[3];
  lonlatz_from_xp(state->x,xp);
  if((mode == BACKWARD_STRAIN)||(mode == FORWARD_STRAIN)){
    fprintf(stderr,"bs: i1:%2i i2:%2i i3:%5i i4:%2i t:%11g dt:%11g loc: %11g, %11g, %11g, os: %11g MS: %11g b:%i%i%i\n",
	    i1,i2,i3,i4,state->t,dt,xp[FSTRACK_X],xp[FSTRACK_Y],xp[FSTRACK_Z],overshoot,mstrain,
	    bailout[0],bailout[1],bailout[2]);
  }else if(mode == BACKWARD_DEPTH){
    fprintf(stderr,"bd: i1:%2i i2:%2i i3:%5i i4:%2i t:%11g dt:%11g loc: %11g, %11g, %11g, os: %11g MS: %11g b:%i%i%i\n",
	    i1,i2,i3,i4,state->t,dt,xp[FSTRACK_X],xp[FSTRACK_Y],xp[FSTRACK_Z],overshoot,mstrain,
	    bailout[0],bailout[1],bailout[2]);

  }
}
