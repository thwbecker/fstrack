#include "fstrack_flow.h"
/*

  Lyapunov exponent by going foward in time until threshold for convergence
  reached 

  on output, each tracer will have only one state, corresponding to the initial 
  position and the time it took for convergence

  $Id: calc_lyapunov.c,v 1.8 2004/04/22 14:10:44 becker Exp becker $


*/
void calc_lyapunov(struct mod *model,int atracer)
{
  // time step and max  time
  EXT_FTRN_PREC ti=0.0,tf= MDP->tmax;
  struct stt *state;
  state=(struct stt *)calloc(1,sizeof(struct stt));
  /* 
     allocate and initialize the state 
  */
  init_state(state,0.0,MDP);
  /*  */
#ifdef FSTRACK_DEBUG
  if(MDP->remove_strain){
    PE("calc_lyapunov: error: remove strain set during Lyapuov mode");
    exit(-1);
  }
#endif
  //
  // start with zero sum of log(alpha), 
  //
  /* this will be used internally for the Exponents, will it? */
  zero_vector(MDP->rlyap,3);
  //
  // use initial position
  //
  copy_state(&model->tracer[atracer].state[0],state,
	     FALSE,"tracer 0","local initial");
  //
  // advect tracer until Lyapunov bailout criteria is reached
  //
  rk_wrapper(state->x,ti,tf,0,TRUE,TRUE,MDP->rkeps,MDP,FALSE);
  //
  if(fabs(state->t - tf)<EPS_PREC){
    fprintf(stderr,"calc_lyapunov: warning: tracer %i: bailout at max time %g\n",
	    atracer,MDP->tmax);
    model->tracer[atracer].discarded = TRUE;
  }else{
#ifdef FSTRACK_DEBUG
    fprintf(stderr,"calc_lyapunov: tracer %i: bailout at time %g\n",
	    atracer,state->t);
#endif
    ;
  }
  //
  // copy the Lyapunov exponent sums
  //
  copy_3dvector(MDP->rlyap,model->tracer[atracer].state[0].isa);
  model->tracer[atracer].state[0].t = state->t;
  free(state);
}
