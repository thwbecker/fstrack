#include "fstrack_flow.h"

/*

  calculate if local_state is above a threshold
  if so, overshoot will hold the relative amount of overshoot,
  if not, overshoot will be zero

  bailout(3): boolean vector with 

  bailout[0]: strain
  bailout[1]: depth
  bailout[2]: time 
  
  bailout flags

  if strain is calculated, will return the local max FSE strain

  $Id: bailoutcrit.c,v 1.16 2010/12/31 18:59:17 becker Exp $

*/
my_boolean bailout_crit_fulfilled(struct stt *local_state, int mode, 
				  struct mod *model,COMP_PRECISION *overshoot,
				  my_boolean *bailout,COMP_PRECISION *mstrain)
{
  COMP_PRECISION finv[9],dovershoot,tovershoot,sovershoot,eval[3];
  static my_boolean init=FALSE;
  static COMP_PRECISION bttime;// backward target time
  /* max strain upwards */
  *mstrain = 0.0;
  // bailout flags for all criteria
  bailout[0] = bailout[1] = bailout[2] = FALSE;
  /* init overshoot */
  *overshoot = 0.0;
  if(!init){
    bttime = model->itime - model->tf;// this is normally simply -model->tf
    init = TRUE;
  }
  switch (mode){
  case BACKWARD_DEPTH:// if tracer is below model->sdepth, bailout
    if(local_state->x[FSTRACK_R] <= model->sdepth){
      *overshoot = (model->sdepth - local_state->x[FSTRACK_R])/model->sdepth;
      bailout[1] = TRUE;
    }
    break;
  case BACKWARD_STRAIN:
    dovershoot = tovershoot = sovershoot = 0.0;// overshoots for the criteria
    /* 

       max difference between eigenvalues that would be accumulated if
       the tracer were to advect from the local position back to it's
       origin if (*mstrain) is reached, or tracer deeper than
       STRAINMAX_DEPTH_BAILOUT

    */
    /*
      
      check if the tracer is below a certain depth STRAINMAX_DEPTH_BAILOUT
      (e.g. 410km) such that strain should be zeroed out

     */
    if(local_state->x[FSTRACK_R] <= STRAINMAX_DEPTH_BAILOUT){// below depth boundary
      dovershoot = (STRAINMAX_DEPTH_BAILOUT - local_state->x[FSTRACK_R])/
	STRAINMAX_DEPTH_BAILOUT;
      /* set the depth bailout flag */
      bailout[1] = TRUE;
    }
    /*

      check for time history longer than cutoff time bstime (which is always > 0)

    */
#ifdef FSTRACK_DEBUG
    if(local_state->t > 0){
      fprintf(stderr,"bailout_crit: error, expecting negative times in backward strain calcuation\n");
      exit(-1);
    }
#endif
    if(local_state->t <= -model->bstime){// time goes backward
      tovershoot = -(model->bstime+local_state->t)/model->bstime;
      /* set the time bailout flag */
      bailout[2] = TRUE;
    }
    /*
      
      determine the finite strain we would have if this tracer were to 
      advect back to its initial position forward in time
      
      use the inverse of the deformation matrix which would transform the
      local state back to a final state with deformation F^-1
      
    */
    invert3x3c((local_state->x+3),finv);
    // obtain max EV of L2 strain
    *mstrain = max_strain_from_def(finv,eval);
    if(*mstrain >= model->maxstrain){
      // strain is above threshold
      sovershoot = (*mstrain - model->maxstrain)/model->maxstrain;
      /* 
	 set the strain bailout flag 
      */
      bailout[0] = TRUE;
    }
    if(bailout[0] || bailout[1] || bailout[2]){
      // select the largest overshoot assuming that non-threshold types are zero
      *overshoot = MAX(sovershoot, tovershoot);
      *overshoot = MAX(*overshoot, dovershoot);
    }
    break;
  case FORWARD_STRAIN:
    dovershoot = tovershoot = sovershoot = 0.0;// overshoots for the criteria
    /* 
       same for forward strain
    */
    if(local_state->x[FSTRACK_R] <= STRAINMAX_DEPTH_BAILOUT){// below depth boundary
      dovershoot = (STRAINMAX_DEPTH_BAILOUT - local_state->x[FSTRACK_R])/
	STRAINMAX_DEPTH_BAILOUT;
      bailout[1] = TRUE;
    }
#ifdef FSTRACK_DEBUG
    if(local_state->t < 0){
      fprintf(stderr,"bailout_crit: error, expecting positive times in forward strain calculation\n");
      exit(-1);
    }
#endif
    if(local_state->t > model->bstime){
      tovershoot = (-model->bstime+local_state->t)/model->bstime;
      bailout[2] = TRUE;
    }
    a_equals_b_vector(finv,(local_state->x+3),9);
    *mstrain = max_strain_from_def(finv,eval);
    if(*mstrain >= model->maxstrain){
      sovershoot = (*mstrain - model->maxstrain)/model->maxstrain;
      bailout[0] = TRUE;
    }
    if(bailout[0] || bailout[1] || bailout[2]){
      *overshoot = MAX(sovershoot, tovershoot);
      *overshoot = MAX(*overshoot, dovershoot);
    }
    break;			/* end forward strain */
  case BACKWARD_TIME:
#ifdef FSTRACK_DEBUG
    if(local_state->t > model->itime){
      fprintf(stderr,"bailout_crit_fulfilled: error: time (%g) should be smaller than itime (%g)\n",
	      local_state->t,model->itime);
      exit(-1);
    }
#endif
    if(local_state->t <= bttime){// time goes backward
      *overshoot = (bttime - local_state->t)/model->tf;
#ifdef FSTRACK_DEBUG
      if(*overshoot < 0){
	fprintf(stderr,"bailout_crit_fulfilled: error: backwar time: overshoot sign error\n");
	exit(-1);
      }
#endif
      bailout[2]=TRUE;
    }
    break;
  default:
    fprintf(stderr,"bailout_crit_fulfilled: error: mode %i undefined\n",mode);
    exit(-1);
    break;
  }
  if(bailout[0] || bailout[1] || bailout[2])
    return TRUE;
  else
    return FALSE;
  

}
