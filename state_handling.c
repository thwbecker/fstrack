/* 

routines that deal with a state of a tracer

$Id: state_handling.c,v 1.3 2004/04/21 11:41:55 becker Exp becker $

*/

#include "fstrack_flow.h"
/*

  add a state to a tracer structure at a certain time,
  and initialize that state

*/
void add_state(struct trc *tracer,COMP_PRECISION time,
	       struct der_par *dp)
{
  int astate;
  astate = tracer->nstate;// index of the new state
  // increment state counter
  tracer->nstate++;
  /* make sure it'll fit */
  if(tracer->nstate >= MY_USHRT_MAX){
    fprintf(stderr,"add_state: error, too many states: %i, nstate variable limit is %i\n",
	    tracer->nstate, MY_USHRT_MAX);
    exit(-1);
  }
  // alocate more room
  if(tracer->nstate == 1)
    tracer->state = (struct stt *)malloc(sizeof(struct stt));
  else
    tracer->state = (struct stt *)realloc(tracer->state,
					  (size_t)tracer->nstate*
					  sizeof(struct stt));
  if(!tracer->state)
    MEMERROR("add_state");
  /* 
     initialize the state 
  */
  init_state((tracer->state+astate),time,dp);
}


//
// initialize a tracer state (this can only be called once)
//
void init_state(struct stt *state,COMP_PRECISION time,
		struct der_par *dp)
{
  if(dp->drex_save_pole){
    state->npole = dp->drex_npole;
    my_vecalloc(&state->pdens,state->npole,"init_state");
  }else{
    state->npole = 0;
    state->pdens = NULL;
  }
  clear_state(state,time);
}
/* 
   clear the tracer state, this can be called repeatedly
*/
void clear_state(struct stt *state,COMP_PRECISION time)
{
  int i;
  COMP_PRECISION *zero_flt;
  state->pipar = 0.0;
  state->isa_exists = FALSE;
  for(i=0;i<3;i++){
    state->isa[i] = state->x[i] = 0.0;
  }
  for(i=0;i<36;i++){
    state->Sav[i] = 0.0;
  }
  for(i=0;i<MAX_OL_AXES_ENTRIES;i++)
    state->olivine_axes_stats[i] = 0.0;
  //
  // initialize deformation with identity matrix, f[i] is x[3+i]
  // diagonal components of F are 0, 4, and 8, thus
  //
  unity_matrix((state->x+3),3);
  //
  // initialize strain state (L2) (this is again the unity matrix, but hey)
  calc_l2_left_stretch((state->x+3),state->left_stretch);
  //
  // init time
  //
  state->t = time;
  //
  if(state->npole){	
    /* if we are carrying the pole figure densities around */
    /* 
       for texture 
    */
    zero_flt=(COMP_PRECISION *)calloc(state->npole,sizeof(COMP_PRECISION));
    if(!zero_flt)
      MEMERROR("clear_state: zero_flt");
    copy_vector(zero_flt,state->pdens,state->npole);
    free(zero_flt);
  }
}
/*

  copy a state from a onto state b, i.e.

  b = a

*/
void copy_state(struct stt *a,struct stt *b,
		my_boolean copy_texture,
		char *aname, char *bname)
{
  copy_vector(a->left_stretch,b->left_stretch,9);
  copy_vector(a->x,b->x,12);
  b->t = a->t;
  b->pipar = a->pipar;
  copy_vector(a->isa,b->isa,3);
  b->isa_exists = a->isa_exists;
  copy_vector(a->Sav,b->Sav,36);
  copy_vector(a->olivine_axes_stats,b->olivine_axes_stats,MAX_OL_AXES_ENTRIES);
  b->npole = a->npole;
  if(copy_texture){
    /* copy the ODF for texture, if present */
    copy_vector(a->pdens,b->pdens,b->npole); 
#ifdef SUPER_DEBUG
    if(b->npole){
      fprintf(stderr,"copy_state: copying from %s to %s (<pdens[%i]>: %g)\n",
	      aname,bname,b->npole,std(b->pdens,b->npole));
    }else{
      fprintf(stderr,"copy_state: copying from %s to %s (no pdens)\n",
	      aname,bname);
    }
#endif
  }
}
