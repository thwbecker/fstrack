#include "fstrack_flow.h"
/*

  advect a single tracer by calling rk_wrapper

  nsteps is the number of steps to output (besides the initial one)

  $Id: advect.c,v 1.39 2010/12/31 18:59:17 becker Exp $

*/


//#define SUPER_FSTRACK_DEBUG
void advect(struct trc *tracer, struct mod *model, 
	    int nsteps,int mode)
{
  COMP_PRECISION dt,last_dt = 0.0;
  my_boolean itrack_strain,eval_strain,rk_error=FALSE; /* rk_error is for RK */
  struct stt *state;
  int i,errcode;		/* errcode is for iteration */
  if(nsteps < 1){
    fprintf(stderr,"advect: need at least one step, nstep: %i\n",
	    nsteps);
    exit(-1);
  }
  /* 
     clear the local states 
  */
  state = (struct stt *)calloc(NR_LOC_STATES,sizeof(struct stt));
  for(i=0;i < NR_LOC_STATES;i++){
    init_state((state+i),0.0,MDP);
  }
  /* 
     texture part
  */
  texture_pathline_init(MDP,
			pressure_model(ZDEPTH((tracer->state+0)->x[FSTRACK_R])));
  /* 
     save the state of the tracer on entry
  */
  compute_state_parameters((tracer->state+0),MDP,"tracer 0");
  copy_state((tracer->state+0),(state+ORIG),TRUE,
	     "tracer 0","local orig");  
#ifdef FSTRACK_DEBUG
  fprintf(stderr,"advect: original location x: %g %g %g\n",
	  state[ORIG].x[0],state[ORIG].x[1],state[ORIG].x[2]);
#endif
  /* 

  deal with the various options for tracer advection

  */
  switch(mode){
  case FORWARD_TIME:
    /* 
       

    advect tracer from time=0 to time=model->tf in nsteps+1 
    steps (including starting step) and trace strain
    
    these dt steps are subdivided appropriately in the RK routine
    
    */
    itrack_strain = TRUE;
    //
    // initial position and strain
    //
    tracer->state[0].t = model->itime;
    /* save initial state in final_state, to be modified  */
    copy_state((tracer->state+0),(state+FINAL),TRUE,
	       "tracer 0","local final");
    /* time step */
    dt = model->tf / (COMP_PRECISION)(nsteps);
    if(!finite(dt)){
      fprintf(stderr,"advect: input error, tf: %g dt: %g steps: %i\n",model->tf,
	      dt,nsteps);
      exit(-1);
    }
#ifdef FSTRACK_DEBUG
    fprintf(stderr,"advect: advecting tracer from %g to %g in %i steps\n",
	   model->itime,model->tf,nsteps);
#endif
    for(i=0;i < nsteps;i++){
      //
      // advect and trace strains from t_i=t to t_f=t+dt
      // and possibly evolve texture, while saving in final_state
      //
      rk_error = rk_wrapper(state[FINAL].x,
			    state[FINAL].t,
			    (state[FINAL].t + dt),
			    0,itrack_strain,FALSE,MDP->rkeps,
			    MDP,TRUE);
      if(rk_error)
	break;
      //
      // update time and strain parameters
      //
      state[FINAL].t += dt;
      compute_state_parameters((state+FINAL),MDP,
			       "local final");
      /* 
	 add new state to tracer and copy the current state to it
      */
      add_state(tracer,state[FINAL].t,MDP);
      copy_state((state+FINAL),(tracer->state+i+1),
		 TRUE,"local final","tracer i");
    }
    // if this loop terminates after one step, t = tf
#ifdef FSTRACK_DEBUG
    if((!rk_error)&&(fabs(state[FINAL].t - model->tf) > MDP->rkeps)){
      fprintf(stderr,"advect: error, final times do not match: %g vs. %g, diff: %e\n",
	      state[FINAL].t,model->tf,state[FINAL].t - model->tf);
      exit(-1);
    }
#endif
    break;
  case BACKWARD_TIME:		/* 
				   move back in time by some interval, and
				   then forward to initial position, 
				   accumulating strain
				*/
  case BACKWARD_TIME_DEPTH:	/* same as above, but also use max depth */
  case BACKWARD_DEPTH:		/* 
				   go back to some depth, and then forward,
				   accumulating strain 
				*/
  case BACKWARD_STRAIN:		/* 
				   go back to some depth and time, so that the final
				   strain is identical to a critical strain
				*/
  case FORWARD_STRAIN:		/* 
				   go forward to some depth and time,
				   so that the final strain is
				   identical to a critical strain
				*/
#ifdef FSTRACK_DEBUG
    fprintf(stderr,"advect: in advect_and_search for type branch\n");
#endif
    /*

      advect tracer to sdepth or to max strain and depth, backward in time, or
      backward to depth and trace them back to origin

    */
    if((mode == BACKWARD_STRAIN)||(mode == FORWARD_STRAIN))
      itrack_strain = TRUE;  // need the strain during iteration for bailout
    else
      itrack_strain = FALSE; /* 
				other modes don't need strain and
				should be quicker without 
			     */
    /* prepare the first state */
    tracer->state[0].t = model->itime;
    /* copy over to intial_state */
    copy_state((tracer->state+0),(state+INITIAL),TRUE,
	       "tracer 0","local initial");
    /* 
       
       go back in time in return as state[FINAL] the origin of the tracer 
       which we will need to follow to arrive at the orig_state such that 
       the tracers fulfills the bailout criteria

    */
#ifdef FSTRACK_DEBUG
    errcode = advect_and_search_for_state(model,mode,(state+INITIAL),
					  (state+FINAL),TRUE,
					  itrack_strain,&last_dt);
#else
    errcode = advect_and_search_for_state(model,mode,(state+INITIAL),
					  (state+FINAL),FALSE,
					  itrack_strain,&last_dt);
#endif
    if(errcode < 0){
      /* 
	 
      iteration error, couldn't find a bailout state for this tracer, discard 
      
      */
      deal_with_discarded_tracer(tracer,state[INITIAL].x);
      return;
    }
    /* 
       all texture modes need to be recomputed, same thing if we need
       a specific number of intermediate output steps. if not, we may 
       use the original path
    */
    if((MDP->texture_mode != NO_TEXTURE) || 
       (!itrack_strain) || (nsteps > 1))
      eval_strain = TRUE;
    else
      eval_strain = FALSE;
    if(eval_strain){
      /*

      if the search for the bailout state didn't involve tracking the
      strain or if we want intermediate tracer states along the way,
      or if we would like to follow texture development
      
      go forward in time for this tracer from newly determined initial
      state while tracing the accumulating strain
      
      */
      if(mode != FORWARD_STRAIN){
	state[INITIAL].t = state[FINAL].t;
	/* new time of initial tracer state */
	clear_state(&state[INITIAL],state[INITIAL].t);
	/* new origin of tracer */
	copy_3dvector(state[FINAL].x,state[INITIAL].x);
	/* 
	   time step
	*/
	dt = (model->itime - state[FINAL].t) / (COMP_PRECISION)(nsteps);
      }else{
	/* 
	   forward strain 
	*/
	state[INITIAL].t = model->itime;
	dt = (-model->itime + state[FINAL].t) / (COMP_PRECISION)(nsteps);
	clear_state(&state[INITIAL],state[INITIAL].t);
	copy_state((state+ORIG),(state+INITIAL),TRUE,
		   "state orig","state initial");
      }
      if(dt <= 0){
	fprintf(stderr,"advect: error: dt negative = %g (%g - %g) in nstep > 1 tracking\n",
		dt,model->itime,state[FINAL].t);
	exit(-1);
      }
      /* 
	   texture part
      */
      /* (re-)compute new strains at initial position */
      compute_state_parameters((state+INITIAL),MDP,"local initial");
      /* the new initial state will be the first tracer state */
      copy_state((state+INITIAL),(tracer->state+0),TRUE,"local initial","tracer 0");
      /* 
	 now perform the forward computation, marchign initial_state
      */
#ifdef FSTRACK_DEBUG
      fprintf(stderr,"advect: re-forward advection branch from t: %g to %g, x0: %g, %g, %g\n",
	      state[INITIAL].t,state[FINAL].t,state[INITIAL].x[0],state[INITIAL].x[1],
	      state[INITIAL].x[2]);
#endif
      for(i=0;i < nsteps;i++){
	//                             t_i           t_f
	rk_error = rk_wrapper(state[INITIAL].x,state[INITIAL].t,
			      (state[INITIAL].t+dt),
			      fabs(last_dt), TRUE,FALSE,MDP->rkeps,MDP,TRUE);
	if(rk_error)
	  break;
	/* update and save */
	state[INITIAL].t += dt;
	/* compute strains and such */
	compute_state_parameters((state+INITIAL),MDP,"local initial");
	/* add to tracer */
	add_state(tracer,state[INITIAL].t,MDP);
	copy_state((state+INITIAL),(tracer->state+i+1),TRUE,
		   "local initial","tracer i");
      }
      if(rk_error)
	break;
    }else{
      /* 
	 
	 we are only interested in the initial and final steps of this
	 tracer (nsteps == 1) and have been tracking its strain
	 already (itrack_strain == TRUE) therefore, we can just flip
	 the states, taking care to invert the deformation matrix


	 this trick doesn't work for texture
	 
      */
#ifdef FSTRACK_DEBUG
      if(tracer->nstate != 1){
	fprintf(stderr,"advect: error, in flip mode trick: nstate != 1 (%i)\n",
		tracer->nstate);
	exit(-1);
      }
#endif
      if(mode != FORWARD_STRAIN){
	/* 
	   
	backward parts

	*/
	// origin at depth and zero strain
	clear_state(&tracer->state[0],state[FINAL].t);
	// start location
	copy_3dvector(state[FINAL].x,tracer->state[0].x);
	// final position at original location and cumulative strain
	add_state(tracer,0.0,MDP);
	copy_3dvector(state[ORIG].x,tracer->state[1].x);
	/* 
	   inverted deformation of state[FINAL] at depth will be 
	   the deformation of the last tracer state at the original 
	   shallow depth 
	*/
	invert3x3c((state[FINAL].x+3),(tracer->state[1].x+3));
      }else{
	/* 
	   forward strain part
	*/
	/* start location */
	copy_3dvector(state[INITIAL].x,tracer->state[0].x);
	// final position 
	add_state(tracer,0.0,MDP);
	a_equals_b_vector(tracer->state[1].x,state[FINAL].x,12);
      }
      /* 
	 compute strain  and other parameters
      */
      compute_state_parameters((tracer->state+1),MDP,"tracer second");
    }
    break;
  case ISA:
    fprintf(stderr,"advect: error: don't call with ISA parameter, call calc_isa striaght\n");
    exit(-1);
    break;
  default:
    fprintf(stderr,"advect: mode %i undefined\n",mode);
    exit(-1);
    break;
  }
  if(rk_error){
    /* 
       we have encountered an internal RK error somewhere along the way
    */
    deal_with_discarded_tracer(tracer,state[ORIG].x);
    return;
  }
#ifdef FSTRACK_DEBUG
  print_tracer_stats(tracer,0,stderr); /* first */
  print_tracer_stats(tracer,1,stderr); /* last */
#endif
  /* discard the local state variables */
  if(MDP->drex_save_pole){
    for(i=0;i<NR_LOC_STATES;i++)
      free(state[i].pdens);
  }
  free(state);
}



/*
  
  wrapper for fortran Runge Kutta routine
  
  return TRUE, if error occured, else FALSE, if everything is OK

  if activate_texture is set, this will operate on drex->x, and copy x[12]
  to this vector

*/
my_boolean rk_wrapper(COMP_PRECISION *x,// state vector x[calc_strain?12:3]
		      COMP_PRECISION ti, // initial time
		      COMP_PRECISION tf, // final time
		      COMP_PRECISION dt,/* suggested absolute value
					   of initial timestep. ifset
					   to zero, will pick
					    something that might or
					    might not be reasonable
					*/
		       
		      my_boolean calc_strain,   // track strains? (or only position) 
		      my_boolean calc_lya, /* calc Lyapunov exponent */
		      COMP_PRECISION rkeps, // eps precision for RK routine
		      struct der_par *dp, /* derivative structure */
		      my_boolean activate_texture) /* if desired, activate texture 
						       computations */
{
  static my_boolean warned=FALSE;
  /* 
     maximum timestep for initial guess and during 
     RK advection 
  */
  my_boolean error = FALSE,reassign=FALSE;
  COMP_PRECISION alpha[3],*xloc=NULL,tspan,dtmin,ttmp1,ttmp2;
  int n,i,err_code;
#ifdef FSTRACK_DEBUG
  static my_boolean dumped_kr_init_state=FALSE;
  FILE *out;
  if(x[FSTRACK_R] > 1.0){
    fprintf(stderr,"rk_wrapper: x_r > 1 (%g) on entry\n",x[FSTRACK_R]);
    exit(-1);
  }
  if(activate_texture && (!calc_strain))
    PEE("rk_wrapper: error: texture but no strain computation");
#endif
  if(dp->remove_strain && (!warned)){
    PE("");
    PE("rk_wrapper: WARNING: removing symmetric part of the vel. grad. matrix");
    PE("");
    warned = TRUE;
  }
  /* 
     check the advection timespan
  */
  tspan = tf - ti;
  if(fabs(tspan) < EPS_COMP_PREC){
    fprintf(stderr,"rk_wrapper: error: on input: ti: %g tf: %g ti-tf: %e\n",
	    ti,tf,tspan);
    return -111;
  }
  if(dt < 0){
    /* dt should be >= 0 */
    fprintf(stderr,"rk_wrapper: error: dt is negative on entry: %g\n",
	   dt);
    exit(-1);
  }
  if(fabs(dt) < EPS_PREC){
    /* 
       adjust the initial guess dt automatically 
    */
    ttmp1 = fabs(tspan);
    dt = ttmp1     / 50.0;      /* (tf-ti)/50 is the first guess */
    ttmp2 = ttmp1 / 5.;	        /* this is the lower limit */
    dtmin = MIN(1e-4,ttmp2);	/* or 1e-4, whichever is smaller */
    dt = MAX(dtmin,dt);
  }
  dt = MIN(dp->hmax,dt);	/* limit to maxstep */
  dt *= MY_SIGN(tspan);
  //
  //set parameters and call odeint
  //
  if (!calc_strain){
    //
    // no strain tracing and vel gradient matrix needed, just velocity advection
    //
    n = 3;
    if(calc_lya){
      fprintf(stderr,"rk_wrapper: You cannot calculate lyapunov exponents\n");
      fprintf(stderr,"rk_wrapper: without calculating strain (ifellipse)!\n");
      fprintf(stderr,"rk_wrapper: ifellipse,iflyap: %i, %i\n",
	      calc_strain,calc_lya);
      exit(-1);
    }
  }else{	
    /* need to compute strain */
    n = 12;
  }

  /* 
     do we want nr_odeint to store kmax 
     intermediate results at each dxsav interval?
  */
  dp->kmax = 0;			/* no storage */
  dp->dxsav = 0.1;		/* some time interval */
  /* 

  check some limits

  */
#ifdef FSTRACK_DEBUG
  if((x[FSTRACK_THETA] < dp->dtheta) || (PI-x[FSTRACK_THETA] < dp->dtheta)){
    fprintf(stderr,"rk_wrapper: near polar on entry: lon: %g lat: %g r: %g, dlat: %g\n",
	    PHI2LONGITUDE(x[FSTRACK_PHI]),THETA2LATITUDE(x[FSTRACK_THETA]),
	    ZDEPTH(x[FSTRACK_R]),PHI2LONGITUDE(dp->dtheta));
  }
#endif
  /* 

  are we dealing with texture?

  */
  /* regular setting, will use the state for forward integration */
  if(activate_texture){
    switch(dp->texture_mode){
    case NO_TEXTURE:
      reassign = FALSE;
      xloc = x;
      break;
    case KR_TEXTURE:
      if(!calc_strain){
	PEE("rk_wrapper: error: need strain computation for texture");
      }
      reassign = TRUE;
      n = dp->drex->bsize;		/* this number include odf,odf_ens,acs,acs_ens */
      /* 
	 copy the state x to the first 12 entries of drex->x array
      */
      a_equals_b_vector(dp->drex->x,x,12);
      /* 
	 use the large drex array for integration 
      */
      xloc = dp->drex->x;
#ifdef FSTRACK_DEBUG
      if(!dumped_kr_init_state){
	out = myopen("krinit.dat","w","rk_wrapper");
	print_vector(xloc,n,out);
	fclose(out);
	fprintf(stderr,"rk_wrapper: written the %i component KR init state to krinit.dat\n",
		n);
	fprintf(stderr,"rk_wrapper: indices: x: %i F: %i odf: %i odf_ens: %i acs: %i acs_ens: %i\n",
		0,3,dp->drex->podf,dp->drex->podf_ens,
		dp->drex->pacs,dp->drex->pacs_ens);
	dumped_kr_init_state = TRUE;
      }
#endif
      /* end kaminski and ribe init */
      break;
    default:
      PEE("rk_wrapper: wrong texture mode");
      break;
    }    
  }else{			/* no texture computation */
    xloc = x;
  }
  /* 
     make sure that tracer coordinates are OK, 
     possibly also make the texture parameters physical
  */
  rk_check_phys_limit(xloc,n,dp,FALSE);
  //
  // integrate to obtain position and deformation matrix at t2
  //
#ifdef FSTRACK_DEBUG
  fprintf(stderr,"rk_wrapper: ti: %11.4e tf: %10.4e dt: %10.3e cstrain: %i clya: %i atexture: %i eps: %g dt_max: %g\n",
	  ti,tf,dt,calc_strain,calc_lya,activate_texture,rkeps,
	  dp->hmax);
#endif  
  err_code = nr_odeint(xloc,n,ti,tf,rkeps,dt,EPS_COMP_PREC,
		       dp->hmax,&dp->nok,&dp->nbad,dp);
  //
  if(err_code != 0){
    /* error handling */
    switch(err_code){
    case -2:
      fprintf(stderr,"rk_wrapper: too many steps in nr_odeint\n");
      break;
    case -1:
      fprintf(stderr,"rk_wrapper: stepsize underflow in nr_odeint\n");
      break;
    default:
      fprintf(stderr,"rk_wrapper: error code %i in nr_odeint\n",
	      err_code);
      break;
    }
    error = TRUE;
  }else
    error = FALSE;
  if(!error){
    if (calc_lya){
      /* 
	 Lyapunov part 
      */
      // normalize
      gramschmidt(xloc,alpha);
      // and increment lyapunov exponents
      for(i=0;i < 3;i++)
	dp->rlyap[i] += log(alpha[i]);
    }
    //
    //     adjust the location vector, x(1...3)
    //
    rk_check_phys_limit(xloc,n,dp,FALSE);
    /* 
       copy back to state, since that gets passed out of this routine
       as output
    */
    if(reassign)
      a_equals_b_vector(x,xloc,12);
#ifdef FSTRACK_DEBUG
    if(fabs(dp->tbailout - tf) > EPS_PREC){
      fprintf(stderr,"rk_wrapper: timing error, final time (%g) != tbailout (%g)\n",
	      tf,dp->tbailout);
      exit(-1);
    }
#endif
  }
  return error;
}

/* 
   
given a state, compute some of the parameters that derive from it, such
as the left-stretch strain, and possibly texture derived parameters


*/
//#define SUPER_FSTRACK_DEBUG
void compute_state_parameters(struct stt *state,struct der_par *dp,
			      char *vname)
{
#ifdef SUPER_FSTRACK_DEBUG
  COMP_PRECISION strain,val[3];
  int i,j;
  FILE *out;
  char fname[200];
#endif
  COMP_PRECISION z,dummy = 0.0;
  static my_boolean warned[3] = {FALSE,FALSE,FALSE};
  /*  
      compute the left stretch strain
  */
  calc_l2_left_stretch((state->x+3),state->left_stretch);
  if(dp->texture_mode != NO_TEXTURE){
    /* 
       TEXTURE COMPUTATION
    */
    /* 
       compute the average stiffness matrix from volume averaging the
       texture grains. sav gets returned C style

    */
    switch(dp->drex_module){
    case EMOD_PT:
      /* 
	 pressure and temperature inferred from depth of tracer
      */
      if(!warned[0]){
	fprintf(stderr,"advect: WARNING: using T, p varying modules\n");
	warned[0] = TRUE;
      }
      /* get inferred temperature and compute tensors at depth and
	 such */
      z = ZDEPTH(state->x[FSTRACK_R]);
      drex_vhr(dp->drex,state->Sav,DREX_CTP_ESTEY,
	       temperature_model(z,1),pressure_model(z),
	       dp->vhr_favg);
      break;
    case EMOD_PT_NEW:
      /* 
	 pressure and temperature inferred from depth of tracer using
	 the new derivatives
      */
      if(!warned[0]){
	fprintf(stderr,"advect: WARNING: using T, p varying modules with new derivatives\n");
	warned[0] = TRUE;
      }
      /* get inferred temperature and compute tensors at depth and
	 such */
      z = ZDEPTH(state->x[FSTRACK_R]);
      drex_vhr(dp->drex,state->Sav,DREX_CTP_NEW,
	       temperature_model(z,1),pressure_model(z),
	       dp->vhr_favg);
      break;
    case EMOD_FIXPT:
      /* fixed pressure, temperature */
      if(!warned[0]){
	fprintf(stderr,"advect: WARNING: using fixed T: %gK , p: %g Gpa varying modules\n",
		dp->drex_module_t,dp->drex_module_p);
	warned[0] = TRUE;
      }
      drex_vhr(dp->drex,state->Sav,DREX_CTP_ESTEY,dp->drex_module_t,
	       dp->drex_module_p,dp->vhr_favg);
      break;
    case EMOD_CONST:
      /* constant, surface (reference) conditions */
      drex_vhr(dp->drex,state->Sav,DREX_CTP_CONST_KR,
	       dummy,dummy,dp->vhr_favg);
      break;
    default:
      fprintf(stderr,"advect: error, elastic module mode %i undefined\n",
	      dp->drex_module);
      exit(-1);
    }
    /* 
       should we compute avg axes and the olivine J index? 
    */
    if(dp->drex_compute_avg_axes){
      if(!warned[1]){
	if(!dp->rotate_grain_coord_sys)
	  fprintf(stderr,"compute_state_parameters: WARNING: not rotating [100] stats with xp, global Cartesian\n");
	else
	  fprintf(stderr,"compute_state_parameters: rotating [100] stats with xp, local Cartesian\n");
	warned[1] = TRUE;
      }
      drex_compute_avg_axes(dp->drex,state->x,dp->rotate_grain_coord_sys);
      /* 
	 save the stats for the olivine [100] unweighted and weighted axes 
      */
      copy_vector(dp->drex->avg_axes,state->olivine_axes_stats,6); 
      copy_vector(dp->drex->dev_axes,(state->olivine_axes_stats+6),2); 
      /* olivine [100] J index */
      state->olivine_axes_stats[8] = dp->drex->olivine_jindex;
    }
    /* if needed, save the ODF for the pole figures */
    if(dp->drex_save_pole){
      /* 
	 
      compute and save an orientation density function for a pole figure
      from the current DREX structure in a tracer state
      
      */
      /* 
	 in general, we want to rotate the pole figures with the flow such
	 that north is always up, and east is right. for constant velocity
	 gradient matrices, the global Caretsian framework might be more
	 appropriate
      */
#ifdef FSTRACK_DEBUG
      if(dp->constant_vgm_for_debugging){
	if(dp->rotate_grain_coord_sys == TRUE)
	  fprintf(stderr,"compute_state_parameters: WARNING: constant VGM overrides rotation of poles\n");
	dp->rotate_grain_coord_sys = FALSE;
      }
#endif
      if(!warned[2]){
	if(!dp->rotate_grain_coord_sys)
	  fprintf(stderr,"compute_state_parameters: WARNING: not rotating ODF with xp, global Cartesian\n");
	else
	  fprintf(stderr,"compute_state_parameters: rotating ODF with xp, local Cartesian\n");
	warned[2] = TRUE;
      }
      /* 
	 compute the pole figure densities for three axes and two
	 components at location x in spherical system
      */
      drex_compute_pdens(dp->drex,state->x,state->t,
			 dp->rotate_grain_coord_sys,
			 dp->print_o1xyz,dp->lpo_tracer_count,
			 FALSE);
      /* 
	 save the orientation density function
      */
      if(state->npole != dp->drex->npole){
	fprintf(stderr,"compute_state_parameters: error: state->npole: %i vs. drex->npole: %i\n",
		state->npole,dp->drex->npole);
	exit(-1);
      }
      copy_vector(dp->drex->pdens,state->pdens,dp->drex->npole); 
#ifdef SUPER_FSTRACK_DEBUG
      strain = max_strain_from_left_stretch(state->left_stretch,val);
      fprintf(stderr,"compute_state_parameters: computed pdens[%i] for %s at time: %g strain: %g\n",
	      dp->drex->npole,vname,state->t,strain);
      /* write pdens to file */
      sprintf(fname,"o1pdens.2.%0g.data",state->t);
      fprintf(stderr,"compute_state_parameters: writing olivinite [100] at time %g to %s\n",
	      state->t,fname);
      out=fopen(fname,"w");
      for(i=0;i  < dp->drex->np[1];i++){/* lat loop */
	for(j=0;j < dp->drex->np[0];j++){ /* lon loop */
	  fprintf(out,"%.3e\n",state->pdens[i*dp->drex->np[0]+j]);
	}
      }
      fclose(out);
#endif
    }
  }
}



  
void deal_with_discarded_tracer(struct trc *tracer,COMP_PRECISION *orig_location)
{
  COMP_PRECISION xp[3];
  /* set the tracer's discarded flag */
  tracer->discarded = TRUE;
  /* 
     output 
  */
  lonlatz_from_xp(orig_location,xp);
  fprintf(stderr,"advect: RK error, tracer discarded: lon0: %g lat0: %g z0: %g\n",
	  xp[FSTRACK_X],xp[FSTRACK_Y],xp[FSTRACK_Z]);
#ifdef FSTRACK_DEBUG
  PEE("advect: exiting since in debug mode");
#endif 
}

/* 
   
initialize the texture parameters for each new computation

this routine should be called anew each time a new tracer is tracked


*/
void texture_pathline_init(struct der_par *dp,
			   COMP_PRECISION pressure)
{
  static my_boolean init=FALSE;
  COMP_PRECISION tau_def[4];
 
  /* 

     texture part

  */
  switch(dp->texture_mode){
  case NO_TEXTURE:		/* 
				   no action, if no texture is to be
				   followed
				*/
    ;
    break;
  case KR_TEXTURE:
    if(!init){
      /* 
	 INITIALIZE KR VARIABLES
      */
      /* 

	 first time parameter init for Kaminski and Ribe code

	 drex_assign_
      */
      /* get slip system */
      drex_assign_tau_from_system(tau_def,dp->drex_type,
				  dp->drex_rss_tau,pressure,
				  dp->drex_sstp);
      /* init */
      drex_initialize(&dp->drex,0,dp->drex_size3,
		      dp->drex_Xol,tau_def,DREX_TAU_ENS_DEF,
		      dp->drex_Mob,dp->drex_chi,dp->drex_lambda,
		      dp->drex_np,dp->drex_start_grains_oriented,
		      FALSE);
      fprintf(stderr,"init_der_par: initialized Kaminski & Ribe texture mode\n");
      fprintf(stderr,"init_der_par: KR parameters: size3: %i lambda: %g Mob: %g chi: %g\n",
	      dp->drex->size3,dp->drex->lambda,dp->drex->Mob,dp->drex->chi);
      fprintf(stderr,"init_der_par: slip systems: %g %g %g %g\n",
	      tau_def[0],tau_def[1],tau_def[2],tau_def[3]);
      fprintf(stderr,"init_der_par: Voigt Reuss Hill averging: %g (%s)\n",
	      dp->vhr_favg,(dp->vhr_favg==0)?("Voigt"):(dp->vhr_favg==1?("Reuss"):("VHR")));
      if(dp->drex_save_pole)
	fprintf(stderr,"init_der_par: saving [abc] axes orientations\n");
      init = TRUE;
    } /* end first time call */
    /* 
       
    the following gets executed each time this routine is called:
    
    assign the random starting values and initial density distribution
    for both olivine and enstatite components

    */
    drex_init_pathline(dp->drex);
    dp->lpo_tracer_count++;
    break;
  default:
    fprintf(stderr,"advect: error: texture mode %i undefined\n",
	    dp->texture_mode);
    exit(-1);
  }
 
}
