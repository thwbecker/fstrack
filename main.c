#include "fstrack_flow.h"
/*


  reads velocities, traces particles, and keeps track of finite strain 
  or texture following Kaminski and Ribe approach

  uses GMT grd files to read in velocities at steady state, or at different
  timesteps 
  
  performs Runge Kutta integration to follow stream lines and compute
  finite strain 
  
  type 'fstrack -h' for man page and explanations

  $Id: main.c,v 1.33 2010/12/31 18:59:17 becker Exp $

*/

int main(int argc, char **argv)
{
  int i,j,nout,itout;
  char ofname[STRLEN],prefix[400];
  struct mod *model;
  char    timebuf[2000];
  time_t  t_now,t_start;
  my_boolean zero_boundary_vr = TRUE;
  fprintf(stderr,"%s: version  $Id: main.c,v 1.33 2010/12/31 18:59:17 becker Exp $\n",
	  argv[0]);
  fprintf(stderr,"%s: compiled on %s at %s\n",
	  argv[0],__DATE__,__TIME__);
  model=(struct mod *)calloc(1,sizeof(struct mod));// allocate with zeroes!
  //
  // set the default parameters
  set_defaults(model);
  //
  // check for command line parameters
  check_for_parameters(argc,argv,model);
  //
  // read velocities (for this to work, the two routines above have to been called)
  
  read_vel_grids(model,model->verbose,zero_boundary_vr);

  // test input
  // alayer=20;
  // write_ascii_layer(model->vt,alayer,model,"vt","vt.dat",model->velscale);
  // initialize the tracer field
  initialize_tracers(model);
  if(!model->ntracer){
    PE("main: no tracers generated, exciting");
    exit(-1);
  }
  if(MDP->strain_fraction_from_gamma)
    sprintf(prefix,"tracer.er");	/* using er.i.grd alpha grids */
  else
    sprintf(prefix,"tracer");

  itout = (int)(model->ntracer / 10.0)+1;	/* output of time each itout */
  // output of some tracer paths each nout
  if(!model->ntout){
    fprintf(stderr,"main: increasing tracer sampling to unity\n");
    model->ntout = 1;
  }
  nout = model->ntracer/model->ntout;
  if(nout == 0)nout = 1;	/* make sure we're not dividing by zero */
  /* 
     switch modes
  */
  switch(model->amode){
  case FORWARD_TIME:
  case BACKWARD_DEPTH:
  case BACKWARD_STRAIN:
  case FORWARD_STRAIN:
  case BACKWARD_TIME:
    t_start = time(NULL);strftime(timebuf,79,"%a, %F, %H:%M:%S",localtime(&t_start));
    fprintf(stderr,"\nmain: starting main tracer loop at %s\n",timebuf);
    for(i=0;i < model->ntracer;i++){
      /*
	
	advect tracer to bailout criterion, and save nsteps during the advection
	path
	
      */
      advect((model->tracer+i),model,model->nsteps,model->amode);
      if(i%nout==0){
	/* 

	output of tracer specific quantities along path

	*/
	// write history
	sprintf(ofname,"hist.%i.%s",i,model->ostring);
	write_tracer_history((model->tracer+i),ofname,LOCATION_LLZ,
			     model->verbose,model);
	// deformation tensor
	sprintf(ofname,"def.%i.%s",i,model->ostring);
	write_tracer_history((model->tracer+i),ofname,DEFORMATION_TENSOR,
			     model->verbose,model);
	// left stretch tensor components
	sprintf(ofname,"strain.%i.%s",i,model->ostring); 
	write_tracer_history((model->tracer+i),ofname,STRAIN_COMP,
			     model->verbose,model);
	if(0){
	  // cartesian left stretch strain components
	  sprintf(ofname,"xyztL2c.%i.%s",i,model->ostring); 
	  write_tracer_history((model->tracer+i),ofname,
			       XYZ_STRAIN_CART_COMP,
			       model->verbose,model);
	}
	if(1){
	  // velocity gradient matrix
	  sprintf(ofname,"vgm.%i.%s",i,model->ostring); 
	  write_tracer_history((model->tracer+i),ofname,LOCATION_VGM_RTPTVGM,
			       model->verbose,model);
	}
	if(1){
	  // velocity gradient matrix, rtp system for coords, matrix cartesian
	  sprintf(ofname,"vgm.sc.%i.%s",i,model->ostring); 
	  write_tracer_history((model->tracer+i),ofname,LOCATION_VGM_RTPTVGM_CART,
			       model->verbose,model);
	}
	if(0){
	  // velocity gradient matrix, xyz system for coords, matrix cartesian
	  sprintf(ofname,"vgm.cc.%i.%s",i,model->ostring); 
	  write_tracer_history((model->tracer+i),ofname,LOCATION_VGM_XYZTVGM_CART,
			       model->verbose,model);
	}
	// sqrt(eigenvalues of left stretch strain) and eigenvectors
	sprintf(ofname,"seval.%i.%s",i,model->ostring);
	write_tracer_history((model->tracer+i),ofname,STRAIN_EVAL,
			     model->verbose,model);
	if(0){
	  // sqrt(eigenvalues of left stretch strain) and eigenvectors, all in cartesian system
	  sprintf(ofname,"seval.c.%i.%s",i,model->ostring);
	  write_tracer_history((model->tracer+i),ofname,STRAIN_EVAL_CART,
			       model->verbose,model);
	}
	if(1){
	  /* strain-rate and vorticity  */
	  sprintf(ofname,"sv.%i.%s",i,model->ostring);
	  write_tracer_history((model->tracer+i),ofname,STRAINRATE_VORTICITY,
			       model->verbose,model);

	}


	if(0){
	  /* 
	     compute the ISA axes for this tracer 
	  */
	  for(j=0;j<model->tracer[i].nstate;j++){
	    /* loop through all states */
	    calc_isa((model->tracer+i),j,model);
	  }
	  /* output */
	  sprintf(ofname,"isa.%i.%s",i,model->ostring);	  
	  write_tracer_history((model->tracer+i),ofname,ISA_AXES,
			       model->verbose,model);
	}
	if(MDP->texture_mode != NO_TEXTURE){
	  /* 
	     
	  texture was computed
	  
	  */

	  if(model->sav_out){
	    /* 
	       stiffness tensor 
	    */
	    if(MDP->drex_module !=  EMOD_CONST)
	      /* variable version, this changes the filename */
	      sprintf(ofname,"savd.%i.%s",i,model->ostring); 
	    else
	      sprintf(ofname,"sav.%i.%s",i,model->ostring); 
	    write_tracer_history((model->tracer+i),ofname, 
				 TEXTURE_SAV,model->verbose,model); 
	  }
	  /* 
	     
	  best fit axis for tranverse isotropy, percentage 
	  anisotropy derived from SAV, and l/n ratio (vsv^2/vsh^2)
	  
	  */
	  if(MDP->drex_module != EMOD_CONST)
	    sprintf(ofname,"tid.%i.%s",i,model->ostring);
	  else
	    sprintf(ofname,"ti.%i.%s",i,model->ostring);
	  write_tracer_history((model->tracer+i),ofname,
			       TEXTURE_TI,model->verbose,model);
	  /* 
	     Rayleigh wave 2phi and 4phi terms 
	  */
	  if(MDP->drex_module != EMOD_CONST)
	    sprintf(ofname,"rphid.%i.%s",i,model->ostring);
	  else
	    sprintf(ofname,"rphi.%i.%s",i,model->ostring);
	  write_tracer_history((model->tracer+i),ofname,
			       TEXTURE_RAYLEIGH_RPHI,
			       model->verbose,model);
	  if(MDP->drex_compute_avg_axes){
	    /* statistics for a,b,c axes */
	    sprintf(ofname,"polestat.%i.%s",i,model->ostring);
 	    write_tracer_history((model->tracer+i),ofname, 
 				 TEXTURE_STATS,model->verbose,model); 
	  }
	  if(MDP->drex_save_pole){
	    /* 
	       pole figure information 
	    */
	    sprintf(ofname,"pole.%i.%s",i,model->ostring);
 	    write_tracer_history((model->tracer+i),ofname, 
 				 TEXTURE_ODF,model->verbose,model); 
	  }
	}
      }
      if((i%itout) == 0){
	t_now = time(NULL);
	strftime(timebuf,79,"%a, %H:%M",localtime(&t_now));
	fprintf(stderr,"main: done with tracer %6i out of %6i at %s, avg %11g s/tracer\n",
		i+1,model->ntracer,timebuf,
		difftime(t_now,t_start)/(COMP_PRECISION)(i+1));
      }
    } /* end i tracer loop */
    /* 



    main output section


    */
    //
    // print out tracers by depth levels (layers)
    //
    for(j=0,i=0;i < model->ntlevels;i++){
      //
      // write initial tracer positions to file
      //
      sprintf(ofname,"%s.i.l.%g.%s",prefix,ZDEPTH(model->tlevel[i]),model->ostring);
      write_tracer_field(model,j,j+model->ntd[i]-1, 0,ofname,ALL_TRACER_LLZAT,
			 model->verbose);
      //
      // write final tracer positions 
      //
      sprintf(ofname,"%s.f.l.%g.%s",prefix,ZDEPTH(model->tlevel[i]),model->ostring);
      write_tracer_field(model,j,j+model->ntd[i]-1,-1,ofname,ALL_TRACER_LLZAT,
			 model->verbose);
      //
      // final tracer left stretch upper right half strain components
      //
      sprintf(ofname,"%s.%s.%g.%s",prefix,TSFILE,ZDEPTH(model->tlevel[i]),model->ostring);
      write_tracer_field(model,j,j+model->ntd[i]-1,-1,ofname,
			 ALL_TRACER_STRAIN_COMP,model->verbose);
      if(0){
	//
	// final tracer L strain eigenvalues and vectors as angles
	//
	sprintf(ofname,"%s.f.e.%g.%s",prefix,ZDEPTH(model->tlevel[i]),model->ostring);
	write_tracer_field(model,j,j+model->ntd[i]-1,-1,ofname,
			   ALL_TRACER_STRAIN_EIGEN,model->verbose);
      }
      //
      // final tracer all F deformation matrix 
      //
      sprintf(ofname,"%s.%s.%g.%s",prefix,TFFILE,ZDEPTH(model->tlevel[i]),model->ostring);
      write_tracer_field(model,j,j+model->ntd[i]-1,-1,ofname,
			 ALL_TRACER_DEFORMATION,model->verbose);
    

      if((MDP->texture_mode != NO_TEXTURE)&&(MDP->rotate_grain_coord_sys)){
	/* 

	LPO fabrics were computed and we are rotating into a
	geographic coordinate system (local S-E-U system). if not, we
	are likely debugging and don't want to confuse things


	*/
	/* 
	   best fit transverse isotropy axis from texture SAV and vsh/vsv ratio
	*/
	if(MDP->drex_module != EMOD_CONST)
	  sprintf(ofname,"%s.tid.%g.%s",prefix,ZDEPTH(model->tlevel[i]),model->ostring);
	else
	  sprintf(ofname,"%s.ti.%g.%s",prefix,ZDEPTH(model->tlevel[i]),model->ostring);
	write_tracer_field(model,j,j+model->ntd[i]-1,-1,ofname,
			   ALL_TRACER_TI,model->verbose);
	if(MDP->drex_module != EMOD_CONST)
	  sprintf(ofname,"%s.rphid.%g.%s",prefix,ZDEPTH(model->tlevel[i]),model->ostring);
	else
	  sprintf(ofname,"%s.rphi.%g.%s",prefix,ZDEPTH(model->tlevel[i]),model->ostring);
	write_tracer_field(model,j,j+model->ntd[i]-1,-1,ofname,
			   ALL_TRACER_RPHI,model->verbose);
	if(model->sav_out){
	    /* 
	       stiffness tensor 
	    */
	  if(MDP->drex_module != EMOD_CONST)
	    sprintf(ofname,"%s.savd.%g.%s",prefix,ZDEPTH(model->tlevel[i]),model->ostring);
	  else
	    sprintf(ofname,"%s.sav.%g.%s",prefix,ZDEPTH(model->tlevel[i]),model->ostring);
	  write_tracer_field(model,j,j+model->ntd[i]-1,-1,ofname,
			     ALL_TRACER_SAV,model->verbose);
	}
      }
      // increment tracer offset counter for each layer
      j += model->ntd[i];
    }
    if(model->ntracer < 10000){
      // all tracer histories
      sprintf(ofname,"%s.hist.%s",prefix,model->ostring);
      write_all_tracer_history(model,ofname,TRUE);
      // all tracer histories as xyz
      sprintf(ofname,"%s.hist.xyz.%s",prefix,model->ostring);
      write_all_tracer_history(model,ofname,FALSE);
    }else{
      PE("main: skipping tracer history output");
    }
    break;
  case ISA:
    /* 
       compute the ISA axis
    */
    PE("main: running main tracer loop for ISA");
    for(i=0; i < model->ntracer;i++){ 
      /* compute the ISA for all tracers, use state 0*/
      calc_isa((model->tracer+i),0,model);
    }
    for(j=0,i=0;i<model->ntlevels;i++){	/* output  */
      sprintf(ofname,"%s.%s.%g.%s",prefix,TISAFILE,
	      ZDEPTH(model->tlevel[i]),model->ostring);
      write_tracer_field(model,j,j+model->ntd[i]-1,-1,ofname,
			 ALL_TRACER_ISA,model->verbose);
      j += model->ntd[i];
    }
    PE("main: ISA done");
    break;
  case CALC_LYAPUNOV:
    if(MDP->texture_mode != NO_TEXTURE)
      PEE("main: error: cannot compute Lyapunov exponents and texture at same time");
    // trace Lyapunov exponents, in this case L2 of each tracer will have the 
    // exponents
    for(i=0;i<model->ntracer;i++)
      calc_lyapunov(model,i);
    
    for(j=0,i=0;i<model->ntlevels;i++){
      // write tracer positions and lyapunov exponents to file
      sprintf(ofname,"tracer.f.lya.%g.%s",ZDEPTH(model->tlevel[i]),model->ostring);
      write_tracer_field(model,j,j+model->ntd[i]-1, 0,ofname,ALL_TRACER_LYA,
			 model->verbose);
      j += model->ntd[i];
    }
    break;
  default:
    PEE("main: opmode error");
    break;
  }
  PE("main: done");
  return (0);
}



