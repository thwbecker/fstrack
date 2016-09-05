#include "fstrack_flow.h"
/*
  
  deal with default parameters and command line options to change
  those 

  $Id: input_para.c,v 1.37 2010/12/29 19:09:55 twb Exp becker $

*/
//
// set the defaults, can be changed later by command line parameters
//
void set_defaults(struct mod *model)
{
  static my_boolean init=FALSE;
  if(init){
    fprintf(stderr,"set_defaults: error: called more than once\n");
    exit(-1);
  }
  /* 
     
  initialize the derivative structure , see deriv.c for defaults
  
  */
  init_der_par(&model->dp);
  /*  

  initialize the rest of the model structure
  
  */
  model->tf=1.0;// default exit time  for forward and backward time modes
  model->itime=0.0;// default initial time for forward advection
  model->sdepth = P410_RL;// exit depth for backward depth mode
  model->maxstrain = MSTRAIN_BAILOUT_DEF;/* strain exit threshold */
  model->bstime = MTIME_BAILOUT_DEF;// max time bailout for backward strain
  /*  */
  model->amode=BACKWARD_TIME;// default advection mode
  model->tinitmode=DIST_EVEN_AREA;// default tracer init mode

  model->ntracer=DEF_NR_TRACERS;// nr of tracers
  model->nsteps=1;// steps of advection for tracers whose path is output
  model->ntout = 2;// output of ntour tracer paths
  model->tdepth_from_file=FALSE;
  model->init_tracer_depth = INIT_TRACER_DEPTH_DEF;
  /* use an initial strain? (normally off) */
  model->use_initial_strain = FALSE;
  // io 
  model->read_gmt=TRUE;
  sprintf(model->ostring,"dat"); /* output file suffix */
  sprintf(model->tdepthfilename,"%s.%s",TDFILE,model->ostring);// tracer depth levels
  model->sav_out = FALSE;	/* default for stiffness tensor */
 
  sprintf(model->sw_sens_file,"%s",SW_SENS_FILE);// without $HOME/
 
  //
  // ISA related, see Kaminski & Ribe (2002) and structures.h
  //
  model->pi_crit = ISA_PI_CRIT_DEF;// change of evec factor
#ifdef FSTRACK_DEBUG
  model->verbose = TRUE;
#else
  model->verbose = FALSE;// default is not verbose for output
#endif
  init = TRUE;
}
/*

  check for input options, syntax is explained in help page printout further down

*/
#define INCHECK(x) {if(x != 1){fprintf(stderr,"check_for_parameters: error reading parameter after option %s\n",argv[i-1]);exit(-1);}}

void check_for_parameters(int argc, char **argv, struct mod *model)
{
  int i,j,k;
  COMP_PRECISION tmp_dbl;
  FILE *in;
  for(i=1;i<argc;i++){
    if((strcmp(argv[i],"-h")==0) || (strcmp(argv[i],"-?")==0) || 
       (strcmp(argv[i],"-help")==0)){
      phelp(argv,model);
    }else if(strcmp(argv[i],"-ft")==0)
      model->amode = FORWARD_TIME;
    else if(strcmp(argv[i],"-l")==0)
      model->amode = CALC_LYAPUNOV;
    else if(strcmp(argv[i],"-time")==0){
      advance_argument(&i,argc,argv);
      INCHECK(sscanf(argv[i],FLT_FORMAT,&model->tf));
    }else if(strcmp(argv[i],"-itime")==0){
      advance_argument(&i,argc,argv);
      INCHECK(sscanf(argv[i],FLT_FORMAT,&model->itime));
    }else if(strcmp(argv[i],"-bstime")==0){
      advance_argument(&i,argc,argv);
      INCHECK(sscanf(argv[i],FLT_FORMAT,&model->bstime));
    }else if(strcmp(argv[i],"-bs")==0){
      model->amode = BACKWARD_STRAIN;
    }else if(strcmp(argv[i],"-fs")==0){
      model->amode = FORWARD_STRAIN;
    }else if(strcmp(argv[i],"-v")==0){
      model->amode = ONLY_VEL_STATS;
    }else if(strcmp(argv[i],"-verbose")==0){
      model->verbose = TRUE;
    }else if(strcmp(argv[i],"-bt")==0)
      model->amode = BACKWARD_TIME;
    else if(strcmp(argv[i],"-th")==0)
      MDP->history = TRUE;
    else if(strcmp(argv[i],"-uis")==0)
      model->use_initial_strain = TRUE;
    else if(strcmp(argv[i],"-o1xyz")==0)
      MDP->print_o1xyz = 1;
    else if(strcmp(argv[i],"-o3xyz")==0)
      MDP->print_o1xyz = 2;
    else if(strcmp(argv[i],"-remove_symm_strain")==0)
      MDP->remove_strain = TRUE;
    else if(strcmp(argv[i],"-nt")==0){
      advance_argument(&i,argc,argv);
      INCHECK(sscanf(argv[i],"%i",&model->ntracer));
    }else if(strcmp(argv[i],"-wt")==0){
      advance_argument(&i,argc,argv);
      INCHECK(sscanf(argv[i],"%i",&model->ntout));
    }else if(strcmp(argv[i],"-os")==0){
      advance_argument(&i,argc,argv);
      INCHECK(sscanf(argv[i],"%s",model->ostring));
    }else if(strcmp(argv[i],"-tim")==0){
      advance_argument(&i,argc,argv);
      INCHECK(sscanf(argv[i],"%i",&model->tinitmode));
    }else if(strcmp(argv[i],"-ns")==0){
      advance_argument(&i,argc,argv);
      INCHECK(sscanf(argv[i],"%i",&model->nsteps));
    }else if(strcmp(argv[i],"-vs")==0){
      advance_argument(&i,argc,argv);
      INCHECK(sscanf(argv[i],"%lf",&tmp_dbl));
      fprintf(stderr,"check_for_parameters: multiplying char. velocity by %g\n",
	      tmp_dbl);
      MDP->velscale *= tmp_dbl;
      fprintf(stderr,"check_for_parameters: total velocity scale: %g\n",
	      MDP->velscale);
    }else if(strcmp(argv[i],"-strain")==0){
      advance_argument(&i,argc,argv);
      INCHECK(sscanf(argv[i],FLT_FORMAT,&model->maxstrain));
    }else if(strcmp(argv[i],"-pc")==0){
      advance_argument(&i,argc,argv);
      INCHECK(sscanf(argv[i],FLT_FORMAT,&model->pi_crit));
    }else if(strcmp(argv[i],"-tm")==0){
      advance_argument(&i,argc,argv);
      INCHECK(sscanf(argv[i],"%i",&MDP->texture_mode));
    }else if(strcmp(argv[i],"-pole")==0){
      MDP->drex_save_pole = TRUE;
    }else if(strcmp(argv[i],"-sav_out")==0){
      model->sav_out = TRUE;
    }else if(strcmp(argv[i],"-bd")==0)
      model->amode = BACKWARD_DEPTH;
    else if(strcmp(argv[i],"-isa")==0)
      model->amode = ISA;
    else if(strcmp(argv[i],"-td")==0)
      model->tdepth_from_file = TRUE;
    else if(strcmp(argv[i],"-tdd")==0){
      advance_argument(&i,argc,argv);
      INCHECK(sscanf(argv[i],FLT_FORMAT,&model->init_tracer_depth));
    }else if(strcmp(argv[i],"-vhr")==0){
      advance_argument(&i,argc,argv);
      INCHECK(sscanf(argv[i],FLT_FORMAT,&MDP->vhr_favg));
    }else if(strcmp(argv[i],"-g")==0)
      model->read_gmt = TRUE;
    else if(strcmp(argv[i],"-b")==0)
      model->read_gmt = FALSE;
    else if(strcmp(argv[i],"-nrgcs")==0) /* grain coordinate system */
      MDP->rotate_grain_coord_sys = FALSE;
    else if(strcmp(argv[i],"-type")==0){
      advance_argument(&i,argc,argv);
      INCHECK(sscanf(argv[i],"%i",&MDP->drex_type));
      if(MDP->drex_type == DREX_FABRIC_FREE){
	for(j=0;j<4;j++){
	  advance_argument(&i,argc,argv);
	  INCHECK(sscanf(argv[i],FLT_FORMAT,&MDP->drex_rss_tau[j]));
	}
      }
    }
    else if(strcmp(argv[i],"-sav_noTp")==0) /* constant stiffness tensors,
					       using refernce tensors
					    */
      MDP->drex_module = EMOD_CONST;
    else if(strcmp(argv[i],"-sav_newTp")==0) /* use the new derivatives */
      MDP->drex_module = EMOD_PT_NEW;
    else if(strcmp(argv[i],"-sav_fixTp")==0){ /* fixed at certain p, T */
      MDP->drex_module = EMOD_FIXPT;
      advance_argument(&i,argc,argv);
      INCHECK(sscanf(argv[i++],FLT_FORMAT,&MDP->drex_module_t));
      INCHECK(sscanf(argv[i],FLT_FORMAT,&MDP->drex_module_p));
    }
    else if(strcmp(argv[i],"-er_diff")==0){ /* diffusion creep factors */
      MDP->strain_fraction_from_gamma = TRUE;
    }else if(strcmp(argv[i],"-depth")==0){
      advance_argument(&i,argc,argv);
      INCHECK(sscanf(argv[i],FLT_FORMAT,&model->sdepth));
      model->sdepth=ND_RADIUS(model->sdepth);
    }else if(strcmp(argv[i],"-mob")==0){
      advance_argument(&i,argc,argv);
      INCHECK(sscanf(argv[i],FLT_FORMAT,&MDP->drex_Mob));
    }else if(strcmp(argv[i],"-chi")==0){
      advance_argument(&i,argc,argv);
      INCHECK(sscanf(argv[i],FLT_FORMAT,&MDP->drex_chi));
    }else if(strcmp(argv[i],"-kr_start_oriented")==0){
      MDP->drex_start_grains_oriented = TRUE;
    }else if(strcmp(argv[i],"-xol")==0){
      advance_argument(&i,argc,argv);
      sscanf(argv[i],FLT_FORMAT,&MDP->drex_Xol);
    }else if(strcmp(argv[i],"-sstp")==0){
      advance_argument(&i,argc,argv);
      sscanf(argv[i],FLT_FORMAT,&MDP->drex_sstp);
    }else if(strcmp(argv[i],"-lambda")==0){
      advance_argument(&i,argc,argv);
      sscanf(argv[i],FLT_FORMAT,&MDP->drex_lambda);
    }else if(strcmp(argv[i],"-eps")==0){
      advance_argument(&i,argc,argv);
      INCHECK(sscanf(argv[i],FLT_FORMAT,&MDP->rkeps));
      fprintf(stderr,"check_for_parameters: setting RK eps to %e\n",
	      MDP->rkeps);
    }else if(strcmp(argv[i],"-hmax")==0){
      advance_argument(&i,argc,argv);
      INCHECK(sscanf(argv[i],FLT_FORMAT,&MDP->hmax));
      fprintf(stderr,"check_for_parameters: setting RK max step to %e\n",
	      MDP->hmax);
    }else if(strcmp(argv[i],"-drex_eps")==0){
      advance_argument(&i,argc,argv);
      INCHECK(sscanf(argv[i],FLT_FORMAT,&MDP->drex_epsfac));
      fprintf(stderr,"check_for_parameters: setting DREX eps to %e of RK eps\n",
	      MDP->drex_epsfac);
      if(MDP->drex_epsfac < 0)
	fprintf(stderr,"check_for_parameters: this implies max timestep is %g/strain_rate\n",
		-MDP->drex_epsfac);
    }else if(strcmp(argv[i],"-size3")==0){
      advance_argument(&i,argc,argv);
      INCHECK(sscanf(argv[i],"%i",&MDP->drex_size3));
    }else if(strcmp(argv[i],"-sws")==0){ /* filename for surface wave sensitivity*/
      advance_argument(&i,argc,argv);
      INCHECK(sscanf(argv[i],"%s",model->sw_sens_file));
    }else if(strcmp(argv[i],"-cvgm")==0){ /* constant velocity gradient matrix */
      MDP->constant_vgm_for_debugging = TRUE;
      /* read in S, T, W */
      for(j=0;j < 3;j++){
	advance_argument(&i,argc,argv);
	INCHECK(sscanf(argv[i],FLT_FORMAT,(MDP->cvfd_stw+j)));
      }
    }else{
      fprintf(stderr,"check_for_parameters: can not use parameter %s, use -h for help page\n",
	      argv[i]);
      exit(-1);
    }
  }
  //
  // select the right input file for tracer init modes
  //
  sprintf(model->tdepthfilename,"%s.%s",TDFILE,model->ostring);// tracer depth levels
  switch(model->tinitmode){
  case SPOTTED_LAT_FROM_FILE:
    sprintf(model->tlatfilename,"%s.%s",TRLFILE,model->ostring);// tracer lon lat init filename
    break;
  case SPOTTED_3D_FROM_FILE:
    sprintf(model->tlatfilename,"%s.%s",TRLZFILE,model->ostring);// tracer lon lat z init filename
    break;
  case SPOTTED_3D_WITH_ATTR:
    sprintf(model->tlatfilename,"%s.%s",TRLZAFILE,model->ostring);// tracer lon lat z attr init filename
    break;
  default:
    strcpy(model->tlatfilename,"dummy");// should not be used anywhere in this case
    break;
  }
  //
  // print a message on what we are intending to do
  //
  switch(model->amode){
  case FORWARD_TIME:
    fprintf(stderr,"check_for_parameters: advecting tracers forward in time from %g to %g\n",
	    model->itime,model->itime+model->tf);
    break;
  case BACKWARD_TIME:
    fprintf(stderr,"check_for_parameters: advecting tracers back in time by %g and forward to init. pos.\n",
       model->tf);
    break;
  case BACKWARD_STRAIN:
    PE("check_for_parameters: advecting tracers back to mstrain and forward to init. pos.");
    break;
  case FORWARD_STRAIN:
    PE("check_for_parameters: advecting tracers forward to mstrain");
    break;
  case BACKWARD_DEPTH:
    PE("check_for_parameters: advecting tracers back to sdepth and forward to init. pos.");
    break;
  case CALC_LYAPUNOV:
    PE("check_for_parameters: calculating Lyapunov exponents");
    break;
  case ONLY_VEL_STATS:
    PE("check_for_parameters: calculating RMS stats of velocity field");
    break;
  case ISA:
    PE("check_for_parameters: calculating infinite strain axes");
    break;
  default:
    fprintf(stderr,"check_for_parameters: mode %i undefined\n",model->amode);
    exit(-1);
    break;
  }
  if(model->bstime < 0 ){
    fprintf(stderr,"check_for_parameters: backward strain time threshold should be > 0 and in Myrs\n");
    exit(-1);
  }
  if(model->ntout < 0)
    PEE("check_for_parameters: ntout< 0 doesn't make sense");
  if((model->amode != CALC_LYAPUNOV)&&(model->amode != ISA)){
  fprintf(stderr,"check_for_parameters: thresholds: tf: %g mstrain: %g tdepth: %g\n",
	    model->tf,model->maxstrain,ZDEPTH(model->sdepth));
    if(model->amode == BACKWARD_STRAIN)
      fprintf(stderr,"check_for_parameters: additional constraints for backward strain: d_max: %g t_max: %g\n",
	      ZDEPTH(STRAINMAX_DEPTH_BAILOUT),model->bstime);
  }else if(model->amode == CALC_LYAPUNOV)
    fprintf(stderr,"check_for_parameters: Lyapunov thresholds: renorm: %g tmax: %g rlbailout: %g\n",
	    MDP->renorm,MDP->tmax,MDP->rlbailout);
  else{ // ISA
    fprintf(stderr,"check_for_parameters: ISA parameters: pi_crit: %g \n",
	    model->pi_crit);
  }
  switch(MDP->texture_mode){
  case NO_TEXTURE:			/* no texture computation */
    break;
  case KR_TEXTURE:
    fprintf(stderr,"check_for_parameters: using Kaminski & Ribe computation, %i grains total\n",
	    (int)pow(MDP->drex_size3,3));
#ifdef FSTRACK_DEBUG
    MDP->drex_compute_avg_axes = TRUE; /* compute the pole figure stats  */
#endif
    break;
  default:
    fprintf(stderr,"check_for_parameters: texture mode %i undefined\n",
	    MDP->texture_mode);
    exit(-1);
    break;
  }
  if(MDP->history)
    fprintf(stderr,"check_for_parameters: expecting to read time history from %s\n",
	    THFILE);
  if(MDP->rkeps < EPS_PREC){
    fprintf(stderr,"check_for_parameters: WARNING: eps for RK (%g) smaller than internal precision (%g)\n",
	    MDP->rkeps,EPS_PREC);
  }
  /* deal with texture timestepping */
  if(MDP->drex_epsfac < 0){
    MDP->strain_rate_control = TRUE;
    MDP->eps_strain  = - MDP->drex_epsfac;
    MDP->drex_epsfac *= -1;
  }else{
    /* using all component checking */
    MDP->strain_rate_control = FALSE;
    MDP->eps_strain = 1e20;
  }
  if(MDP->remove_strain)
    PE("check_for_parameters: WARNING: removing symmetric part of VGM matrix!!!");
  if(MDP->strain_fraction_from_gamma){
    /* 

    deal with er.i.grd grids for scaling of diffusion vs. dislocation creep
    the code expects 

    */
    
    PE("check_for_parameters: determining symmetric part of VGM matrix from er.i.grds");
#ifdef FSTRACK_USE_GGRD			/* initialize er.i.grd for gamma */
    if(ggrd_grdtrack_init_general(TRUE,GAMMA_FILE,DFILE,"-fg",MDP->ggrd_alpha,TRUE,TRUE,FALSE)){
      PE("check_for_parameters: could not initialize er grids for gamma factors");
      exit(-1);
    }
    ggrd_grdtrack_rescale(MDP->ggrd_alpha,FALSE,TRUE,FALSE,0.0);	/* take beta = log10(gamma) */
    /* 
       rescale to alpha = 1/(1+1/beta) which is ~0 for diffusion, and
       ~1 for dislocation creep */
    for(i=0;i < (MDP->ggrd_alpha)->nz;i++){
      k = i * (MDP->ggrd_alpha)->mm;
      for(j=0;j < (MDP->ggrd_alpha)->mm;j++,k++){
	(MDP->ggrd_alpha)->f[k] = 1.0/(1.0+1.0/(MDP->ggrd_alpha)->f[k]);
	if(!finite((MDP->ggrd_alpha)->f[k])){
	  PE("check_for_parameters: error: alpha init not finite");
	  exit(-1);
	}
      }
    }
#else
    PE("check_for_parameters: need to compile with hc/ggrd package for this option");
#endif
  }
  if(model->use_initial_strain){
    fprintf(stderr,"check_for_parameters: WARNING: reading initial strain (spherical system) from %s\n",
	    INIT_STRAIN_FILE);
    in = myopen(INIT_STRAIN_FILE,"r","input_para");
    for(i=0;i<9;i++)
      if(fscanf(in,FLT_FORMAT,(model->initial_strain+i))!=1){
	fprintf(stderr,"input_para: read error: file %s for init strain, component %i\n",
		INIT_STRAIN_FILE,i+1);
	exit(-1);
      }
    fclose(in);
  }
  if(MDP->drex_save_pole){
    fprintf(stderr,"input_para: saving pole figure [] axes\n");
    /* also compute averages in this case */
    MDP->drex_compute_avg_axes = TRUE;
    if(MDP->texture_mode == NO_TEXTURE){ /* test, if texture actually computed */
      fprintf(stderr,"input_para: logic error: saving pole densities but no texture mode set (-tm)\n");
      exit(-1);
    }
    
  } 
  model->use_sw_sens = TRUE;	/* in this case, read also 
				   the sensitivity kernels */
  /* 
     read in sensitivity kernels 
  */
  read_sw_sens(model);
}
// check, if we can read in another values for option argv[i]
void advance_argument(int *i,int argc, char **argv)
{
  if(argc <= *i + 1){// no arguments left
    fprintf(stderr,"%s: input parameters: error: option ""%s"" needs a value\n",
	    argv[0],argv[*i]);
    exit(-1);
  }
  *i += 1;
}

/*
  
  print out help page

*/
void phelp(char **argv,struct mod *model)
{
  PE("");
  PE("reads 3-D spherical flow fields in binary, either from");
  PE("vr/t/p.*.bin files or netcdf vr.i.grd, vt.i.grd, and vp.i.grd.");
  if(model->read_gmt){
    PE("The default is to read GMT (or netdf) grd files.");
  }else{
    PE("The default is to read our native binary format.");
  }
  PE("These are velocities in [cm/yr] in a radial coordinate system.");
  fprintf(stderr,"Internally things are scaled such that the characteristic time is %g Myr\nand i goes from 1",
	  TIMESCALE/1e6);
  fprintf(stderr," to N. N is the number of entries in file \"%s\" that holds the depths\n",
	  DFILE);
  PE("of the velocity layers in km (positive numbers). Program advects");
  PE("tracers while keeping track of finite strain.");
  PE("");
  PE("Options (defaults are in parantheses):");
  PE("");
  PE("I/O related:");
  fprintf(stderr,"\t-b\t\tread native binary files instead of netcdf GRD files %s\n",(!model->read_gmt)?("(default)"):(""));
  fprintf(stderr,"\t-g\t\tread netcdf GRD files instead on native binary format %s\n",(model->read_gmt)?("(default)"):(""));
  PE("");
  PE("\t-vs\tvalue\tmultiply velocity scale by factor value (1)");
  PE("");
  fprintf(stderr,"\t-th\t\tread velocity history at time intervals as specified in \"%s\",\n",
	  THFILE);
  fprintf(stderr,"\t\t\twhich holds the start and end time in units of %g Myr of each velocity stage\n",
	  TIMESCALE/1e6);
  fprintf(stderr,"\t\t\tthis will look for the velocity grids in directories 1/ ... n/\n");
  fprintf(stderr,"\t\t\twhere n is the number of times that were specified in \"%s\".\n",
     THFILE);
  PE("");
  fprintf(stderr,"\t-v\t\tread velocities, print statistics to \"%s\", then exit.\n",VSFILE);
  PE("");
  fprintf(stderr,"\t-os\tstring\twill append string to all output/input files (%s)\n",
     model->ostring);
  fprintf(stderr,"\t-sws\tstring\twill read surface wave kernels from $HOME/string.50.dat for 50s period (%s)\n",
	  model->sw_sens_file);

  PE("");
  PE("\t-verbose\tbe verbose in reporting progress");
  PE("");
  PE("tracer advection related:");
  PE("");
  fprintf(stderr,"The default operation mode of %s is backward-in-time, which can be modified\n",
	  argv[0]);
  PE("by the following switches (default values given in parentheses):\n");
  fprintf(stderr,"\t-time\tvalue\tsets the absolute time interval for advection (%g)\n",
	  model->tf);
  PE("\t-ft\t\trun advection up to specified time (which can be negative), default is to run -bt style");
  fprintf(stderr,"\t-itime\tvalue\tinitial time for forward, final time for backward advection, default: %g\n",
	  model->itime);
  PE("\t-bt\t\trun advection backward in time (default), then forward to initial location");
  PE("\t\t\t(This is different from -ft with a negative time)");
  PE("\t-bs\t\tfollow tracers back until a strain threshold is reached, then go forward to init pos");
  PE("\t\t\tthe threshold strain is defined in the routine mstrain and selected in fstrack.h");
  fprintf(stderr,"\t\t\tNOTE: will also bailout if tracer is deeper than %g km before strain  accumulation\n",
	  ZDEPTH(STRAINMAX_DEPTH_BAILOUT));
  fprintf(stderr,"\t\t\t      or if it takes more than %g absolute time units to accumulate the strain\n",
	  model->bstime);
  fprintf(stderr,"\t\t\t      (but see -depth and -bstime).\n");
  fprintf(stderr,"\t-strain\tvalue\tset threshold strain to value (%g)\n",
	  model->maxstrain);
  PE("\t-fs\t\tfollow tracers forward in time until a strain threshold is reached (mostly for testing)");
#ifdef  MAX_EVAL_STRAIN
  PE("\t\t\tStrain threshold is calculated from E1-1, to change see fstrack.h");
#elif defined MAX_FRAC_EVAL_STRAIN
  PE("\t\t\tStrain threshold is calculated from max(log(E1/E2),log(E2/E3)), to change see fstrack.h");
#else
  PE("ERROR: no criterion defined during compilation, see fstrack.h");
  exit(-1);
#endif
  PE("\t\t\twhere E1 > E2 > E3 are the eigenvalues of the stretching matrix L");
  fprintf(stderr,"\t-bstime\tvalue\tmodify additional absolute time threshold for -bs mode from default of %g\n",
	  model->bstime);
  PE("\t-bd\t\tadvect backward until threshold depth is reached, then go to init");
  fprintf(stderr,"\t-depth\tvalue\tset threshold depth to value (%g)\n",
	  ZDEPTH(model->sdepth));
  fprintf(stderr,"\t-uis\t\tread in an initial strain matrix from %s (off)\n",
	 INIT_STRAIN_FILE);
  PE("\t\t\twhich will be added to the unity matrix normally used to initiale F");
  PE("");
  PE("\t-l\t\tcalculate the Lyapunov exponents for each tracer");
  fprintf(stderr,"\t\t\tin this case, tracers are advected until tmax %g, renormalizing at %g\n",
	  MDP->tmax,MDP->renorm);
  fprintf(stderr,"\t\t\tand bailing out if changes in one lyapunov exponent are less than %g\n",
	  MDP->rlbailout);
  PE("");
  PE("texture related:");
  PE("\t-isa\t\tcalculate the infinite strain axis and Gamma (GOL) at each point");
  fprintf(stderr,"\t-pc\tvalue\tset pi_crit for determining steady-state for ISA calculation (%g)\n",
	  model->pi_crit);
  fprintf(stderr,"\t-tm\tvalue\tadditionally compute texture (%i)\n",NO_TEXTURE);
  fprintf(stderr,"\t\t\tset to %i for no texture, %i for Kaminski&Ribe approach\n",
	  NO_TEXTURE,KR_TEXTURE);
  fprintf(stderr,"\t-size3\tvalue\tnumber of grains in each direction (%i)\n",
	  MDP->drex_size3);
  fprintf(stderr,"\t-kr_start_oriented\t\tstart with an oriented distribution of grains\n\t\t\t(only for testing purposes, off)\n");
  fprintf(stderr,"\t-pole\t\tsave the pole figure density (off)\n");
  fprintf(stderr,"\t-sav_out\tsave the Voigt averaged stiffness tensor (off)\n");
  fprintf(stderr,"\t-sav_noTp\t\tuse moduli that are constant with T and p, i.e. depth (varying default)\n");
  fprintf(stderr,"\t-sav_fixTp\tvalueT valuep\tuse moduli at T = valueT [K] and p = valuep [Gpa]\n");
  fprintf(stderr,"\t-sav_newTp\t\tuse new p,T derivatives (old default)\n");
  fprintf(stderr,"\t-mob\tvalue\tgrain boundary mobility (%g)\n",
	  MDP->drex_Mob);
  fprintf(stderr,"\t-chi\tvalue\tgrain boundary sliding parameter (%g)\n",
	  MDP->drex_chi);
  fprintf(stderr,"\t-lambda\tvalue\tgrain nucleation parameter (%g)\n",
	  MDP->drex_lambda);
  fprintf(stderr,"\t-xol\tvalue\tolivine fraction in %% (rest enstatite) (%g)\n",
	  MDP->drex_Xol);
  fprintf(stderr,"\t-type\tvalue\tuse type \"type\" slip systems (%i)\n",
	  MDP->drex_type);
  fprintf(stderr,"\t\t\t%i: low  stress, low  water (A-type) slip system\n",DREX_FABRIC_ATYPE);
  fprintf(stderr,"\t\t\t%i: high stress, high water (B-type) slip system\n",DREX_FABRIC_BTYPE);
  fprintf(stderr,"\t\t\t%i: mod. stress, high water (C-type) slip system\n",DREX_FABRIC_CTYPE);
  fprintf(stderr,"\t\t\t%i: high stress, low  water (D-type) slip system\n",DREX_FABRIC_DTYPE);
  fprintf(stderr,"\t\t\t%i: mod. stress, mod. water (E-type) slip system\n",DREX_FABRIC_ETYPE);
  fprintf(stderr,"\t\t\t%i: high pressure     ([001] active) slip system\n",DREX_FABRIC_HIGHP);
  fprintf(stderr,"\t\t\t%i: free format, reads in RSS for (010)[100], (001)[100] (010)[001] (100)[001] \n\t\t\tas the next three parameters\n",
	  DREX_FABRIC_FREE);
  fprintf(stderr,"\t\t\t%i: transition between A and high pressure at pressure sstp (%g) GPa\n",
	  DREX_FABRIC_SWITCH_A_HIGHP,MDP->drex_sstp);
  fprintf(stderr,"\t-sstp\tvalue\tslip system transition pressure (GPa) (%g)\n",
	  MDP->drex_sstp);

  fprintf(stderr,"\t-nrgcs\t\tdo not rotate the grain orientations into an E-N-U system but\n");
  fprintf(stderr,"\t\t\trather leave them in a Cartesian system\n");
  fprintf(stderr,"\t-o1xyz\t\tprint the olivine [100] directions cosines and ODF for each grain\n");
  fprintf(stderr,"\t\t\tat each timestep.\n");
  fprintf(stderr,"\t-o3xyz\t\tprint the olivine [100] [010] and [001] directions cosines and ODF for each grain\n");
  fprintf(stderr,"\t-vhr\tvalue\taveraging: Voigt (0) Reuss (1) VRH (0.5) (%g)\n",
	  MDP->vhr_favg);
  
  PE("");

  PE("strain computation related:");
  PE("\t-remove_symm_strain\n\t\t\twill remove the symmetric part of the velocity gradient matrix");
  PE("\t\t\taccording to rules in fse_deriv.F. You will be left with rotations only! OFF by default.");
  PE("\t-cvgm\tS T W\tuse a constant velocity gradient tensor for debugging. This will read in ");
  PE("\t\t\tthree float values for the S, T, and W parameters as in McKenzie & Jackson.\n");
  PE("\t-er_diff\tuse only a part alpha of the strain-rate tensor, where alpha is determined");
  PE("\t\t\tfrom er.i.grd gamma ratios for log10(eps_dislocation/eps_diffusion)\n");

  PE("");
  PE("\t\t\tWarning: if more than one mode flag is specified, the last will override.");
  PE("");

  PE("RK related:");
  fprintf(stderr,"\t-eps\tvalue\tset Runge Kutta precision to value (%g)\n",
	  MDP->rkeps);
  fprintf(stderr,"\t-hmax\tvalue\tmaximum timestep independent of error in char. time (%g)\n",
	  MDP->hmax);
  fprintf(stderr,"\t-drex_eps\tvalue\tset texture eps to value fraction of Rk eps (%.4e)\n",
	  MDP->drex_epsfac);
  fprintf(stderr,"\t\t\tif value is < 0, will not test for texture components but limit\n");
  fprintf(stderr,"\t\t\t   the maximum RK timestep to -value * 1/max_strain_rate\n");
  PE("");

  PE("tracer related:");
  fprintf(stderr,"\t-wt\tvalue\tset the number of tracers whose path is output to value (%i)\n",
	  model->ntout);
  fprintf(stderr,"\t-nt\tvalue\tset the number of tracers (approximately) to value (%i)\n",
	  model->ntracer);
  fprintf(stderr,"\t-ns\tvalue\tset the number of intervals (w/o initial state) for tracers whose path is output(%i)\n",
	  model->nsteps);
  fprintf(stderr,"\t-td\t\tread the tracer depth levels from \"%s\", else one level at %g (see -tdd)\n",
	  model->tdepthfilename,model->init_tracer_depth);
  fprintf(stderr,"\t-tdd\tvalue\tchange the single layer depth to value (%g)\n",
	  model->init_tracer_depth);
  fprintf(stderr,"\t-tim\tvalue\tset the tracer initialization mode (%i)\n",
	  model->tinitmode);
  PE("\n\t\t\ttracers are distributed at each spefified depth horizontally according to:");
  fprintf(stderr,"\t\t\t%i: distribute \"evenly\", ie. for roughly equal area boxes\n",
	  DIST_EVEN_AREA);
  fprintf(stderr,"\t\t\t%i: distribute with fixed dphi dtheta spacing\n",
	  DIST_EVEN_DX);
  fprintf(stderr,"\t\t\t%i: read lateral tracer (lon lat format) locations from \"%s.%s\"\n",
	  SPOTTED_LAT_FROM_FILE,TRLFILE,model->ostring);
  PE("\n\t\t\ttracers locations are specified in 3-D (this overrides -td settings)");
  fprintf(stderr,"\t\t\t%i: read 3-D tracer (lon lat z format) location from \"%s.%s\"\n",
	  SPOTTED_3D_FROM_FILE,TRLZFILE,model->ostring);
  fprintf(stderr,"\t\t\t%i: read tracers with %i attribute(s) (lon lat z attr...)  from \"%s.%s\"\n",
	  SPOTTED_3D_WITH_ATTR,NR_T_ATTR,TRLZAFILE,model->ostring);
  PE("\t\t\t\tnote that attributes are simply tracer labels and could be dealt with externally");
  PE("");
  PEE("");
}

