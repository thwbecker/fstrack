#include "fstrack_flow.h"
//
// init the tracers
//
// $Id: init.c,v 1.19 2016/09/05 04:44:47 becker Exp $
//

/*

  distribute tracers over sphere

  modes can be DIST_EVEN_DX, DIST_AREA_EVEN, or SPOTTED_FROM_FILE

  tracer parameters are defined in init_parameter

  different depth levels get stacked up with model->ntd[i] tracers in each 
  level

*/

void initialize_tracers(struct mod *model)
{
  int i,j,k,l,j3,n,m,tcnt,ntperlayer,nt3;
  my_boolean bailout=FALSE;
  COMP_PRECISION x[3],*tloc,dx[3],dphi0,phimax,thetamax,lon,lat,z,*attr,add_strain[9];
  FILE *in;
  static my_boolean init=FALSE;
  tloc=NULL;phimax=FLT_MIN;attr=NULL;
  if(init){
    fprintf(stderr,"initialize_tracers: WARNING: called more than once, untested\n");
    model->ntracer=0;
    // clear the tracers
    free(model->tracer);
    model->tracer = (struct trc *)calloc(1,sizeof(struct trc));
    add_state((model->tracer),0.0,MDP);
  }
  // initialize the tracer depth levels
  init_tracer_depths(model);
  if(model->tinitmode == SPOTTED_LAT_FROM_FILE){
    //
    // read lateral tracer distribution from file and redefine the number of 
    // tracers this way
    //
    in=myopen(model->tlatfilename,"r","distribute_tracers");
    model->ntracer=0;
    nt3=0;
    my_vecalloc(&tloc,3,"distribute_tracers: 0");
    while(fscanf(in,TWO_FLT_FORMAT,&lon,&lat)==2){
      tloc[nt3 + FSTRACK_PHI]=  LONGITUDE2PHI(lon);
      tloc[nt3 + FSTRACK_THETA]=LATITUDE2THETA(lat); 
      model->ntracer++;nt3 += 3;
      my_vecrealloc(&tloc,(nt3+3),"distribute_tracers: 1");
    }
    fclose(in);
    fprintf(stderr,"distribute_tracers: read %i lateral tracer locations from %s\n",
	    model->ntracer,model->tlatfilename);
    model->ntracer *= model->ntlevels;// make room for other layers with depth
  }else if(model->tinitmode == SPOTTED_3D_FROM_FILE){
    // read 3-D tracer distribution from file and redefine the number of 
    // tracers this way
    in=myopen(model->tlatfilename,"r","distribute_tracers");
    model->ntracer=0;
    nt3=0;
    my_vecalloc(&tloc,3,"distribute_tracers: 2");
    while(fscanf(in,THREE_FLT_FORMAT,&lon,&lat,&z)==3){
      tloc[nt3 + FSTRACK_PHI]=  LONGITUDE2PHI(lon);
      tloc[nt3 + FSTRACK_THETA]=LATITUDE2THETA(lat); 
      tloc[nt3 + FSTRACK_R]=ND_RADIUS(z);
      check_tracer_location((tloc+nt3));
      model->ntracer++;nt3 += 3;
      my_vecrealloc(&tloc,(nt3+3),"distribute_tracers: 3");
    }
    fclose(in);
    fprintf(stderr,"distribute_tracers: read %i 3-D tracer locations from %s\n",
	    model->ntracer,model->tlatfilename);
  }else if(model->tinitmode == SPOTTED_3D_WITH_ATTR){
    // read 3-D tracer distribution with attributes
    in=myopen(model->tlatfilename,"r","distribute_tracers");
    model->ntracer=0;nt3=0;
    my_vecalloc(&tloc,3,"distribute_tracers: 4");
    my_vecalloc(&attr,NR_T_ATTR,"distribute_tracers: 5");
    while(fscanf(in,THREE_FLT_FORMAT,&lon,&lat,&z)==3){// read in locations
      for(l=0;l<NR_T_ATTR;l++)// read in attribute(s)
	if(fscanf(in,FLT_FORMAT,(attr+model->ntracer*NR_T_ATTR+l))!=1){
	  PE("init: read error");exit(-1);}
      tloc[nt3 + FSTRACK_PHI]=  LONGITUDE2PHI(lon);
      tloc[nt3 + FSTRACK_THETA]=LATITUDE2THETA(lat); 
      tloc[nt3 + FSTRACK_R]=ND_RADIUS(z);
      check_tracer_location((tloc+nt3));
      model->ntracer++;nt3 += 3;
      my_vecrealloc(&tloc,(nt3+3),"distribute_tracers: 6");
      my_vecrealloc(&attr,(model->ntracer+1)*NR_T_ATTR,"distribute_tracers: 7");
    }
    fclose(in);
    fprintf(stderr,"distribute_tracers: read %i 3-D tracers and attributes from %s\n",
	    model->ntracer,model->tlatfilename);
  }
  //
  // allocate ntracer tracers and initialize with zeroes
  //
  allocate_and_clear_tracers(model,model->ntracer);
  ntperlayer = model->ntracer/model->ntlevels;
  for(tcnt=i=0;i<model->ntlevels;i++){// loop through depths
    //
    // start assigning locations
    //
    model->ntd[i]=0;
    switch(model->tinitmode){
    case SPOTTED_LAT_FROM_FILE:
      for(j=j3=0;j<ntperlayer;j++,j3+=3){
	model->tracer[tcnt].state[0].x[FSTRACK_PHI] =   tloc[j3+FSTRACK_PHI];
	model->tracer[tcnt].state[0].x[FSTRACK_THETA] = tloc[j3+FSTRACK_THETA];
	model->tracer[tcnt].state[0].x[FSTRACK_R] = model->tlevel[i];
	model->ntd[i]++;tcnt++;
      }
      fprintf(stderr,"init_tracers: generated %i from file at z: %g\n",model->ntd[i],
	      ZDEPTH(model->tlevel[i]));
      break;
    case SPOTTED_3D_FROM_FILE:
    case SPOTTED_3D_WITH_ATTR:
      for(j=j3=0;j<ntperlayer;j++,j3+=3){
	if(model->tinitmode == SPOTTED_3D_WITH_ATTR)// attribute(s)
	  for(l=0;l<NR_T_ATTR;l++)
	    model->tracer[tcnt].attr[l] = attr[j*NR_T_ATTR+l];
	for(k=0;k<3;k++)// locations
	  model->tracer[tcnt].state[0].x[k] =   tloc[j3+k];
	if((model->amode == BACKWARD_STRAIN)&&
	   (model->tracer[tcnt].state[0].x[FSTRACK_R] <= model->sdepth)){
	  fprintf(stderr,"init_tracer_depths: error: mode is backward strain with limit depth %g\n",
		  ZDEPTH(model->sdepth));
	  fprintf(stderr,"init_tracer_depths: but tracer depth level init at %g attempted\n",
		  ZDEPTH(model->tracer[tcnt].state[0].x[FSTRACK_R]));
	  exit(-1);
	}
	model->ntd[i]++;tcnt++;
      }
      break;
    case DIST_EVEN_AREA:
      /* 
	 evenly (1/sin(theta)) distribute tracers of the surface of a sphere 
      */
      m=(int)(sqrt((COMP_PRECISION)ntperlayer)*1.75);
      dphi0 = TWOPI / (COMP_PRECISION)m;
      dx[FSTRACK_THETA]=dphi0;
      x[FSTRACK_THETA]= dx[FSTRACK_THETA];// start theta loop
      thetamax=PI - dx[FSTRACK_THETA] + 1.0e-3;
      while((x[FSTRACK_THETA] <= thetamax)&&(!bailout)){
	dx[FSTRACK_PHI] = dphi0/(COMP_PRECISION)sin(x[FSTRACK_THETA]);// phi loop
	x[FSTRACK_PHI]=dx[FSTRACK_PHI]/2.0;
	phimax=TWOPI - dx[FSTRACK_PHI]/2.0 + 1.0e-3;
	while((x[FSTRACK_PHI]<=phimax)&&(!bailout)){
	  model->tracer[tcnt].state[0].x[FSTRACK_PHI] = x[FSTRACK_PHI];
	  model->tracer[tcnt].state[0].x[FSTRACK_THETA] = x[FSTRACK_THETA];
	  model->tracer[tcnt].state[0].x[FSTRACK_R] = model->tlevel[i];
	  model->ntd[i]++;tcnt++;
	  if(model->ntd[i] == ntperlayer)
	    bailout=TRUE;
	  x[FSTRACK_PHI] += dx[FSTRACK_PHI];
	}
	x[FSTRACK_THETA] += dx[FSTRACK_THETA];
      }
      if(bailout){
	fprintf(stderr,"WARNING: bailout at theta/phi: %g%%/%g%%\n",
		x[FSTRACK_THETA]/thetamax*100.0,x[FSTRACK_PHI]/phimax*100.0);
      }
      fprintf(stderr,"init_tracers: generated %i even area type at z: %g\n",model->ntd[i],
	      ZDEPTH(model->tlevel[i]));
      break;
    case DIST_EVEN_DX:
      /*
	grid with fixed dtheta dphi distribution from dtheta <= theta <= pi - dtheta
	and 0 <= phi <= 2pi -dphi
      */
      n=(int)(sqrt(((COMP_PRECISION)ntperlayer)/2.0));
      m=2*n;
      dx[FSTRACK_PHI] = TWOPI / ((COMP_PRECISION)m);
      dx[FSTRACK_THETA]= PI / ((COMP_PRECISION)n);
      x[FSTRACK_THETA]=dx[FSTRACK_THETA];// start theta loop
      thetamax=PI-dx[FSTRACK_THETA]+1.0e-3;
      phimax=TWOPI-dx[FSTRACK_PHI]/2.0+1.0e-3;
      while((x[FSTRACK_THETA]<=thetamax)&&(!bailout)){
	x[FSTRACK_PHI]=dx[FSTRACK_PHI]/2.0;
	while((x[FSTRACK_PHI]<=phimax)&&(!bailout)){
	  model->tracer[tcnt].state[0].x[FSTRACK_PHI] = x[FSTRACK_PHI];
	  model->tracer[tcnt].state[0].x[FSTRACK_THETA] = x[FSTRACK_THETA];
	  model->tracer[tcnt].state[0].x[FSTRACK_R] = model->tlevel[i];
	  model->ntd[i]++;tcnt++;
	  if(model->ntd[i] == ntperlayer)
	    bailout=TRUE;
	  x[FSTRACK_PHI] += dx[FSTRACK_PHI];
	}
	x[FSTRACK_THETA] += dx[FSTRACK_THETA];
      }
      if(bailout){
	fprintf(stderr,"WARNING: bailout at theta/phi: %g%%/%g%%\n",
		x[FSTRACK_THETA]/thetamax*100.0,x[FSTRACK_PHI]/phimax*100.0);
      }
      fprintf(stderr,"init_tracers: generated %i fixed dphi dtheta tracers at z: %g\n",
	      model->ntd[i],ZDEPTH(model->tlevel[i]));
      break;
    default:
      fprintf(stderr,"init_tracers: mode %i undefined\n",model->tinitmode);
      exit(-1);
      break;
    }
  }
  fprintf(stderr,"init_tracers: allocated %i tracers in total\n",tcnt);
  //
  // assign real number of tracers and shrink
  model->ntracer=tcnt;
  model->tracer=realloc(model->tracer,sizeof(struct trc)*model->ntracer);
  if((model->tinitmode == SPOTTED_LAT_FROM_FILE)||// free tloc if reading from file
     (model->tinitmode == SPOTTED_3D_FROM_FILE))
    free(tloc);
  if(model->tinitmode == SPOTTED_3D_WITH_ATTR){
    free(tloc);free(attr);
  }
  //
  // check locations 
  for(i=0;i<model->ntracer;i++){
    check_phys_lim_tracer(model->tracer[i].state[0].x,model->tracer[i].state[0].x);
    if(close_to_pole(model->tracer[i].state[0].x)){
      fprintf(stderr,"init_tracers: tracer %i too close to poles: lon/lat: %g/%g\n",
	      i,PHI2LONGITUDE(x[FSTRACK_PHI]),
	      THETA2LATITUDE(x[FSTRACK_THETA]));
      exit(-1);
    }
  }
  /* 
     possibly add initial strain to tracers
  */
  if(model->use_initial_strain){
    fprintf(stderr,"init_tracers: WARNING: adding strain: (%g,%g,%g,%g,%g,%g,%g,%g,%g) (spherical system)\n",
	    model->initial_strain[0],model->initial_strain[1],model->initial_strain[2],
	    model->initial_strain[3],model->initial_strain[4],model->initial_strain[5],
	    model->initial_strain[6],model->initial_strain[7],model->initial_strain[8]);
    for(i=0;i<model->ntracer;i++){
      /* convert to internal, Cartesian system */
      polar_to_cart_mat_at_r(model->initial_strain,add_strain,model->tracer[i].state[0].x);
      for(j=0;j<9;j++)
	model->tracer[i].state[0].x[3+j] += add_strain[j];
    }
  }
  init = TRUE;
}


//
// allocate ntracer tracers with one state and clear
//
void allocate_and_clear_tracers(struct mod *model,int ntracer)
{
  int i;
  my_boolean init=FALSE;
  if(ntracer < 0){
    fprintf(stderr,"allocate_tracers: error: ntracer: %i\n",ntracer);
    exit(-1);
  }
  if(!init){// first initialization
    model->tracer = (struct trc *)calloc(ntracer,sizeof(struct trc));
    init=TRUE;
  }else{// reallocation
    model->tracer = (struct trc *)realloc(model->tracer,ntracer*sizeof(struct trc));
  }
  if(!model->tracer)
    MEMERROR("allocate_tracers");
  for(i=0;i < model->ntracer;i++)
    clear_tracer((model->tracer+i));
  model->ntracer = ntracer;
  for(i=0;i<model->ntracer;i++)
    add_state((model->tracer+i),0.0,MDP);// generate a state, initialized with t=0
 
}
/* 
   clear a tracer 
*/
void clear_tracer(struct trc *tracer)
{
  int i;
  tracer->discarded = FALSE;
  if(tracer->nstate){
    for(i=0;i<tracer->nstate;i++){
      if(tracer->state[i].npole){
	free(tracer->state[i].pdens);
	tracer->state[i].npole = 0;
      }
    }
    free(tracer->state);
  }
  tracer->nstate = 0;		/* reset state counter */
  for(i=0;i < NR_T_ATTR;i++)
    tracer->attr[i] = 0.0;
}

void init_tracer_depths(struct mod *model)
{
  FILE *in;
  int i;
  static my_boolean init=FALSE;
  if(init){
    free(model->tlevel);
    free(model->ntd);
  }
  model->tlevel = (COMP_PRECISION *)malloc(sizeof(COMP_PRECISION));
  model->ntd =    (unsigned int *)calloc(1,sizeof(unsigned int));
  if(!model->tlevel || !model->ntd)
    MEMERROR("init_tracer_depth");
  if((model->tinitmode == SPOTTED_3D_FROM_FILE)||
     (model->tinitmode == SPOTTED_3D_WITH_ATTR)){// depths will be read from file, too
    model->ntlevels=1;// fake depth level
    model->tlevel[0]=ND_RADIUS(-1);
  }else{
    if(model->tdepth_from_file){
      in = myopen(model->tdepthfilename,"r","init_tracer_depths");
      // read in the depth levels
      i=0;
      while(fscanf(in,FLT_FORMAT,&model->tlevel[i])==1){
	model->tlevel[i]=ND_RADIUS(model->tlevel[i]);
	if((model->amode == BACKWARD_STRAIN)&&
	   (model->tlevel[i] <= model->sdepth)){
	  fprintf(stderr,"init_tracer_depths: error: mode is backward strain with limit depth %g\n",
		  ZDEPTH(model->sdepth));
	  fprintf(stderr,"init_tracer_depths: but tracer depth level init at %g attempted\n",
		  ZDEPTH(model->tlevel[i]));
	  exit(-1);
	}
	model->ntd[i]=0;
	i++;
	model->tlevel=(COMP_PRECISION *)realloc(model->tlevel,
						sizeof(COMP_PRECISION)*(i+1));
	model->ntd=(unsigned int *)realloc(model->ntd,sizeof(unsigned int)*(i+1));
	if(!model->tlevel || !model->ntd)
	  MEMERROR("init_tracer_depth");
      }
      model->ntlevels=i;
      if(model->ntlevels==0){
	fprintf(stderr,"init_tracer_depths: did not read any depth levels from \"%s\n",
		model->tdepthfilename);
	exit(-1);
      }else{
	fprintf(stderr,"init_tracer_depths: read %i depth levels from \"%s\"\n",
		model->ntlevels,model->tdepthfilename);
      }
      fclose(in);
    }else{
      // assign only one standard depth level
      model->ntlevels=1;
      model->tlevel[0]=ND_RADIUS(model->init_tracer_depth);
      if((model->amode == BACKWARD_STRAIN)&&
	 (model->tlevel[0] <= model->sdepth)){
	fprintf(stderr,"init_tracer_depths: error: mode is backward strain with limit depth %g\n",
		  ZDEPTH(model->sdepth));
	fprintf(stderr,"init_tracer_depths: but tracer depth level init at %g attempted\n",
		ZDEPTH(model->tlevel[0]));
	exit(-1);
      }
    }
  }
  init = TRUE;
}

/*

  check if coordinate is close to South or North pole
  input is a 3-D vector in polar coordinates
*/

#define CPOLE_EPS 1.0e-4
my_boolean close_to_pole(COMP_PRECISION *x)
{
  static COMP_PRECISION sp = PI - CPOLE_EPS;
  if((x[FSTRACK_THETA] <  CPOLE_EPS) || (x[FSTRACK_THETA] > sp)){
    return TRUE;
  }else 
    return FALSE;

}
void check_tracer_location(COMP_PRECISION *x)
{
  if(x[FSTRACK_R]<0 || x[FSTRACK_R]>1.0){
    fprintf(stderr,"check_tracer_location: radial tracer component (r: %g z: %g) out of bounds\n",
	    x[FSTRACK_R],ZDEPTH(x[FSTRACK_R]));
    exit(-1);
  }
}
