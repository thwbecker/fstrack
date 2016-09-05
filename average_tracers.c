#include "fstrack.h"
/*

  average tracer strain files with depth by doing a weighted sum (as a function
  of depth and/or age) of the individual L tensor components

  if flag use_weights_file is set, will read in file with depth
  dependent weights in format z [km, 0<z<2800, say] w
  
  if flag use_age is set, will weigh by age 
  
  for weighting, see routines below


  if use_isa flag is set, will assume infinite strain axes

  $Id: average_tracers.c,v 1.25 2006/08/22 20:02:34 becker Exp $

*/
int main(int argc, char **argv)
{
  int i,j,nzero,ntracer=0,otracer,ntlevels,start_search,found;
  my_boolean use_age=FALSE,use_weights_file=FALSE,use_isa=FALSE,skip;
  FILE *in;
  char filename[STRLEN],ostring[STRLEN];
  COMP_PRECISION w,strain_fac,tmp,dz,cage,evec[9],
    eval[3],picrit,depth,*tlevel;
  struct str *strain,tmp_strain[1];
  struct dw dweights;
  picrit = ISA_PI_CRIT_DEF;// Grain orientation lag cutoff factor
  cage = 50.;// critical age for age weighting, in Myr
  sprintf(ostring,"%s","dat");
  /*

    check for command line arguments

  */
  if(argc >= 2){
    if((strcmp(argv[1],"-h")==0)||(strcmp(argv[1],"-help")==0))
      argc=900;// bailout for help 
  }
  switch(argc){
  case 1:
    break;
  case 2:
    sscanf(argv[1],"%i",&i);
    if(i>0)use_weights_file=TRUE;
    break;
  case 3:
    sscanf(argv[1],"%i",&i);
    if(i>0)use_weights_file=TRUE;
    sscanf(argv[2],"%i",&i);
    if(i>0)use_age=TRUE;
    break;
  case 4:
    sscanf(argv[1],"%i",&i);
    if(i>0)use_weights_file=TRUE;
    sscanf(argv[2],"%i",&i);
    if(i>0)use_age=TRUE;
    sscanf(argv[3],"%i",&i);
    use_isa=(my_boolean)i;
    break;
  case 5:
    sscanf(argv[1],"%i",&i);
    if(i>0)use_weights_file=TRUE;
    sscanf(argv[2],"%i",&i);
    if(i>0)use_age=TRUE;
    sscanf(argv[3],"%i",&i);
    use_isa=(my_boolean)i;
    sscanf(argv[4],"%lf",&picrit);
    break;
  default:
    fprintf(stderr,"%s [use_weights_file, 0] [use_age, 0] [use_isa, 0] [picrit, %g]\n\n",
	    argv[0],picrit);
    fprintf(stderr,"averages strains in %s.ddd.dat files where the depths ddd are given by %s\n",
	    TSFILE,TDFILE);
    fprintf(stderr,"   a file that has a list of depths [km]\n");
    fprintf(stderr,"\nif use_weights_file is set to unity, will use file %s for reading in weights\n",
	    WSFILE);
    fprintf(stderr,"   this file is in format z[km] w\n");
    fprintf(stderr,"\nif use_age          is set to unity, will use age of tracer for averaging\n");
    fprintf(stderr,"   right now, the weighting function cuts off at %g Ma\n\n",cage);
    fprintf(stderr,"\nuse_isa          0: f.s finite strain tensor format\n");
    fprintf(stderr,"\n                 1: expects ISA format\n");
    fprintf(stderr,"\n                 2: expect  TI format\n");
    
    fprintf(stderr,"\npicrit              is the cutoff value, if the pi_par change rate is > picrit, will not use this tracer\n");
    exit(-1);
    break;
  }
  if(use_weights_file){
    read_weights_file(&dweights, argv);
    fprintf(stderr,"%s: read in %i weights from file %s\n",
	    argv[0],dweights.n,WSFILE);
  }
  if(use_isa && use_age){
    fprintf(stderr,"%s: error, can't use age and ISA or TI format\n",
	    argv[0]);
    exit(-1);
  }
  if(use_age)
    fprintf(stderr,"%s: using tracer ages for weights\n",argv[0]);
  if(use_isa==1)
    fprintf(stderr,"%s: expecting ISA format, critical pi: %g\n",
	    argv[0],picrit);
  else if(use_isa==2)
    fprintf(stderr,"%s: expecting TI format\n",argv[0]);
  //
  // read in depth levels from the tdepth.dat file, later we will
  // expect those depth levels to be reflected in the tracer.f.s.*.dat 
  // files
  //
  sprintf(filename,"%s.%s",TDFILE,ostring);
  in=myopen(filename,"r",argv[0]);
  ntlevels = 0;
  tlevel = (COMP_PRECISION *)malloc(sizeof(COMP_PRECISION));
  if(!tlevel)MEMERROR(argv[0]);
  while(fscanf(in,FLT_FORMAT,(tlevel+ntlevels))==1){
    ntlevels++;
    tlevel = (COMP_PRECISION *)
      realloc(tlevel,(ntlevels+1)*sizeof(COMP_PRECISION));
    if(!tlevel)MEMERROR(argv[0]);
  }
  fclose(in);
  if(ntlevels == 0){
    fprintf(stderr,"%s: error, no depth layers read\n",argv[0]);
    exit(-1);
  }else{
    fprintf(stderr,"%s: expecting %i different tracer layers as specified by %s\n",
	    argv[0],ntlevels,TDFILE);
  }
  /*

    read in strains from different depths

  */
  strain = (struct str *)malloc(sizeof(struct str));
  if(!strain)
    MEMERROR(argv[0]);
  /*

  start i - loop for files

  */
  otracer = 0;
  for(i=0;i < ntlevels;i++){	/* depth loop */
    sprintf(filename,"%s.%s.%g.%s","tracer",TSFILE,tlevel[i],ostring);
    in=myopen(filename,"r",argv[0]);
    ntracer = found = 0;
    while(((use_isa == 0)&&
	   //
	   // >>> lon lat depth Lrr Lrt Lrp Ltt Ltp Lpp time dz dx <<< is the normal format
	   // the last two items aren't used here!
	   //
	   (fscanf(in,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		   &tmp_strain->x[FSTRACK_X],&tmp_strain->x[FSTRACK_Y],&depth,
		   (tmp_strain->l+RR),(tmp_strain->l+RT),
		   (tmp_strain->l+RP),(tmp_strain->l+TT),
		   (tmp_strain->l+TP),(tmp_strain->l+PP),
		   &tmp_strain->age,eval,(eval+1)) == 12)) || 
	  //
	  // >>> lon lat depth ISAr ISAt ISAp use pi <<< is the ISA format
	  //
	  ((use_isa == 1)&&
	   (fscanf(in,"%lf %lf %lf %lf %lf %lf %lf %lf",
		  &tmp_strain->x[FSTRACK_X],&tmp_strain->x[FSTRACK_Y],&depth,
		   &evec[FSTRACK_E1*3+FSTRACK_R],&evec[FSTRACK_E1*3+FSTRACK_THETA],&evec[FSTRACK_E1*3+FSTRACK_PHI],
		   &tmp_strain->age,&tmp_strain->pipar) == 8)) ||
	  //
	  // TI format lon lat depth TIr TIt TIp azi_amp
	  //
	  ((use_isa == 2)&&
	   (fscanf(in,"%lf %lf %lf %lf %lf %lf %lf",
		   &tmp_strain->x[FSTRACK_X],&tmp_strain->x[FSTRACK_Y],&depth,
		   &evec[FSTRACK_E1*3+FSTRACK_R],&evec[FSTRACK_E1*3+FSTRACK_THETA],&evec[FSTRACK_E1*3+FSTRACK_PHI],
		   &strain_fac) == 7))
	  ){
      /* 
	 check if this tracer has really been advected to the target
	 depth, allow for a generous mismatch
      */
      if((fabs(depth) > EPS_COMP_PREC)&&(fabs((depth - tlevel[i])/depth) > 8e-2)){
	fprintf(stderr,"%s: depth mismatch: expected %g, but tracer %i has z: %g\n",
		argv[0],tlevel[i],ntracer+1,depth);
	exit(-1);
      }
      /* 

      prepare the new tracer to be added to the average

      */

      // obtain depth and age dependent general weight 
      w = averaging_weight(tmp_strain->age,cage,tlevel[i],
			   use_age,use_weights_file,
			   dweights,use_isa,tmp_strain->pipar,
			   picrit);
      //
      // average of largest stretching part
      //
      if(fabs(w) > EPS_PREC){
	// if w is zero, we don't have to calculate anything
	// for non_isa:
	//
	// calculate horizontal projection and direction of largest eigenvector
	// of L azimuth will be in degress
	//
	// for ISA: 
	// calculate the horizontal projection of the ISA eigenvector
	//
	if(!use_isa){// regular FSE ellipsoid
	  calc_eigensystem_sym(tmp_strain->l,eval,evec,TRUE);
	  //
	  // the reasoning here is: anisotropy scales with \zeta
	  // (\zeta = \ln(\lambda_1/\lambda_2) so that the strength of
	  // horizontal projection of e1 vector is log(e1/e2) *
	  // sqrt(1.0-e1r*e1r)
	  //
	  if((strain_fac = eval[FSTRACK_E1]/eval[FSTRACK_E2]) < 0){
	    fprintf(stderr,"%s: eigenvalue error, e1/e2<0, file: %s, ntracer: %i\n",
		    argv[0],filename,ntracer+1);
	    fprintf(stderr,"%s: is this an ISA format file, did you remember the switch?\n",
		    argv[0]);
	    exit(-1);
	  }
	  // scale with the log based on Ribe 1992
	  strain_fac = log(strain_fac);
	}else if(use_isa == 1){// ISA case
	  strain_fac = 1.0;// the zeta stretch ratio is meaningless 
	  tmp_strain->l[RR] = evec[FSTRACK_E1*3+FSTRACK_R];/* 
					      assign the r-component of the
					      largest eigenvector as radial
					      L_rr value for later
					      averaging
					   */
	}else if(use_isa == 2){	/* TI */
	  /* the scaling factor was already read in as amplitude */
	  tmp_strain->l[RR] = evec[FSTRACK_E1*3+FSTRACK_R];
	}
	//
	// scale the largest stretching with the horizontal projection of the 
	// e1-eigenvector, i.e. the horizontal component
	//
	strain_fac *= save_sqrt(1.0 - evec[FSTRACK_E1*3+FSTRACK_R]*evec[FSTRACK_E1*3+FSTRACK_R]);
	//
	// azimuth (CW from North) of horizontal projection, in degrees
	//
	tmp_strain->es[STR_ESI_AZI] = RAD2DEG(atan2(evec[FSTRACK_E1*3+FSTRACK_PHI],
						    -evec[FSTRACK_E1*3+FSTRACK_THETA]));
	//
	// we are only interested in orientations, not directions. hence, multiply
	// by two
	//
	tmp_strain->es[STR_ESI_AZI] *= 2.0;
	fix_deg_angle(&tmp_strain->es[STR_ESI_AZI]);
	//
	// x and y parts from azimuth
	calc_xy_from_azi_deg(tmp_strain->es[STR_ESI_AZI],
			     (tmp_strain->es+STR_ESI_AX),
			     (tmp_strain->es+STR_ESI_AY));
	//
	// linear strain scaling
	tmp_strain->es[STR_ESI_AX] *= strain_fac;
	tmp_strain->es[STR_ESI_AY] *= strain_fac;
	tmp_strain->es[STR_ESI_S_H] = strain_fac;
	//
	// radial part
	tmp_strain->es[STR_ESI_LRR] = tmp_strain->l[RR];
      }else{
	// for safety, initialize with zeroes in case w == 0
	for(j=0;j<9;j++)
	  tmp_strain->l[j] = 0.0;
	for(j=0;j<STR_N_ESI;j++)
	  tmp_strain->es[j] = 0.0;
	strain_fac = 0.0;
      }
      
      /* 

      here's where we add the strain to the sum

      */

      if(i == 0){
	//
	// init the strains tracers one by one since we are in the  first file
	//
	strain = (struct str *)realloc(strain,(otracer+1)*
				       sizeof(struct str));
	if(!strain)
	  MEMERROR("");
	/*

	initialize the ntracer - th  tracer

	*/
	// lon-lat location of tracer
	strain[otracer].x[FSTRACK_X] = tmp_strain->x[FSTRACK_X];
	strain[otracer].x[FSTRACK_Y] = tmp_strain->x[FSTRACK_Y];
	if(!use_isa)
	  for(j=0;j<9;j++)// l strain components
	    strain[otracer].l[j] = tmp_strain->l[j] * w;
	for(j=0;j < STR_N_ESI;j++)// horizontal eigensystem, see above
	  strain[otracer].es[j] = tmp_strain->es[j] * w;
	strain[otracer].wss = strain_fac * w; // strain dependent weight
	strain[otracer].ws  = w;// depth and age dependent weight
	strain[otracer].dazi = 0.0;// absolute difference in azimuths
	strain[otracer].oldazi = tmp_strain->es[STR_ESI_AZI];
	/* 

	the first file determines the maximum number of tracers

	*/
	otracer++;		/* original tracer counter */
	/* 

	end first depth level
	
	*/
      }else{// add the strains to existing strain tracer
	/* find the proper tracer to add things to */
	start_search = ntracer;
	while((!colocated((strain+ntracer),tmp_strain,depth,20))&&(ntracer < otracer))
	  ntracer++;
	if(ntracer == otracer){
	  fprintf(stderr,"%s: WARNING: no match found for %g, %g in file %s\n",
		  argv[0],tmp_strain->x[FSTRACK_X],tmp_strain->x[FSTRACK_Y],filename);
	  skip = TRUE;
	  ntracer = start_search; /* reset the tracer counter */
	}else{
	  skip = FALSE;
	  found++;
	}
	if(!skip){
	  /* add to the depth average */
	  if(!use_isa)
	    for(j=0;j<9;j++)
	      strain[ntracer].l[j]  += tmp_strain->l[j] * w;
	  for(j=0;j < STR_N_ESI;j++)
	    strain[ntracer].es[j] += tmp_strain->es[j] * w;
	  strain[ntracer].wss += strain_fac * w;
	  strain[ntracer].ws  += w;
	  // calculate angular 2theta deviation
	  tmp = fabs(tmp_strain->es[STR_ESI_AZI]-
		     strain[ntracer].oldazi);
	  if(tmp > 180)
	    tmp = 360 - tmp;
	  strain[ntracer].dazi += tmp;
	  // 2 theta, for the summing
	  strain[ntracer].oldazi = tmp_strain->es[STR_ESI_AZI];
	}
      }

    } /* end of tracer in file loop */
    fclose(in);
    if(i != 0){
      if(found != otracer){
	fprintf(stderr,"%s: WARNING: nr of point mismatch: %i (now) vs %i (before), file: %s\n",
		argv[0],found,otracer,
		filename);
      }
    }else{
      if(!otracer){
	fprintf(stderr,"%s: error, didn't read any tracers from %s\n",
		argv[0],filename);
	exit(-1);
      }
    }
    fprintf(stderr,"%s: read %i strain tracers from %s, depth dependent weight: %g\n",
	    argv[0],ntracer,filename,
	    averaging_weight(tmp_strain->age,cage,tlevel[i],
			     use_age,use_weights_file,dweights,
			     FALSE,tmp_strain->pipar,picrit));
  }
  //
  // average
  //
  fprintf(stderr,"%s: averaging %i tracers\n",argv[0], otracer);
  for(nzero=i=0;i < otracer;i++){
    if((fabs(strain[i].ws) < EPS_PREC)||(fabs(strain[i].wss) < EPS_PREC))
      nzero++;// count the tracers without real entries
    /*
      
    in case no valid values were read in (w==0), this will 
    result in a bunch of NaNs, but that's OK

    */
    if(!use_isa)
      for(j=0;j<9;j++) 
	strain[i].l[j]  /= strain[i].ws; 
    /*
      
    replace the azimuth average by the one estimated 
    from ax and ay
    
    */
    strain[i].es[STR_ESI_AX] /= strain[i].wss; 
    strain[i].es[STR_ESI_AY] /= strain[i].wss; 
    calc_azi_from_axay(strain[i].es[STR_ESI_AX], 
		       strain[i].es[STR_ESI_AY],  
		       &strain[i].es[STR_ESI_AZI]);  
    strain[i].es[STR_ESI_AZI] /= 2.0;/* correct angle  */
    fix_deg_angle((strain[i].es+STR_ESI_AZI));  
    strain[i].es[STR_ESI_LRR] /= strain[i].ws;      /* L_rr part  */ 
    strain[i].es[STR_ESI_S_H] = strain[i].wss / strain[i].ws;  
    strain[i].dazi /= 2.0;  
  }
  if(nzero)
    fprintf(stderr,"%s: WARNING: %i out of %i tracers weren't assigned any averaged values\n",
	    argv[0],nzero,otracer);
  if(!use_isa){
    //
    // output of average tensor components for regular averaging, not ISA nor TI
    //
    sprintf(filename,"%s.%s.avg.%s","tracer",TSFILE,ostring);
    in=myopen(filename,"w",argv[0]);
    fprintf(stderr,"%s: writing averaged tensor components to %s\n",argv[0],filename);
    for(i=0;i < otracer;i++)
      fprintf(in,"%11.4e %11.4e %11.4e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %10g %11g %11g\n",
	      strain[i].x[0],strain[i].x[1],0.0,
	      strain[i].l[RR],strain[i].l[RT],strain[i].l[RP],
	      strain[i].l[TT],strain[i].l[TP],strain[i].l[PP],
	      0.0,0.0,0.0);
    fclose(in);
  }
  //
  // output of average eigensystem
  //
  sprintf(filename,"%s.%s.es.avg.%s","tracer",TSFILE,ostring);
  in=myopen(filename,"w",argv[0]);
  fprintf(stderr,"%s: writing averaged horizontal eigensystem to %s\n",argv[0],filename);
  for(i=0;i < otracer;i++)// format is lon lat azi eh*log(e1/e2) Lrr
    fprintf(in,"%11.4e %11.4e %11g %13.6e %13.6e\n",
	    strain[i].x[0],strain[i].x[1],
	    strain[i].es[STR_ESI_AZI],strain[i].es[STR_ESI_S_H],
	    strain[i].es[STR_ESI_LRR]);
  fclose(in);
  //
  // output of total absolute fast angular difference per length interval
  //
  dz = tlevel[ntlevels-1] - tlevel[0];
  fprintf(stderr,"%s: total layer thickness determined as %g km\n",argv[0],dz);
  sprintf(filename,"%s.%s.dazi.%s","tracer",TSFILE,ostring);
  in=myopen(filename,"w",argv[0]);
  fprintf(stderr,"%s: writing total absolute fast angular differences to %s\n",
	  argv[0],filename);
  for(i=0;i < otracer;i++)
    fprintf(in,"%11.4e %11.4e %11.4e \n",
	    strain[i].x[0],strain[i].x[1],strain[i].dazi/dz);
  fclose(in);
  return 0;
}
/*

  weights routine

  age:       age of tracer
  cage:      cutoff age for weights
  use_age:   use age weighting
  z:         depth
  use_depth: use depth weighting 
  dwegihts:  depth weight structure
  use_isa: input is in the ISA format, in this case "age" is the "use" flag
  pi: for ISA, pi parameter
  picrit: the critical pi change level for bailout (<< 1)

*/
COMP_PRECISION averaging_weight(COMP_PRECISION age, COMP_PRECISION cage, 
				COMP_PRECISION z,int use_age,int use_depth,
				struct dw dweights, my_boolean use_isa,
				COMP_PRECISION pi,COMP_PRECISION picrit)
{
  COMP_PRECISION w;
  w = 1.0;
  if(use_age){// weight by the age of the tracer
    if(use_isa){
      fprintf(stderr,"average_tracers: weight: error: ISA/TI and use_age doesn't make sense\n");
      exit(-1);
    }
    w *= (age <= cage)?(1.0-(age/cage)):(0.0);
  }
  if(use_depth){
    // weight by the depth dependence as given by the kernels
    w *= interpolate_dw(z,dweights);
  }
  if(use_isa == 1){
    // this is a ISA tracer, in this case "age" is the "use" flag
    if(fabs(age) < EPS_PREC)// don't use this tracer
      w = 0.0;// set the weight to zero
    else{
      // we can probably use this tracer
      if(fabs(age) > EPS_PREC){
	fprintf(stderr,"average_tracers: weight: error: use_isa is set, but use flag (%g) is neither zero nor unity\n",
		age);
	exit(-1);
      }
      if(pi > picrit)// too large a change level for gamma
	w = 0.0;
    }
  }
  return w;
}
//
//
// linearly interpolate weighting function to depth z
//
// WARNING: the routine expects that dweights doesn't change during
// program run
//
COMP_PRECISION interpolate_dw(COMP_PRECISION z, 
			      struct dw dweights)
{
  static my_boolean init=FALSE;
  static int n1;
  int i,j;
  COMP_PRECISION fac1,fac2;
  if(!init){// check if  depth array is sorted
    for(i=1;i<dweights.n;i++)
      if(dweights.z[i]<=dweights.z[i-1]){
	fprintf(stderr,"interpolate_dw: depth levels for weights not in ascending order\n");
	for(i=0;i<dweights.n;i++)
	  fprintf(stderr,"i: %3i z: %11g w: %11g\n",i,dweights.z[i],dweights.w[i]);
	exit(-1);
      }
    if(dweights.n<2){
      fprintf(stderr,"interpolate_dw: need two datapoints for interpolation\n");
      exit(-1);
    }
    n1 = dweights.n-1;
    init=TRUE;
  }
  get_lin_int_weights(z,dweights.z,n1,&i,&j,&fac1,&fac2);
  return fac1 * dweights.w[i] + fac2 * dweights.w[j];
}
/* 

  read in depth dependent weights, e.g. normalized kernels
  
*/
void read_weights_file(struct dw *dweights, char **argv)
{
  FILE *in;
  size_t scnt;
  in=myopen(WSFILE,"r",argv[0]);
  dweights->n = 0;
  scnt=sizeof(COMP_PRECISION);
  dweights->z=(COMP_PRECISION *)malloc(scnt);
  dweights->w=(COMP_PRECISION *)malloc(scnt);
  while(fscanf(in,"%lf %lf",(dweights->z+dweights->n),
	       (dweights->w+dweights->n)) == 2){
    if(dweights->n > 0){// check if two depths are same, shift then
      if(fabs(dweights->z[dweights->n] - dweights->z[dweights->n-1])<EPS_PREC){
	dweights->z[dweights->n] += 1.0e-6;
	dweights->z[dweights->n-1] -= 1.0e-6;
	fprintf(stderr,"%s: splitting weights %4i and %4i to depths %11e and %11e\n",
		argv[0],dweights->n-1,dweights->n,
		dweights->z[dweights->n-1],dweights->z[dweights->n]);
      }else if(dweights->z[dweights->n]<dweights->z[dweights->n-1]){//ascending?
	fprintf(stderr,"%s: can not deal with descreasing depth ordering\n",
		argv[0]);
	fprintf(stderr,"%s: weights %i and %i are at depths %g and %g\n",
		argv[0],dweights->n-1,dweights->n,dweights->z[dweights->n-1],
		dweights->z[dweights->n]);
	exit(-1);
      }
    }
    dweights->n++;scnt+=sizeof(COMP_PRECISION);
    dweights->z=(COMP_PRECISION *)realloc(dweights->z,scnt);
    dweights->w=(COMP_PRECISION *)realloc(dweights->w,scnt);
    if((!dweights->z) || (!dweights->w))MEMERROR(argv[0]);
  }
  fclose(in);
}
/* 
   expects lon/lat and depth in km
 */
my_boolean colocated(struct str *a,struct str *b, COMP_PRECISION depth,
		     COMP_PRECISION eps_km)
{
  COMP_PRECISION ap[3],bp[3];
  ap[FSTRACK_R] = bp[FSTRACK_R] =  ND_RADIUS(fabs(depth)); /*  */
  ap[FSTRACK_PHI] = LONGITUDE2PHI(a->x[FSTRACK_X]);bp[FSTRACK_PHI]=LONGITUDE2PHI(b->x[FSTRACK_X]);
  ap[FSTRACK_THETA] = LATITUDE2THETA(a->x[FSTRACK_Y]); bp[FSTRACK_THETA] = LATITUDE2THETA(b->x[FSTRACK_Y]);
  if(dist_on_sphere(ap,bp) <= eps_km)
    return TRUE;
  else
    return FALSE;
}
