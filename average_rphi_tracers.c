
#include "fstrack.h"
/*

usage:

average_rphi depth.dat nr_period [restrict_ani,0] [kernel_dir]


--- averages the 2phi and 4phi terms for Rayleigh waves ---



if nr_period > 0

RPHI MODE: uses precomputed kernel averages from rphi file


average rphi type tracer files which are in format

lon lat z L N T_1 azi_2p^T1 amp_2p^T1 azi_4p^T1 amp_4p^T1 T2 azi_2p^T2 amp_2p^T2 azi_4p^T2 amp_4p^T2...

L = \rho v_sv^2
N = \rho v_sh^2

azi is CW from North in degrees


if nr_period < 0

SAV MODE: reads sav files, and computes average for several periods



DOES NOT DEAL WITH LOVE WAVES



*/

struct rphi{
  COMP_PRECISION lon,lat;
  COMP_PRECISION T[N_SW_SENS],ax[N_SW_SENS][2][2];
};


int main(int argc, char **argv)
{
  FILE *in,*out;
  int ntlevels,nperiod,ntracer=0,i,j,k,max_ntracer=0,restrict_ani=0,rphi_mode = 1;
  COMP_PRECISION *tlevel,T[N_SW_SENS],l_fac,n_fac,scale,sav[36],loc_azi[2],
    amp[N_SW_SENS][2],lon,lat,depth,xp[3],loc_afactor[4],loc_bfactor[4],
    kernel[5],ba[2],hf[2],gl[2],ca[2],cn[2],a_fac,c_fac,f_fac;
  char filename[STRLEN];
  struct rphi *tracer=NULL;
  struct mod *model;
  /* for sensitityv */
  model= (struct mod *)calloc(1,sizeof(struct mod));
  
  if(argc < 3){
    fprintf(stderr,"%s depth.dat nr_periods [restrict_ani,0] [kernel_dir, %s]\n",
	    argv[0],SW_SENS_FILE);
    fprintf(stderr,"nr_period < 0 switches to sav mode, else expects rphi files\n");
    exit(-1);
  }
  //
  // read in depth levels from a file, later we will
  // expect those depth levels to be reflected in the tracer.f.s.*.dat 
  // files
  //
  sscanf(argv[2],"%i",&nperiod);
  if(argc > 3)
    sscanf(argv[3],"%i",&restrict_ani);
  if(nperiod < 0){
    /* sav mode  */
    rphi_mode = 0;
    /* directory for sensitivity kernels */
    if(argc > 4)
      sprintf(model->sw_sens_file,"%s",argv[4]);
    else
      sprintf(model->sw_sens_file,"%s",SW_SENS_FILE);
    /* initialize the sensitivity kernels */
    read_sw_sens(model);

    nperiod = N_SW_SENS;

  }else{

    if((nperiod > N_SW_SENS)||(nperiod < 1)){
      fprintf(stderr,"%s: nperiod (%i) not valid, max is %i\n",
	      argv[0],nperiod,N_SW_SENS);
      exit(-1);
    }
  }
  
  /* read depth file */
  in=myopen(argv[1],"r",argv[0]);
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
  }
  if(rphi_mode)
    fprintf(stderr,"%s: expecting %i depth layers and %i periods in rphi files, restrict: %i\n",
	    argv[0],ntlevels,nperiod,restrict_ani);
  else
    fprintf(stderr,"%s: expecting %i depth layers in sav files, interpolatin for several periods, restrict: %i\n",
	    argv[0],ntlevels,restrict_ani);
  for(i=0;i < ntlevels;i++){
    /* 
       
    start depth loop

    */
    if(rphi_mode)
      sprintf(filename,"tracer.rphi.%g.dat",tlevel[i]);
    else
      sprintf(filename,"tracer.sav.%g.dat",tlevel[i]);
    fprintf(stderr,"%s: reading from %s\n",argv[0],filename);
    in = myopen(filename,"r",argv[0]);
    ntracer = 0;
    while(fscanf(in,THREE_FLT_FORMAT,&lon,&lat,&depth) == 3){
      if(depth < 0){
	fprintf(stderr,"%s: error: expected positive depths, read %g\n",
		argv[0],depth);
	exit(-1);
      }
      if(fabs(depth -tlevel[i]) > MAX(0.15,fabs(0.02*tlevel[i]))){
	fprintf(stderr,"%s: depth mismatch: read %g, expected: %g in file: %s, line: %i\n",
		argv[0],depth,tlevel[i],filename,ntracer+1);
	exit(-1);
      }
      if(lon >= 360)
	lon -= 360.0;
      /* 
	 check if we need to restrict the anisotropy 
      */
      scale = restrict_ani_factor(restrict_ani,depth);

      /* 
	 
      decide which mode to use

      */
      if(rphi_mode){
	/* 
	   
	REGULAR MODE, READ IN PRECOMPUTED KERNEL AVERAGES (WE NOW PREFER TO RECOMPUTE, IN CASE)
	
	
	*/
	/* 
	   read in the l, n factors
	   l=\rho v_sv^2
	   n=\rho v_sh^2
	   
	*/
	if(fscanf(in,TWO_FLT_FORMAT,&l_fac,&n_fac) != 2){
	  fprintf(stderr,"%s: read error l/n factors\n",argv[0]);
	  exit(-1);
	}
	for(j=0;j < nperiod;j++){
	  /* 
	     
	  read in the local fast azimuths and strengths
	  
	  */
	  if(fscanf(in,"%lf %lf %lf  %lf %lf",(T+j),
		    &loc_azi[0],&amp[j][0],&loc_azi[1],&amp[j][1])!=5){
	    fprintf(stderr,"%s: read error period data\n",argv[0]);
	    exit(-1);
	  }
	  if((amp[j][0] < 0 )||(amp[j][1] < 0)){
	    fprintf(stderr,"%s: tracer %i depth %g amplitude error: %g %g\n",
		    argv[0],ntracer,tlevel[i],amp[j][0],amp[j][1]);
	    exit(-1);
	  }
	}
	if(i == 0){	
	  /* 
	     first time around 
	  */
	  tracer = (struct rphi *)realloc(tracer,sizeof(struct rphi)*(ntracer+1));
	  if(!tracer)
	    MEMERROR("");
	  tracer[ntracer].lon = lon;tracer[ntracer].lat=lat;
	  for(j=0;j < nperiod;j++){
	    tracer[ntracer].T[j] = T[j];
	    for(k=0;k < 2;k++){
	      /* convertt azimuth */
	      loc_azi[k] = DEG2RAD(((double)(2*(k+1))) * loc_azi[k]);
	      /* sum abcd factors */
	      tracer[ntracer].ax[j][k][0] = sin(loc_azi[k]) * (amp[j][k]*scale);
	      tracer[ntracer].ax[j][k][1] = cos(loc_azi[k]) * (amp[j][k]*scale);
	    }
	  }
	}else{
	  /* add to existing */
	  if((fabs(tracer[ntracer].lon-lon)> MAX(0.02,fabs(0.02*lon)))||
	     (fabs(tracer[ntracer].lat -lat) > MAX(0.02,fabs(0.02*lat)))){
	    fprintf(stderr,"%s: location error: %g %g %g %g tracer %i\n",
		    argv[0],tracer[ntracer].lon,lon,tracer[ntracer].lat,lat,i+1);
	    exit(-1);
	  }
	  for(j=0;j < nperiod;j++){
	    if(fabs(tracer[ntracer].T[j] - T[j]) > EPS_PREC){
	      fprintf(stderr,"%s: period error: j: %i T: %g %g tracer %i\n",
		      argv[0],j+1,tracer[ntracer].T[j],T[j],i+1);
	      exit(-1);
	    }
	    for(k=0;k<2;k++){
	      /* convert azimuth */
	      loc_azi[k] = DEG2RAD(((double)(2*(k+1))) * loc_azi[k]);
	      /* sum abcd factors */
	      tracer[ntracer].ax[j][k][0] += sin(loc_azi[k]) * (amp[j][k]*scale);
	      tracer[ntracer].ax[j][k][1] += cos(loc_azi[k]) * (amp[j][k]*scale);
	    }
	  }
	}
      }else{
	/* 
	   
	SAV MODE: READ IN TENSORS AND COMPUTE AVERAGE HERE
	
	*/
	/* 
	   read in sav tensor 
	*/
	if(!read_sym_6by6(sav,in)){
	  fprintf(stderr,"%s: sav read error\n",argv[0]);
	  exit(-1);
	}
	/* compute rphi terms but do not rotate the tensor as SAV
	   tensors are already rotated  */
	xp_from_lonlatz(lon,lat,depth,xp);
	for(j=0; j < nperiod ;j++){
	  /* 
	     get a factors and such. need to rotate from the regular
	     cartesian system of the sav tensors to the SW system
	  */
	  compute_phi_terms_from_Sav(sav, model->swsens[j].p,xp,loc_afactor,loc_bfactor,
				     &l_fac,&n_fac,
				     &a_fac,&c_fac,&f_fac,
				     (kernel),(kernel+1),(kernel+3),(kernel+4),(kernel+5),ba,hf,gl,ca,cn,
				     model,DREX_REG2SW_CONVENTION);
	  if(i==0){	
	    /*
	      first time around 
	    */
	    if(j == 0){
	      tracer = (struct rphi *)realloc(tracer,sizeof(struct rphi)*(ntracer+1));
	      if(!tracer)
		MEMERROR("");
	      tracer[ntracer].lon = lon;tracer[ntracer].lat=lat;
	    }
	    tracer[ntracer].T[j] =  model->swsens[j].p;
	    for(k=0;k<2;k++){
	      tracer[ntracer].ax[j][k][0] = loc_afactor[k*2+0] * scale;
	      tracer[ntracer].ax[j][k][1] = loc_afactor[k*2+1] * scale;
	    }
	  }else{
	    /* add to existing */
	    if((fabs(tracer[ntracer].lon-lon)> MAX(0.02,fabs(0.02*lon)))||
	       (fabs(tracer[ntracer].lat -lat) > MAX(0.02,fabs(0.02*lat)))){
	      fprintf(stderr,"%s: location error: %g %g %g %g tracer %i\n",
		      argv[0],tracer[ntracer].lon,lon,tracer[ntracer].lat,lat,i+1);
	      exit(-1);
	    }
	    if(fabs(tracer[ntracer].T[j] - model->swsens[j].p) > EPS_PREC){
	      fprintf(stderr,"%s: period error: j: %i T: %g %g tracer %i\n",
		      argv[0],j+1,tracer[ntracer].T[j],model->swsens[j].p,i+1);
	      exit(-1);
	    }
	    for(k=0;k < 2;k++){
	      tracer[ntracer].ax[j][k][0] += loc_afactor[k*2+0] * scale;
	      tracer[ntracer].ax[j][k][1] += loc_afactor[k*2+1] * scale;
	    }
	  }
	} /* end period loop */
	/* end sav mode branch */
      }
      ntracer++;
    } /* end tracer loop */
    if(i == 0){
      max_ntracer = ntracer;
    }else{
      if(max_ntracer != ntracer){
	fprintf(stderr,"%s: number of tracer mismatch (%i vs %i) at file %i, depth %g\n",
		argv[0],ntracer,max_ntracer,i+1,depth);
	exit(-1);
      }
    }
  } /* end depth loop */
  /* 

  average and output 

  */

  out = myopen("tracer.rphi.avg.dat","w","");
  for(i=0;i<ntracer;i++){
    fprintf(out,"%7.2f %7.2f ",
	    tracer[i].lon,tracer[i].lat);
    for(j=0;j < nperiod;j++){
      fprintf(out,"%5.1f ",tracer[i].T[j]);
      for(k=0;k < 2;k++){
	loc_azi[k] = RAD2DEG(atan2(tracer[i].ax[j][k][0],
				   tracer[i].ax[j][k][1])/((double)(2*(k+1))));
	while(loc_azi[k] >= 360.0)
	  loc_azi[k] -= 360.0;
	while(loc_azi[k] <= 0.0)
	  loc_azi[k] += 360.0;
	fprintf(out,"%6.1f %8.3e ",loc_azi[k],
		hypot(tracer[i].ax[j][k][0],tracer[i].ax[j][k][1]));
      }
    }
    fprintf(out,"\n");
  }
  fclose(out);
  fprintf(stderr,"%s: written to tracer.rphi.dat\n",argv[0]);
  return 0;
}
