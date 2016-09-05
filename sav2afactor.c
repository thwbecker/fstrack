
#include "fstrack.h"
/*

usage:

sav2afactor period [kernel_dir, ...] [ani_restrict, 0] [norot, 0]

reads in sav file and computes the 2phi and 4phi contribution for a
certain period surface wave using stored kernels

Love is not fully implemented yet 


*/

int main(int argc, char **argv)
{

  int n,i,period,coord_convention,
    norot = 0,			/* typically, rotate */
    restrict_ani = 0;		/* use all of the anisotropy */
  COMP_PRECISION afactor[4],bfactor[4],lon,lat,depth,lmod,nmod,xp[3],scale,sav[36],
    kernel[5],ba[2],hf[2],gl[2],ca[2],cn[2],amod,fmod,cmod;
  double rho;
  struct prem_model *prem;
  struct mod *model;
  my_boolean warned = FALSE;
  /* for sensitityv */
  model= (struct mod *)calloc(1,sizeof(struct mod));
  /* for density */
  prem= (struct prem_model *)calloc(1,sizeof(struct prem_model));

  if(argc < 2){
    fprintf(stderr,"%s period [kernel_dir, %s] [ani_restrict, %i] [no_rotation, %i]\n",
	    argv[0],SW_SENS_FILE,restrict_ani,norot);
    exit(-1);
  }
  // period to extract
  sscanf(argv[1],"%i",&period);
  /* directory for sensitivity kernels */
  if(argc > 2)
    sprintf(model->sw_sens_file,"%s",argv[2]);
  else
    sprintf(model->sw_sens_file,"%s",SW_SENS_FILE);
  if(argc > 3)
    sscanf(argv[3],"%i",&restrict_ani);
  if(argc > 4)
    sscanf(argv[4],"%i",&norot);
  /* 
     initialize the sensitivity kernels 
  */
  read_sw_sens(model);
  fprintf(stderr,"%s: computing 2phi/4phi a_i factors for period %i, kernels in %s, restrict ani: %i %s\n",
	  argv[0],period,model->sw_sens_file,restrict_ani,
	  (norot)?("not rotating"):("rotating into SW coordinate system"));
  /* init prem */
  prem_read_model(PREM_MODEL_FILE,prem,FALSE);

  if(norot == 1)
    coord_convention = DREX_NO_ROTATION;
  else if(norot == 2)
    coord_convention = DREX_SW_CART_CONVENTION;
  else{
    /* 
       tensor is in regular system, convert to surface waves , this is
       the regular operational mode
       
       in this case, the location of the input data does, of course,
       not matter

    */
    coord_convention = DREX_REG2SW_CONVENTION; 
  }

  n = 0;
  while(fscanf(stdin,THREE_FLT_FORMAT,&lon,&lat,&depth) == 3){
    if(depth < 0){
      fprintf(stderr,"%s: error: expected positive depths, read %g\n",
	      argv[0],depth);
      exit(-1);
    }
    /* 
       additional depth depending scaling (zero-ing out?)  
    */
    scale = restrict_ani_factor(restrict_ani,depth);
    /* 
       read in sav tensor
    */
    if(!read_sym_6by6(sav,stdin)){
      fprintf(stderr,"%s: sav read error, tracer %i\n",argv[0],n);
      exit(-1);
    }
    prem_get_rho(&rho,1-depth/R_E,prem);

    /* compute rphi terms and rotate the tensor from the regular
       system to the surface wave system  */
    xp_from_lonlatz(lon,lat,depth,xp);
    /* get a factors and such */
    compute_phi_terms_from_Sav(sav,(COMP_PRECISION)period,xp,
			       afactor,bfactor,&lmod,&nmod,&amod,&cmod,&fmod,
			       (kernel),(kernel+1),(kernel+3),(kernel+4),(kernel+5),
			       ba,hf,gl,ca,cn,model,
			       coord_convention);
    /* rescale */
    if(scale != 1.0){
      if(!warned){
	fprintf(stderr,"%s: WARNING: rescaling depth %g with %g\n",argv[0],depth,scale);
	warned = TRUE;
      }
      for(i=0;i<4;i++){
	afactor[i] *= scale;
	bfactor[i] *= scale;
      }
    }
    /* print the BA HF and GL factors */
    //fprintf(stderr,"bca=%8.5f; hcf=%8.5f; gcl=%8.5f;\n",
    //	    ba[0],hf[0],gl[0]);
    //fprintf(stderr,"bsa=%8.5f; hsf=%8.5f; gsl=%8.5f;\n",
    //	    ba[1],hf[1],gl[1]);
    //
    // to compare with montagne & Nataf 1986
    // cat tensors/peselnick_2d.sav | sav2afactor 50 progs/src/fstrack/sw_sens/PREMb/fsens 0 1 | gawk '{printf("%.3f %.3f\n",$12,$13)}'
    // could reproduce table 2
    //
    fprintf(stdout,"%12g %12g %12g\t\t%12.5e %12.5e %12.5e %12.5e\t%8.5f %8.5f %8.5f %8.5f %8.5f %8.5f\t%.7e %.7e %.7e\t%11g %.7e %.7e %.7e %.7e\t%.7e\t%12.5e %12.5e %12.5e %12.5e\t%.5e %.5e %.5e %.5e\n",
	    lon,lat,depth,	/* 1..3 */
	    afactor[0],afactor[1], /* 4..7    Rayleigh wave afactors */
	    afactor[2],afactor[3],
	    /* 8: B_c/A B_s/A 10: H_c/F H_s/F 12: G_c/L G_s/L  */
	    ba[0],ba[1],hf[0],hf[1],gl[0],gl[1], 
	    /* 14 = N/L = xi  15 = C/A =phi   16 = eta */
	    nmod/lmod,cmod/amod,fmod/(amod-2.0*lmod), /* 14...16 */
	    /* 17 = rho [kg/m^3]  18 = A  19 = C   20 = N        21 = L [Pa] 22 T */
	    rho,amod*1e9,cmod*1e9,nmod*1e9,lmod*1e9,temperature_model(depth,1), /* 17 .. 22 */
	    bfactor[0],bfactor[1], /* 23...26    Love wave bfactors */
	    bfactor[2],bfactor[3],
	    /* 27: C_c/N C_s/N 29: C_c/A C_s/A */
	    cn[0],cn[1],ca[0],ca[1]);	/* 27 - 30 */
    n++;
  }
  fprintf(stderr,"%s: read %i tracer savs \n",argv[0],n);
  return 0;
}
