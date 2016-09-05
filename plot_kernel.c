
#include "fstrack.h"
/*

usage:

plot_kernel period 


*/

int main(int argc, char **argv)
{

  int i,j,period,restrict_ani = 0,n;
  COMP_PRECISION afactor[4],bfactor[4],
    lmod,nmod,amod,fmod,cmod,xp[3],scale,sav[36],sav_ens[36],sav_oli[36],depth,savr[36],
    ol_frac=0.7,vhrf = 0.0;	/* olivine fraction, averaging */
  COMP_PRECISION euler[3] = ENSTATITE_EULER; /* for en rotation */
  COMP_PRECISION RAK, RFK, RLK,LLK,LNK, ba[2],hf[2],gl[2],ca[2],cn[2];
  my_boolean warned = FALSE;
  struct mod *model;
  int tmode = 2;		/* 0: ol/en at 200km conditions 
				   1: Ji et al tensor
				   2: peselnick tensor

				*/

  /* for sensitityv */
  model= (struct mod *)calloc(1,sizeof(struct mod));
  if(argc < 2){
    fprintf(stderr,"%s period [kernel_dir, %s] [ani_restrict, %i]\n",
	    argv[0],SW_SENS_FILE,restrict_ani);
    exit(-1);
  }
  // period to extract
  sscanf(argv[1],"%i",&period);
  fprintf(stderr,"%s: computing kernels for period %is\n",
	  argv[0],period);
  /* directory for sensitivity kernels */
  if(argc > 3)
    sprintf(model->sw_sens_file,"%s",argv[2]);
  else
    sprintf(model->sw_sens_file,"%s",SW_SENS_FILE);
  if(argc > 4)
    sscanf(argv[3],"%i",&restrict_ani);
  /* initialize the sensitivity kernels */
  read_sw_sens(model);
  /*  */
  fprintf(stderr,"%s: computing 2phi/4phi a_i factors for period %i, kernels in %s\n",
	  argv[0],period,model->sw_sens_file);
  /* tensors */
  switch(tmode){
  case 0:			/* ol/en mix */
    drex_get_sav_constants(sav_oli,DREX_CTP_OL_ESTEY,1573.0,5.0); /* T, p (200km) */
    drex_get_sav_constants(sav_ens,DREX_CTP_EN_ESTEY,1573.0,5.0); /* T, p @ room */
    /* combine */
    drex_rotate_6x6_deg_ftrn(sav_ens,savr,(euler),(euler+1),(euler+2)); /* rotate ens like below */
    zero_small_entries(savr,36);
    drex_st_vhr_avg_ftrn(sav_oli,savr,&ol_frac,&vhrf,sav); /* get VHR average */
    break;
  case 1:
  case 2:
    if(tmode == 1){		/* Ji et al tensor */
      sav[0*6+0] = 2.2444; sav[0*6+1] = 0.6727; sav[0*6+2] = 0.7071; sav[0*6+3] = -0.0028; 
      sav[0*6+4] = -0.0033; sav[0*6+5] = -0.0094; sav[1*6+1] = 1.9313; sav[1*6+2] = 0.6869; 
      sav[1*6+3] = 0.0029; sav[1*6+4] = 0.0019; sav[1*6+5] = -0.0088; sav[2*6+2] = 2.1052; 
      sav[2*6+3] = 0.0057; sav[2*6+4] = -0.0025; sav[2*6+5] = -0.0025; sav[3*6+3] = 0.6552; 
      sav[3*6+4] = -0.0065; sav[3*6+5] = 0; sav[4*6+4] = 0.7303; sav[4*6+5] = -0.0005; 
      sav[5*6+5] = 0.6865;
    }else if(tmode == 2){	/* peselnick tensor */
      sav[0*6+0] = 2.3648; sav[0*6+1] = 0.7253; sav[0*6+2] = 0.7228; sav[0*6+3] = -0.0008; 
      sav[0*6+4] = -0.0196; sav[0*6+5] = 0; sav[1*6+1] = 2.2081; sav[1*6+2] = 0.7187; sav[1*6+3] = 0.0169; 
      sav[1*6+4] = 0.0164; sav[1*6+5] = -0.004; sav[2*6+2] = 2.2016; sav[2*6+3] = 0.0182; 
      sav[2*6+4] = -0.0024; sav[2*6+5] = 0.0041; sav[3*6+3] = 0.7494; sav[3*6+4] = -0.0058; 
      sav[3*6+5] = -0.0103; sav[4*6+4] = 0.7921; sav[4*6+5] = 0.0128; sav[5*6+5] = 0.7877; 
    }
    for(i=0;i<6;i++)		/* fill in lower triangle part */
      for(j=0;j<i;j++)
	sav[i*6+j] = sav[j*6+i];
    break;
  default:
    fprintf(stderr,"%s: tensor mode %i undefined\n",
	    argv[0],tmode);
    exit(-1);
    break;
  }
  print_6by6_nice(sav,stderr);
  /* 


  THIS EXPECTS THAT THE SAV TENSORS ARE IN THE REGULAR CARTESIAN SYSTEM X = S, Y = E, Z = U
  AND WILL ROTATE THEM INTO THE SURFACE WAVE SYSTEM


  */



  n = 0;
  for(depth=0.0;depth <= 410;depth+=1){
    /* additional depth depending scaling (zero-ing out?)  */
    scale = restrict_ani_factor(restrict_ani,depth);
    /* 
       need radius  
    */
    xp_from_lonlatz(0.0,0.0,depth,xp);
    compute_phi_terms_from_Sav(sav,(COMP_PRECISION)period,xp,afactor,bfactor,&lmod,&nmod,
			       &amod,&cmod,&fmod,
			       &RAK,&RFK,&RLK,&LLK,&LNK,ba,hf,gl,ca,cn,
			       model,DREX_REG2SW_CONVENTION);
    if(!n){
      fprintf(stderr,"%s:\tG/L\t\tC/A\t\tB/A\t\tH/F\n",argv[0]);
      fprintf(stderr,"%s:C\t%11g\t%11g\t%11g\t%11g\n",argv[0],gl[0],ca[0],ba[0],hf[0]);
      fprintf(stderr,"%s:S\t%11g\t%11g\t%11g\t%11g\n",argv[0],gl[1],ca[1],ba[1],hf[1]);
      fprintf(stderr,"%s: \t%11g\t%11g\t%11g\t%11g\n",argv[0],
	      hypot(gl[0],gl[1]),hypot(ca[0],ca[1]),
	      hypot(ba[0],ba[1]),hypot(hf[0],hf[1]));
    }
    if(scale != 1.0){
      if(!warned){
	fprintf(stderr,"%s: WARNING: rescaling depth %g with %g\n",argv[0],depth,scale);
	warned = TRUE;
      }
      for(i=0;i<4;i++)
	afactor[i] *= scale;
    }
    /*   1          2   3   4                         5   6   7   8                9     10                */
    /* depth       k_A k_F k_L                       A_1 A_2 A_3 A_4              A2p   a4p*/
    fprintf(stdout,"%11g\t%12.5e %12.5e %12.5e\t\t%12.5e %12.5e %12.5e %12.5e\t%12.5e %12.5e\n",
	    depth,
	    RAK,RFK,RLK,
	    afactor[0],afactor[1],afactor[2],afactor[3],
	    hypot(afactor[0],afactor[1]),
	    hypot(afactor[2],afactor[3]));

    n++;
  }
  return 0;
}
