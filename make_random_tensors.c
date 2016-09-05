#include "fstrack.h"

int main(int argc, char **argv)
{
  COMP_PRECISION depth,pressure,temperature,xp[3],lon,lat,alpha,gamma,beta,sig_dip,x1,x2,x3,x4;
  COMP_PRECISION sav[36],sav_o[36],sav_e[36],sav_tmp[36],sav_r[36],cn[2];
  COMP_PRECISION afactor[4],bfactor[4],ba[2],hf[2],gl[2],ca[2],amod,fmod,cmod,lmod,nmod,kernel[5],scale;
  struct mod *model;
  long seed;		/* for random number generator */
  FILE *in,*out;
  COMP_PRECISION 
    ol_frac = 0.7,		/* olivine fraction */
    vhrf = 0.0;			/* 0: voigt 0.5: VHR 1: Reuss  */
  /* 
     Euler angles in degrees for en rotation before 
     Voigt averaging for debugging modes  1 and 2
  */
  COMP_PRECISION euler[3] = ENSTATITE_EULER;
  int coord_convention = DREX_REG2SW_CONVENTION; /* simply rotate to SW convention */
  COMP_PRECISION period = 50;	/* not needed */

  if(argc < 3){
    fprintf(stderr,"%s rand4.dat aziout.dat \n",argv[0]);
    exit(-1);
  }

  seed = -1;			/*  */

  RAND_GEN(&seed);   		/* init */


  alpha = beta = gamma = 0.0;

  /* for sensitityv */
  model= (struct mod *)calloc(1,sizeof(struct mod));
  sprintf(model->sw_sens_file,"%s",SW_SENS_FILE);
  read_sw_sens(model);



  sig_dip = 10;			/* stanrard deviation for dip axes */

  depth = 50;			/* in km */
  
  pressure = pressure_model(depth);
  temperature = temperature_model(depth,1);
  fprintf(stderr,"depth: %11g T: %11g p: %11g\n",depth,temperature,pressure);

  
  /* get old tensor */
  drex_get_sav_constants(sav_o,DREX_CTP_OL_ESTEY,temperature,pressure); /* olivine */
  drex_get_sav_constants(sav_tmp,DREX_CTP_EN_ESTEY,temperature,pressure); /* enstatite */
  /* rotate enstatite */
  drex_rotate_6x6_deg_ftrn(sav_tmp,sav_e,(euler),(euler+1),(euler+2));zero_small_entries(sav_e,36);
  /* get VHR average */
  drex_st_vhr_avg_ftrn(sav_o,sav_e,&ol_frac,&vhrf,sav); 
  /* tensor is now in sav */

  in  = myopen(argv[1],"r",argv[0]);
  out = myopen(argv[2],"w",argv[0]);
  
  scale = 1.0;

  while(fscanf(in,"%lf %lf %lf %lf %lf %lf",
	       &lon,&lat,&x1,&x2,&x3,&x4) == 6){
    
    /* convert to spherical (not needed, really) */  
      xp_from_lonlatz(lon,lat,depth,xp);
      
      alpha = RAD2DEG(atan2(x1,x2));
      beta =  RAD2DEG(atan2(x3,hypot(x1,x2)));

      
      fprintf(stdout,"%11g %11g %11g\n",alpha,beta,scale);


      /* rotate tensor */
      drex_rotate_6x6_deg_ftrn(sav,sav_r,&alpha,&beta,&gamma);
      zero_small_entries(sav_r,36);
      
      
      
      /* get a factors and such */
      compute_phi_terms_from_Sav(sav_r,(COMP_PRECISION)period,xp,afactor,bfactor,&lmod,&nmod,&amod,&cmod,&fmod,
				 (kernel),(kernel+1),(kernel+3),(kernel+4),(kernel+5),ba,hf,gl,ca,cn,model,
				 coord_convention);
      /* xi = N/L and Gc/Gs anisotorpy */
      fprintf(out,"%11g %11g\t%11g %11g %11g\n",lon,lat,scale*nmod/lmod,scale*gl[0],scale*gl[1]);

  }


  fclose(out);
  fclose(in);
  return 0;
}
