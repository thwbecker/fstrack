#include "fstrack.h"
/* 


read in depth(for pressure) temperature/mean_tempare and angle variations and create a
tensor from those


*/
int main(int argc, char **argv)
{
  COMP_PRECISION depth,pressure,temperature,xp[3],lon,lat,alpha,gamma,beta;
  COMP_PRECISION sav[36],sav_o[36],sav_e[36],sav_tmp[36],sav_r[36];
  COMP_PRECISION T;
  struct mod *model;
  COMP_PRECISION 
    ol_frac = 0.7,		/* olivine fraction */
    vhrf = 0.0;			/* 0: voigt 0.5: VHR 1: Reuss  */
  /* 
     Euler angles in degrees for en rotation before 
     Voigt averaging for debugging modes  1 and 2
  */
  COMP_PRECISION euler[3] = ENSTATITE_EULER;

  if(argc > 1){
    fprintf(stderr,"%s \nread lon lat depth T/Tmean dip from stdin, write tensor to stdout",argv[0]);
    exit(-1);
  }


  alpha = beta = gamma = 0.0;

  /* for sensitityv */
  model= (struct mod *)calloc(1,sizeof(struct mod));
  sprintf(model->sw_sens_file,"%s",SW_SENS_FILE);
  read_sw_sens(model);

  /* read lon lat depth T/Tmean dip[degrees] */
  while(fscanf(stdin,"%lf %lf %lf %lf %lf",&lon,&lat,&depth,&T,&beta)==5){
  
    pressure = pressure_model(depth);
    
    temperature = temperature_model(depth,1) * T;
  
    /* get old tensor */
    drex_get_sav_constants(sav_o,DREX_CTP_OL_ESTEY,temperature,pressure); /* olivine */
    drex_get_sav_constants(sav_tmp,DREX_CTP_EN_ESTEY,temperature,pressure); /* enstatite */
    /* rotate enstatite */
    drex_rotate_6x6_deg_ftrn(sav_tmp,sav_e,(euler),(euler+1),(euler+2));zero_small_entries(sav_e,36);
    /* get VHR average */
    drex_st_vhr_avg_ftrn(sav_o,sav_e,&ol_frac,&vhrf,sav); 
    /* tensor is now in sav */
    
    xp_from_lonlatz(lon,lat,depth,xp);

    /* rotate tensor */
    drex_rotate_6x6_deg_ftrn(sav,sav_r,&alpha,&beta,&gamma);
    zero_small_entries(sav_r,36);

    fprintf(stdout,"%g %g %g\t",lon,lat,depth);
    print_sym_6by6(sav_r,stdout);
    
  }

  return 0;
}
