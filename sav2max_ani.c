/* 

   finds maximum anisotropy from angle scan, prints

   \Phi \xi \Xi
   

*/
#include "fstrack.h"



/* 



cat sav.dat | sav2max_ani

format: 

lon lat depth upper_triangle_voigt


 */
int main(int argc, char **argv)
{
  COMP_PRECISION sav[36],savr[36],alpha,beta,gamma;
  COMP_PRECISION afactor[4],bfactor[4],lmod,nmod,lon,lat,depth,xp[3],
    kernel[5],ba[2],hf[2],gl[2],ca[2],cn[2],amod,fmod,cmod;
  double xi,Xi,Phi,xi_max,Xi_max,Phi_max,da;
  double rho;
  struct prem_model *prem;
  struct mod *model;
  FILE *in = stdin;
  int n,coord_convention;
  int period;
  /* for sensitityv */
  model= (struct mod *)calloc(1,sizeof(struct mod));
  /* for density */
  prem= (struct prem_model *)calloc(1,sizeof(struct prem_model));
  /* 
     defaults
  */
  sprintf(model->sw_sens_file,"%s",SW_SENS_FILE);
  if(argc > 1){
    fprintf(stderr,"%s\n",argv[0]);
    fprintf(stderr,"reads tensor in x y z upper_triangle_6_by_6 format, and maximizes anisotropy\n");
    exit(-1);
  }
  period = 50;
  read_sw_sens(model);
  prem_read_model(PREM_MODEL_FILE,prem,FALSE);
  coord_convention = DREX_NO_ROTATION;

  
  n=0;
  while(fscanf(in,THREE_FLT_FORMAT,&lon,&lat,&depth)==3){
    prem_get_rho(&rho,1-depth/R_E,prem);
    xp_from_lonlatz(lon,lat,depth,xp);
    /* 
       
    read in 6,6 matrix
    
    */
    if(!read_sym_6by6(sav,in)){
      fprintf(stderr,"%s: error\n",argv[0]);
      exit(-1);
    }
    
    
    xi_max=0;
    Phi_max=0;
    Xi_max=0;
    da = 3;
    for(alpha=0;alpha<360;alpha+=da)
      for(beta=0;beta<90;beta+=da)
	for(gamma=0;gamma<90;gamma+=da){
    
	  drex_rotate_6x6_deg_ftrn(sav,savr,&alpha,&beta,&gamma);

	  compute_phi_terms_from_Sav(savr,(COMP_PRECISION)period,xp,
				     afactor,bfactor,&lmod,&nmod,&amod,&cmod,&fmod,
				     (kernel),(kernel+1),(kernel+3),(kernel+4),(kernel+5),
				     ba,hf,gl,ca,cn,model,
				     coord_convention);
	  xi = nmod/lmod;
	  Phi = sqrt(gl[0]*gl[0] + gl[1] * gl[1]);
	  Xi = sqrt(3.)*(sqrt(xi)-1)/sqrt(2+xi);
	  //fprintf(stderr,"%g %g %g %g %g %g\n",alpha,beta,gamma,Phi,xi,Xi);
	  if(xi>xi_max)xi_max=xi;
	  if(Xi>Xi_max)Xi_max=Xi;
	  if(Phi>Phi_max)Phi_max=Phi;
	}
    printf("%10.3f %10.3f %10.3f\n",Phi_max,xi_max,Xi_max);
    n++;
  }
  free(model);
  free(prem);
}
  
