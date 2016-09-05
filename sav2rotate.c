/* 

rotates tensors

$Id: sav2rotate.c,v 1.9 2011/04/12 06:19:09 becker Exp becker $


*/
#include "fstrack.h"



/* 

usage:  
 

cat sav.dat | sav2rotate alpha [beta, 0] [gamma, 0] [aniso_scale, 1]

format: 

lon lat depth upper_triangle_voigt

give angles in degree

 */
int main(int argc, char **argv)
{
  COMP_PRECISION sav[36],savr[36],alpha,beta,gamma,x[3],aniso_scale,cani[36],cmon[36],
    kmod,gmod,vel[9],symm_frac[6],tiaxis[6],ciso[36],chex[36],ctet[36],cort[36],ctri[36],
    sav_scc[36],scc_irot[9];
  my_boolean rescale = FALSE;
  int scca_old_mode = 1;
  FILE *in = stdin;
  int n,i,j,nr;
  /* 
     defaults
  */
  beta = gamma = 0.0;
  aniso_scale = 1.0;		/*  */
  if(argc < 2){
    fprintf(stderr,"%s alpha [beta, %g] [gamma, %g] [aniso_scale, %g](angles in degrees)\n",
	    argv[0],beta,gamma,aniso_scale);
    fprintf(stderr,"reads tensor in x y z upper_triangle_6_by_6 format, rotates, and prints new tensor to stdot\n");
    exit(-1);
  }
  sscanf(argv[1],FLT_FORMAT,&alpha);
  if(argc > 2)
    sscanf(argv[2],FLT_FORMAT,&beta);
  if(argc > 3){
    sscanf(argv[3],FLT_FORMAT,&gamma);
  }
  if(argc > 4)
    sscanf(argv[4],"%lf",&aniso_scale);
  if((aniso_scale < 0)||(aniso_scale>1)){
    fprintf(stderr,"%s: WARNING anisotropy scale should probably be between 0 and 1, but is %g\n",
	    argv[0],aniso_scale);
  }
  if(aniso_scale != 1.0){
    rescale = TRUE;
    fprintf(stderr,"%s: WARNING: only using %g%% of anisotropic component\n",
	    argv[0],aniso_scale * 100);
  }
  /* 
     parameters
  */
  fprintf(stderr,"%s: alpha: %g dip: %g gamma: %g (degrees)\n",
	  argv[0],alpha,beta,gamma);
  n=0;
  while(fscanf(in,THREE_FLT_FORMAT,x,(x+1),(x+2))==3){
    /* 
       
    read in 6,6 matrix
    
    */
    if(!read_sym_6by6(sav,in)){
      fprintf(stderr,"%s: error\n",argv[0]);
      exit(-1);
    }
    /* 
       rotate tensor out of horizontal 
    */
    drex_rotate_6x6_deg_ftrn(sav,savr,&alpha,&beta,&gamma);
    a_equals_b_vector(sav,savr,36);
    zero_small_entries(sav,36);
    if(rescale){
      /* 
	 decompose 
      */
      drex_decsym(sav,&kmod,&gmod,vel,symm_frac,tiaxis,ciso,chex,ctet,cort,cmon,ctri,
		  sav_scc,scc_irot,&scca_old_mode);
      a_equals_b_minus_c_vector(cani,sav,ciso,36);	/* anisotropic part */
      for(i=0;i<36;i++)
	sav[i] = ciso[i] + aniso_scale * cani[i];
    }
    fprintf(stdout,"%lg %lg %lg\t",x[0],x[1],x[2]);
    for(i=0;i<6;i++)
      for(j=i;j<6;j++){
	fprintf(stdout,"%lg ",sav[i*6+j]);
      }
    printf("\n");

    
    n++;
  }
}
  
