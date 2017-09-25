/* 

usage:

sav2cijkl [sav.dat, stdin] [density in g/cm^3, 3.35] [beta, 0] [aniso_scale, 1]

convert Sav stiffness tensor format to normalized Cijkl format

input format:

lon lat depth upper_triangl_C_ij (j fast), i.e.
			  
SAV_11 SAV_12 SAV_13 ... SAV_16 \
       SAV_22 SAV_23 ... SAV_26 \
       ... 
                         SAV_66

output: 

Cijkl normalized by density in 1..3 format

options:

density [g/cm^3]

densities: 3.35  for olivine, approx
           3.353 for olivine/enstatite(30%)


beta: if != 0: rotate tensor with (0,beta,0) angles (degree) before conversion

aniso_scale: output tensor will be iso_comp + anis_scale * aniso_comp, i.e. full 
tensor for aniso_scale = 1, isotropic tensor for aniso_scale = 0


$Id: sav2cijkl.c,v 1.10 2016/09/05 04:44:47 becker Exp $


*/
#include "fstrack.h"

int main(int argc, char **argv)
{
  VERA_PREC sav[36],cij[81],rho,x[3],alpha,beta,gamma,savr[36],aniso_scale,cani[36],cmon[36],
    kmod,gmod,vel[9],symm_frac[6],tiaxis[6],ciso[36],chex[36],ctet[36],cort[36],ctri[36],
    sav_scc[36],scc_irot[9];
  my_boolean rotate = FALSE, rescale = FALSE;
  int i,nr,n,j,iop=FTRN_STDOUT;
  FILE *in;
  int scca_old_mode = 1;
  /* 
     defaults

  */
  /* 
     densities

     3.35  for olivine
     3.353 for olivine/enstatite(30%)

  */
  rho = 3.353;			/* density in g/cm^3 */
  /* rotation angles */
  alpha = beta = gamma = 0.0;
  /* anisotropy factor */
  aniso_scale = 1.0;
  if(argc > 1){
    /* first argument will be filename */
    in = myopen(argv[1],"r",argv[0]);
  }else{
    in = stdin;
  }
  if(argc > 2)
    sscanf(argv[2],"%lf",&rho);
  if(argc > 3)
    sscanf(argv[3],"%lf",&beta);
  if(argc > 4)
    sscanf(argv[4],"%lf",&aniso_scale);
  if(fabs(beta) > EPS_COMP_PREC){
    fprintf(stderr,"%s: WARNING: rotating tensor by %g degrees down from horizontal (beta angle)\n",
	    argv[0],beta);
    rotate = TRUE;
  }
  if((aniso_scale < 0)||(aniso_scale>1)){
    fprintf(stderr,"%s: WARNING: anisotropy scale should typically be between 0 and 1, but is %g\n",
	    argv[0],aniso_scale);
  }
  if(aniso_scale != 1.0){
    rescale = TRUE;
    fprintf(stderr,"%s: WARNING: using %g%% of anisotropic component\n",
	    argv[0],aniso_scale * 100);
  }
  /* 
     input 
  */
  n=0;
  while(fscanf(in,THREE_FLT_FORMAT,x,(x+1),(x+2))==3){
    if(!read_sym_6by6(sav,in)){
      fprintf(stderr,"%s: error\n",argv[0]);
      exit(-1);
    }
    if(rotate){		/* rotate tensor out of horizontal */
      fprintf(stderr,"%s: WARNING: rotating by %g, %g, %g\n",
	      argv[0],alpha,beta,gamma);
      /* rotate, angles are in degree */
      drex_rotate_6x6_deg_ftrn(sav,savr,&alpha,&beta,&gamma);
      a_equals_b_vector(sav,savr,36);
      zero_small_entries(sav,36);
    }
    if(rescale){
      /* decompose */
      drex_decsym(sav,&kmod,&gmod,vel,symm_frac,tiaxis,ciso,chex,ctet,cort,cmon,ctri,
		  sav_scc,scc_irot,&scca_old_mode);
      a_equals_b_minus_c_vector(cani,sav,ciso,36);	/* anisotropic part */
      for(i=0;i<36;i++)
	sav[i] = ciso[i] + aniso_scale * cani[i];
    }
    //print_6by6_nice(sav,stdout);fprintf(stderr,"\n");
    /* print Cij to stderr */
    //print_cij_vera_sorted(sav, stderr);
    /* 
       convert to normalized cijkl format 
    */
    vera_sav_to_cijkl_ftrn(sav,&rho,cij);
    /* print to Cijkl stdout */
    vera_print_cijkl_ftrn(cij,&iop);
    n++;
  }
  fprintf(stderr,"%s: converted %i tensors using rho: %g\n",
	  argv[0],n,rho);
}
  
