/*

  either print table of PREM values (if PRINT_TABLE is defined) 
  or give them at a certain radius (else)

  if USE_DEPTH is defined for not defined PRINT_TABLE, will use
  depth instead of radius


  reads polynomial coefficients of anisotropic PREM from a model
  file specified with -DPREM_MODEL_FILE
  
  
  $Id: readprem.c,v 1.9 2004/03/02 19:56:22 becker Exp becker $

*/

#include "prem.h"

int main(int argc, char **argv)
{

  hc_boolean init = FALSE;
  struct prem_model prem;
  double r,vp,vs,vsv,vsh,vpv,vph,rho,eta,qmu,qkappa;
  char filename[301];
#ifndef PRINT_TABLE
  int header_printed = 0;
  strncpy(filename,PREM_MODEL_FILE,300);
  /* 
     
  this part of the code reads in depth/radius  levels

  */
  if(argc != 1){
#ifdef USE_DEPTH		/* depth mode */
    fprintf(stderr,"%s reads z [km] depths (0=surface...%g) from stdin and \nreturns\n\nv_p[m/s] v_s[m/s] rho[kg/m^3] v_pv[m/s] v_ph[m/s] v_sv[m/s] v_sh[m/s] eta vsVgt\n\nfrom PREM\n",
	    argv[0],PREM_RE_KM);
#else  /* radius mode */
    fprintf(stderr,"%s reads r [km] radii (0...%g=surface) from stdin and \nreturns\n\nv_p[m/s] v_s[m/s] rho[kg/m^3] v_pv[m/s] v_ph[m/s] v_sv[m/s] v_sh[m/s] eta vsVgt\n\nfrom PREM\n",
	    argv[0],PREM_RE_KM);
#endif
    exit(-1);
  }
  while(fscanf(stdin,"%lf",&r)==1){
    if(r < 0 || r > PREM_RE_KM){
#ifdef USE_DEPTH
      fprintf(stderr,"range error: z: %g (gotta be between 0 (surface) and %g)\n",
	      r,PREM_RE_KM);
#else
      fprintf(stderr,"range error: r: %g (gotta be between %g (surface) and 0)\n",
	      r,PREM_RE_KM);
#endif
      exit(-1);
    }
#ifdef USE_DEPTH
    r = PREM_RE_KM - r;
#endif
    r *= 1e3;
#else
    /* 
       this part for the table mode 
    */
  double rmin=0.0,rmax=PREM_RE_KM,dr=50.,rhoint,vint,oldr,oldrho,dv;
  int i;
  oldr = oldrho = 0.0;
  if((argc>5)||((argc>1)&&(strcmp(argv[1],"-h")==0))){
#ifndef USE_DEPTH
    fprintf(stderr,"%s [rmin, %g] [rmax, %g] [dr, %g] [prem_model, %s]\nprints table of PREM values\n",
	    argv[0],rmin,rmax,dr,PREM_MODEL_FILE);
#else
    fprintf(stderr,"%s [zmin, %g] [zmax, %g] [dz, %g] [prem_model, %s]\nprints table of PREM values\n",
	    argv[0],rmin,rmax,dr,PREM_MODEL_FILE);
#endif
    exit(-1);
  }
  if(argc>1)
    sscanf(argv[1],"%lf",&rmin);
  if(argc>2)
    sscanf(argv[2],"%lf",&rmax);
  if(argc>3)
    sscanf(argv[3],"%lf",&dr);
  if(argc>4){
    strncpy(filename,argv[4],300);
  }else
    strncpy(filename,PREM_MODEL_FILE,300);

#ifdef USE_DEPTH
  rmin = PREM_RE_KM - rmin;
  rmax = PREM_RE_KM - rmax;
  oldr = rmin ;rmin = rmax;rmax=oldr; /*  flip */
#endif  
  /* scale to meters */
  rmin *= 1e3;
  rmax *= 1e3;
  if(dr == 0)dr = 1;
  dr   *= 1e3;
#endif


  if(!init){
    // read in PREM model
    fprintf(stderr,"%s: using PREM model in %s\n",argv[0],filename);
    if(prem_read_model(filename,&prem,TRUE)){
      fprintf(stderr,"%s: couldn't init model\n",argv[0]);
      exit(-1);
    }
    init = TRUE;
  }

#ifndef PRINT_TABLE
  /* 
     single depth/radius mode
  */
  if(!header_printed){
    fprintf(stderr,"# v_p [m/s]     v_s [m/s]    rho[kg/m3]  v_pv [m/2]   v_ph [m/s]    v_sv [m/s]   v_sh [m/s]   eta        q_mu         q_kappa        v_s^Vgt [m/s]\n");
    header_printed = 1;
  }
  prem_get_values(&rho,&vs,&vsv,&vsh,&vp,&vpv,&vph,&eta,
		  &qmu,&qkappa,r,&prem);
  fprintf(stdout,"%12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.7f %12.5e %12.5e %12.4f\n",
	  vp,vs,rho,vpv,vph,vsv,vsh,eta,qmu,qkappa,
	  prem_vs_voigt(vsh, vsv, vph, vpv, eta));
  }
#else
  /* 
     
  table mode 

  */
  fprintf(stderr,"#  r [km]    z [km]         v_p [m/s]     v_s [m/s]    rho[kg/m3]  v_pv [m/2]   v_ph [m/s]    v_sv [m/s]   v_sh [m/s]   eta        q_mu         q_kappa        v_s^Vgt [m/s]\n");
  for(vint=rhoint=0.0,i=0,r=rmin;r<=rmax+1.0e-4;r+=dr,i++){
    prem_get_values(&rho,&vs,&vsv,&vsh,&vp,&vpv,&vph,&eta,
		    &qmu,&qkappa,r,&prem);
    fprintf(stdout,"%12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.7f %12.5e %12.5e %12.4f\n",
	    r/1e3,PREM_RE_KM-r/1e3,
	    vp,vs,rho,vpv,vph,vsv,vsh,eta,qmu,qkappa,
	    prem_vs_voigt(vsh, vsv, vph, vpv, eta));
    // integrate density 
    if(i > 0){
      dv = (r-oldr) * pow((r+oldr)/2.0,2);
      rhoint += (rho+oldrho)/2.0 * dv;
      vint += dv;
    }
    oldr=r;
    oldrho=rho;
  }
  //fprintf(stdout,"# M=(\\int \\rho r^2 dr) rm=(\\int \\rho r^2 dr)/(\\int r^2 dr) M/3.19632e+23(M_m)\n# %g %g %g\n",
  //rhoint,rhoint/vint,rhoint/3.19632e23);
#endif

  return 0;
}
