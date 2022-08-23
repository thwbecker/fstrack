#include "fstrack.h"
#include <math.h>
/* 

routines dealing with surface wave sensitivities


$Id: sens_handling.c,v 1.11 2010/12/31 18:59:17 becker Exp $


*/

/* 

compute the 2phi and 4phi sensitivities given a surfave wave of period
p at location xp[r,theta,phi] and an elastic tensor Sav in the global
Cartesian system.

if rotate_tensor is true, the tensor will be rotated into the local
system where z is up, x is south and y is East. phi is the azimuth
clockwise from North

For Rayleigh waves

a[0] = 2 phi cos
a[1] = 2 phi sin
a[2] = 4 phi cos
a[3] = 4 phi sin

For Love waves

b[0] = 2 phi cos
b[1] = 2 phi sin
b[2] = 4 phi cos
b[3] = 4 phi sin


WARNING: Love waves not implemented yet !!!


this also computes the Vsv and Vsh terms l and n

Sav gets passed C-style 

*/
void compute_phi_terms_from_Sav(COMP_PRECISION *Sav, /* C_ij tensor  */
				COMP_PRECISION p, /* period, in seconds */
				COMP_PRECISION *xp, /* location (needed for rotation) */
				COMP_PRECISION *a, /* a[4] factors for 2phi/4phi for Rayleigh */
				COMP_PRECISION *b, /* a[4] factors for 2phi/4phi for Love  */
				COMP_PRECISION *lloc, /* effective L,N,A,C,F moduli */
				COMP_PRECISION *nloc,
				COMP_PRECISION *aloc,
				COMP_PRECISION *cloc,
				COMP_PRECISION *floc,
				COMP_PRECISION *RAK, /* Rayleigh kernels:  A/T dT/dA*/
				COMP_PRECISION *RFK, /* F/T dT/dF  */
				COMP_PRECISION *RLK, /* L/T (dT/dL) */
				COMP_PRECISION *LLK, /* Love L/T (dT/dL) */
				COMP_PRECISION *LNK, /* Love N/T (dT/dN) */
				COMP_PRECISION *ba, /* B_{c,s} / A*/
				COMP_PRECISION *hf, /* H_{c,s} / F */
				COMP_PRECISION *gl, /* G_{c,s} / L */
				COMP_PRECISION *ca, /* C_{c,s} / A */
				COMP_PRECISION *cn, /* C_{c,s} / N */
				struct mod *model,
				int coord_convention) /* should we
							 rotate this
							 tensor into
							 a surface
							 wave system?
							 
							 call from within code with 
							 
							 DREX_REG2SW_CONVENTION

							 if the sav
							 tensor has
							 already been
							 rotated to
							 the
							 DREX_REG_CART_CONVENTION, then use 


							 DREX_REG2SW_CONVENTION

						      */
{
  COMP_PRECISION Savr[36],alpha,beta,gamma;
  if(!model->sw_sens_init){
    fprintf(stderr,"compute_phi_terms_from_Sav: error: sensitivity kernels not initialized\n");
    exit(-1);
  }
  /* 
     
     for Rayleigh 2phi, we need A, F, and L
     for Rayeligh 4phi, we need A 
     obtain the interpolated A, F, and L kernel values at this depth

     for Love wave 2phi, we need L 
     for Love wave 4phi, we need N 


  */
  interpolate_sens(p,xp[FSTRACK_R],model->swsens,RAK,RFK,RLK,LLK,LNK);

  if(coord_convention != DREX_NO_ROTATION){
    /* 
       
    rotate the stiffness tensor into the local cartesian system
    for surface waves, the convention is that x is north, y is 
    east, and z is down
    
    */
    drex_calc_euler_rad_from_xp(xp,&alpha,&beta,&gamma,
				&coord_convention);
    /* 
       
    this assumes that Sav is stored FORTRAN style
    
    */
    drex_rotate_6x6_rad_ftrn(Sav,Savr,&alpha,&beta,&gamma);
    /* 
       obtain the elasticity ratios in the rotated reference frame, Sav
       is FORTRAN style
    */
    drex_compute_swpar_ftrn(Savr,gl,cn,ca,ba,hf,lloc,nloc,
			    aloc,cloc,floc);
  }else{
    /* no rotation */
    drex_compute_swpar_ftrn(Sav,gl,cn,ca,ba,hf,lloc,nloc,
			    aloc,cloc,floc);
  }
  /* 
     get the 2,4 phi parameters for Rayleigh waves
  */
  a[0] =  ba[0] * (*RAK) + hf[0] * (*RFK) + gl[0] * (*RLK); /* cos(2phi) */
  a[1] =  ba[1] * (*RAK) + hf[1] * (*RFK) + gl[1] * (*RLK); /* sin(2phi) */
  a[2] =  ca[0] * (*RAK);		/* cos(4phi) */
  a[3] =  ca[1] * (*RAK);		/* sin(4phi) */
  /* 2phi and 4phi for love waves */
  //fprintf(stderr,"%g %g %g %g %g %g\n",gl[0],gl[1],cn[0],cn[1],*LLK,*LNK);
  b[0] = -gl[0] * (*LLK);	/* cos(2phi) */
  b[1] = -gl[1] * (*LLK);	/* sin(2phi) */
  b[2] = -cn[0] * (*LNK);	/* cos(4phi) */
  b[3] = -cn[1] * (*LNK);	/* sin(4phi) */
  
}
/* 
   given sensitivity structures, find the Rayleigh A, F, L, and Love L and N parameters for period p 
   at radius r (0...1)
*/
void interpolate_sens(COMP_PRECISION p, COMP_PRECISION r,
		      struct swss *sens, COMP_PRECISION *RA,
		      COMP_PRECISION *RF,COMP_PRECISION *RL,
		      COMP_PRECISION *LL, COMP_PRECISION *LN)
{
  int i,j,il;
  COMP_PRECISION fac1,fac2;
  
  for(il=0;il < N_SW_SENS;il++)
    if(fabs(p - sens[il].p)<EPS_PREC)
      break;
  if(il == N_SW_SENS){
    fprintf(stderr,"interpolate_sens: error: period %g not in tables\n",p);
    exit(-1);
  }
  if(r < sens[il].l[0].r){	/* don't interpolate to deeper depths */
    fprintf(stderr,"interpolate_sens: error: radius r of %g is to deep compared to %g for period %g\n",
	    r, sens[il].l[0].r,sens[il].p);
    exit(-1);
  }
  /* linearly interpolate */
  j=0;
  while( (j < sens[il].n) && (sens[il].l[j].r < r))
    j++;
  if(j == 0)
    j = 1;
  i = j - 1;
  fac1 =(r - sens[il].l[i].r)/(sens[il].l[j].r - sens[il].l[i].r);
  fac2 = 1.0 - fac1;
  /* interpolated values */
  /* rayleigh */
  *RA =  fac1  * sens[il].l[j].RA + fac2 * sens[il].l[i].RA;
  *RF =  fac1  * sens[il].l[j].RF + fac2 * sens[il].l[i].RF;
  *RL =  fac1  * sens[il].l[j].RL + fac2 * sens[il].l[i].RL;
  /* love */
  *LL =  fac1  * sens[il].l[j].LL + fac2 * sens[il].l[i].LL;
  *LN =  fac1  * sens[il].l[j].LN + fac2 * sens[il].l[i].LN;

}



/* 


read in the sensitivity kernels in z[km] A F L format for all periods


*/
void read_sw_sens(struct mod *model)
{
  FILE *in;
  int i;
  char fname[STRLEN];
  COMP_PRECISION p[N_SW_SENS] = SW_SENS_PERIODS;
  if(model->sw_sens_init){
    fprintf(stderr,"read_sw_sens: error: call only once\n");
    exit(-1);
  }
  model->swsens=(struct swss *)malloc(N_SW_SENS*sizeof(struct swss));
  if(!model->swsens)
    MEMERROR("");
  if(model->verbose)
    fprintf(stderr,"\nWARNING\nLove kernels not implemented yet\n\n");

  for(i=0;i < N_SW_SENS;i++){
    model->swsens[i].p = p[i];	/* period */
    sprintf(fname,"%s/%s.%g.dat",getenv("HOME"),model->sw_sens_file,model->swsens[i].p);
    in = myopen(fname,"r","read_sw_sens");
    model->swsens[i].n = 0;
    model->swsens[i].l = (struct swss_layer *)malloc(sizeof(struct swss_layer));
    if(!model->swsens[i].l)
      MEMERROR("");
    /* Rayleigh */
    /* expects depth A F L */
    while(fscanf(in,FOUR_FLT_FORMAT,&model->swsens[i].l[model->swsens[i].n].r,
		 &model->swsens[i].l[model->swsens[i].n].RA,
		 &model->swsens[i].l[model->swsens[i].n].RF,
		 &model->swsens[i].l[model->swsens[i].n].RL)==4){
      /* Love */
      model->swsens[i].l[model->swsens[i].n].LL = model->swsens[i].l[model->swsens[i].n].LN = 1.0;
      /* 
	 convert depth in [km] to radius 
      */
      model->swsens[i].l[model->swsens[i].n].r = 
	ND_RADIUS( model->swsens[i].l[model->swsens[i].n].r);
      model->swsens[i].n += 1;
      model->swsens[i].l = (struct swss_layer *)
	realloc(model->swsens[i].l,(model->swsens[i].n+1)*sizeof(struct swss_layer));
      if(!model->swsens[i].l)MEMERROR("");
      if(model->swsens[i].n > 1){
	if(model->swsens[i].l[model->swsens[i].n-1].r < 
	   model->swsens[i].l[model->swsens[i].n-2].r){
	  fprintf(stderr,"read_sw_sens: error: file %s, depths should be sorted ascendingly, in km, 0 = surface\n",
		  fname);
	  exit(-1);
	}
      }
    }
    fclose(in);
  }
  if(model->verbose){
    fprintf(stderr,"read_sw_sens: initialized %i sensitivity kernels with ",N_SW_SENS);
    for(i=0;i < N_SW_SENS;i++)
      fprintf(stderr,"%gs (%i) ",model->swsens[i].p,model->swsens[i].n);
    fprintf(stderr," periods\n");
  }
  
  model->sw_sens_init = TRUE;
}

/* 
   compute the azimuths of the maximum 2 phi and 4 phi (azi[2])
   orientation anisotropy given four afactors
   (degrees)
*/
void compute_azideg_amp_from_afactors(COMP_PRECISION *afactor,
				      COMP_PRECISION *azi,
				      COMP_PRECISION *amp)
{
  int i;
  COMP_PRECISION two_phi, four_phi;
  /* 2phi */

  two_phi = atan2(afactor[1],afactor[0]);
  azi[0] = two_phi / 2.0;
  //amp[0] = afactor[0] * cos(two_phi) + afactor[1] * sin(two_phi);
  /* that should be faster */
  amp[0] = hypot( afactor[0] , afactor[1] );

  /* 4phi */
  four_phi = atan2(afactor[3],afactor[2]);
  azi[1] = four_phi / 4.0;
  //amp[1] = afactor[2] * cos(four_phi) + afactor[3] * sin(four_phi);
  amp[1] = hypot( afactor[2] , afactor[3] );

  for(i=0;i < 2;i++){
    azi[i] *= ONEEIGHTYOVERPI;
    fix_deg_angle((azi+i));
  }
}
