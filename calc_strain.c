#include "fstrack.h"
/*

 $Id: calc_strain.c,v 1.23 2007/02/28 19:12:44 becker Exp $
 
 routines to determine various train matrices and the 
 bailout criterion accoring to mstrain

*/

/* 

given a velocity gradient matrix, compute the strain rate and vorticity,
if 

*/
void calc_strainrate_and_vorticity(COMP_PRECISION *vgm,
				   COMP_PRECISION *strainrate,
				   COMP_PRECISION *vorticity,
				   my_boolean use_units)
{
  COMP_PRECISION e[9],v[9],eval[3];
  /* compute symmetric part of vgm, strain-rate tensor */
  calc_cd_symm_part(vgm,e);
  /* get eigenvalues  */
  char_strain_rate_from_e(e,eval);
  //fprintf(stderr,"e: %g %g %g %g %g %g %g %g %g\n",e[0],e[1],e[2],e[3],e[4],e[5],e[6],e[7],e[8]);
  *strainrate = (eval[2] - eval[0])/2.0;	/* this should be a mean of
					   two very similar numbers e1 = -e3, e2 = 0*/
  //fprintf(stderr,"s: %11g %11g %11g %11g\n",eval[0],eval[1],eval[2],*strainrate);
  /* compute assymmetric part */
  calc_cd_asymm_part(vgm,v);
  char_strain_rate_from_e(v,eval); /* eigenvalues */
  //fprintf(stderr,"v: %g %g %g %g %g %g %g %g %g\n",v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8]);
  *vorticity = (eval[2]-eval[0])/2.0; /*  */
  //fprintf(stderr,"v: %11g %11g %11g %11g\n",eval[0],eval[1],eval[2],*vorticity);

  if(use_units){
    *strainrate /= TIMESCALE_S;
    *vorticity  /= TIMESCALE_S;
  }
}

/* 
   wrapper for the routine below, compute char strain-rate 
   from the symmetric part of the velocity gradient matrix
   vg[9] is non-symmetric and input in C fashion 
   (doesn't matter here)

*/
COMP_PRECISION char_strain_rate_from_vg(COMP_PRECISION *vg)
{
  COMP_PRECISION strain[9],val[3];
  //
  // calculate the symmetric part of vg (vel grad matrix)
  // to obtain the strain rate matrix
  //
  calc_cd_symm_part(vg,strain);
  return char_strain_rate_from_e(strain,val);
}
/* 
   
compute the characteristic strain-rate, defined as the maximum 
of the absolute values of the eigenvalues of the symmetric
strain-rate matrix, which is input. also returns the eigenvalues 
of the strain matrix, dimension those as [3]


*/
COMP_PRECISION char_strain_rate_from_e(COMP_PRECISION *strain, COMP_PRECISION *eval)
{
  COMP_PRECISION vec[9];	/* vec is not actually assigned! */
  //
  // calculate the eigenvalues of the strain rate matrix
  //
  calc_eigensystem_sym(strain,eval,vec,FALSE);
  //
  // the largest strain rate will define the shortest timescale
  //
  // the characteristic strain rate is the max(abs(EVs(strain)))
  return max_abs_vec(eval,3);
}


/*
 calculate Lagrangian strain tensor E^L = .5* (F^T dot F - I) 
 eg. Dahlen and Tromp 2.20

 not good for finite strain comparison
 input is F, output is e, both are 3 x 3

*/
void calc_el_strain(COMP_PRECISION *f, COMP_PRECISION *e)
{
  calc_at_dot_a(f,e);
  e[0] -= 1.0;e[4] -= 1.0;e[8] -= 1.0;
  scale_vector(e,0.5,9);
}
void calc_el_strain_(COMP_PRECISION *f, COMP_PRECISION *e)
{
  calc_el_strain(f,e);
}
/*

  calculate the Cauchy deformation tensor
  B^-1 = (F^-1)^T (F^-1), B (not B^-1) is identical to L2
  (Malvern, 1969)

*/
void calc_cauchy_strain(COMP_PRECISION *f, COMP_PRECISION *b)
{
  COMP_PRECISION finv[9];
  invert3x3c(f,finv);
  calc_at_dot_a(finv,b);
}
/*
  
  calculate the squared left stretch tensor L: L^2 = F dot F^T
  e.g. Dahlen and Tromp, eq 2.28
  
  the eigenvaectors and values of L form the strain ellipsoid at
  the final location of the tracer path

*/
void calc_l2_left_stretch_(COMP_PRECISION *f, COMP_PRECISION *e)
{
  calc_l2_left_stretch(f, e);
}
void calc_l2_left_stretch(COMP_PRECISION *f, COMP_PRECISION *e)
{
  calc_a_dot_at(f,e);
}
/*

  calculate the max stretch 

  input is deformation matrix F, uses mstrain routine for the strain
  criteria 

  also calculates the squared eigenvalues of L, val

*/
COMP_PRECISION max_strain_from_def(COMP_PRECISION *f, COMP_PRECISION *val)
{
  COMP_PRECISION left_stretch[9];
  calc_l2_left_stretch(f,left_stretch);// compute left-stretch L^2 strain from F
  return max_strain_from_left_stretch(left_stretch,val);
}
/* compute the max strain from the left stretch^2 matrix and compute 
   the eigenvalues */
COMP_PRECISION max_strain_from_left_stretch(COMP_PRECISION *l2, 
					    COMP_PRECISION *val)
{
  COMP_PRECISION vec[9];
  calc_eigensystem_sym(l2,val,vec,FALSE);// get eigenvalues of L^2 strain
  return mstrain(val);
}
/*

  returns max strain given that the eigenvalues of B in val are sorted
  in ascending order, ie. val[0] = e3 and val[2] = e1
  
*/
COMP_PRECISION mstrain(COMP_PRECISION *val)
{
#ifdef MAX_EVAL_STRAIN
  // simply use size of largest eigenvalue
  return sqrt(val[2])-1.0;
#elif defined MAX_FRAC_EVAL_STRAIN
  COMP_PRECISION f1,f2,f3;
  //
  // use logarithmic strains
  //
  // \xi = log(l_1/l_2) and 
  //
  // and
  //
  // \chi = log(l_2/l_3)
  //
  // since we are passing the eigenvalues of B, l^2
  // use log(l_1^2/l_2^2) = 2 log(l_1/l_2)
  //
  
  f1=0.5*log(val[2]/val[1]);// l_1/l_2
  f2=0.5*log(val[1]/val[0]);// l_2/l_3
  return MAX(f1,f2);
  //f3=0.5*log(val[2]/val[0]);// l_1/l_3
  //return f3;

#else
  PE("mstrain: error: no criterion defined during compilation, see fstrack.h");
  exit(-1);
#endif
}



