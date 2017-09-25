/*


miscellaneous trigonometry routines dealing with different systems
(Cartesian and spherical) and other functions on a sphere


$Id: trig.c,v 1.18 2016/09/05 04:44:47 becker Exp $



*/
#include "fstrack.h"
#include <math.h>		/* just to be sure */
/* 



COORDINATE CONVERSION



*/
//
// go from spherical to cartesian coordinate system, input
// and output vectors have to be [3]
//
void rtp2xyz(COMP_PRECISION *xp,COMP_PRECISION *xc)
{
  COMP_PRECISION tmpdbl,sp,cp,ct,st;

  my_sincos(xp[FSTRACK_PHI],&sp,&cp);
  my_sincos(xp[FSTRACK_THETA],&st,&ct);

  tmpdbl=st * xp[FSTRACK_R];
  xc[FSTRACK_X]=tmpdbl * cp;
  xc[FSTRACK_Y]=tmpdbl * sp;
  xc[FSTRACK_Z]=ct * xp[FSTRACK_R];
}
/* 
   go from cartesian to spherical 
*/
void xyz2rtp(COMP_PRECISION *xc, COMP_PRECISION *xp)
{
  COMP_PRECISION tmp1,tmp2;
  tmp1  = xc[FSTRACK_X]*xc[FSTRACK_X]  + xc[FSTRACK_Y]*xc[FSTRACK_Y];
  tmp2  = tmp1 + xc[FSTRACK_Z]*xc[FSTRACK_Z];
  if(tmp2 <= 0.0)
    xp[FSTRACK_R] = 0.0;
  else
    xp[FSTRACK_R]     = sqrt(tmp2);
  xp[FSTRACK_THETA] = atan2(sqrt(tmp1),xc[FSTRACK_Z]);
  xp[FSTRACK_PHI]   = atan2(xc[FSTRACK_Y],xc[FSTRACK_X]);
}
/* 


VECTOR CONVERSION


*/
/*

  convert a polar vector at polar system location xp into a cartesian vector

*/
void polar_vec2cart_vec_at_xp(COMP_PRECISION *polar_vec,
			      COMP_PRECISION *cart_vec,
			      COMP_PRECISION *xp)
{
  COMP_PRECISION polar_basis[3][3];
  calc_polar_basis_at_r(polar_basis,xp); /* get the basis vectors in cartesian as matrix */
  abase_vec2bbase_vec(polar_vec,cart_vec,polar_basis);
}
/*
  convert a Cartesian vector at spherical system location xp into a
  spherical system vector
*/
void cart_vec2polar_vec_at_xp(COMP_PRECISION *cart_vec,
			      COMP_PRECISION *polar_vec,
			      COMP_PRECISION *xp)
{
  COMP_PRECISION cart_basis[3][3];
  calc_cart_basis_at_r(cart_basis,xp); /* get the basis vectors in cartesian as matrix */
  abase_vec2bbase_vec(cart_vec,polar_vec,cart_basis);
}
/* FORTRAN wrapper */
#ifdef GCC_USCR 
void polar_vec2cart_vec_at_xp__(COMP_PRECISION *a,COMP_PRECISION *b,
				COMP_PRECISION *xp)
#else
void polar_vec2cart_vec_at_xp_(COMP_PRECISION *a,COMP_PRECISION *b,
			       COMP_PRECISION *xp)
#endif
{polar_vec2cart_vec_at_xp(a,b,xp);}

/* 


MATRIX CONVERSIONS


*/
/* 
   
converts a a[9] matrix in polar system to a matrix b[9] in a cartesian
system, xp are the coordinates in the spherical system

*/
void polar_to_cart_mat_at_r(COMP_PRECISION *a,COMP_PRECISION *b,
			    COMP_PRECISION *xp)
{
  COMP_PRECISION f[3][3],fc[3][3];
  mat9to3x3(a,f);	/* convert a[9] matrix to [3][3] */
  polar_to_cart_mat_at_r3x3(f,fc,xp);
  mat3x3to9(fc,b);	/* reassign */
}
/* the other way round, convert cart to polar for a[9] matrices at spherical
   system location xp */
void cart_to_polar_mat_at_r(COMP_PRECISION *a,COMP_PRECISION *b,
			    COMP_PRECISION *xp)
{
  COMP_PRECISION f[3][3],fc[3][3];
  mat9to3x3(a,f);	/* convert a[9] matrix to [3][3] */
  cart_to_polar_mat_at_r3x3(f,fc,xp);
  mat3x3to9(fc,b);	/* reassign */
}
void rotate_cart_mat(COMP_PRECISION *a,COMP_PRECISION *b,
		     COMP_PRECISION alpha,COMP_PRECISION beta,
		     COMP_PRECISION gamma)
{
  COMP_PRECISION f[3][3],fc[3][3];
  mat9to3x3(a,f);	
  rotate_cart_mat_3x3(f,fc,alpha,beta,gamma);
  mat3x3to9(fc,b);	
}
/* 

   FORTRAN WRAPPERS. WARNING: THESE ROUTINES STILL OPERATE ON C-STYLE 
   MATRICES

*/
#ifdef GCC_USCR 
void polar_to_cart_mat_at_r__(COMP_PRECISION *a,COMP_PRECISION *b,
			      COMP_PRECISION *xp)
#else
void polar_to_cart_mat_at_r_(COMP_PRECISION *a,COMP_PRECISION *b,
			     COMP_PRECISION *xp)
#endif
{polar_to_cart_mat_at_r(a,b,xp);}
#ifdef GCC_USCR 
void cart_to_polar_mat_at_r__(COMP_PRECISION *a,COMP_PRECISION *b,
			     COMP_PRECISION *xp)
#else
void cart_to_polar_mat_at_r_(COMP_PRECISION *a,COMP_PRECISION *b,
			     COMP_PRECISION *xp)
#endif
{cart_to_polar_mat_at_r(a,b,xp);}
/* 
   
converts a a[3][3] matrix in polar system to a matrix b[3][3] in a cartesian
system, xp are the coordinates in the spherical system

*/
void polar_to_cart_mat_at_r3x3(COMP_PRECISION a[3][3],COMP_PRECISION b[3][3],
			       COMP_PRECISION *xp)
{
  COMP_PRECISION rot[3][3];
  /* obtain the rotation matrix at the location in spherical sys */
  calc_polar_basis_at_r(rot,xp); 
  rotate_mat(a,b,rot);	/* rotate tensor */
}
/* the other way around */
void cart_to_polar_mat_at_r3x3(COMP_PRECISION a[3][3],COMP_PRECISION b[3][3],
			       COMP_PRECISION *xp)
{
  COMP_PRECISION rot[3][3];
  /* obtain the rotation matrix at the location in spherical sys */
  calc_cart_basis_at_r(rot,xp); 
  rotate_mat(a,b,rot);	/* rotate tensor */
}
/* 

rotate a Caretsian system a matrix by angle alpha, beta, gamma angles
(radians) as in Dahlen and Tromp, p. 920


*/
void rotate_cart_mat_3x3(COMP_PRECISION a[3][3],COMP_PRECISION b[3][3],
			 COMP_PRECISION alpha, COMP_PRECISION beta,
			 COMP_PRECISION gamma)
{
  COMP_PRECISION rot[3][3];
  /* 
     obtain the rotation matrix at the location in spherical sys 
  */
  calc_rotmat_cart(rot,alpha,beta,gamma); 
  rotate_mat(a,b,rot);	/* rotate tensor */
}

/* 


MISCELLANEOUS FUNCTIONS


*/

/*
  
obtain a general rotation matrix with angles alpha, beta, and gamma
(in radians) as defined in Dahlen and Tromp, p. 921
  
*/
void calc_rotmat_cart(COMP_PRECISION r[3][3],COMP_PRECISION alpha,
		      COMP_PRECISION beta, COMP_PRECISION gamma)
{
  COMP_PRECISION sin_alpha,cos_alpha;
  COMP_PRECISION sin_beta,cos_beta;
  COMP_PRECISION sin_gamma,cos_gamma;

  my_sincos(alpha,&sin_alpha,&cos_alpha);
  my_sincos(beta, &sin_beta, &cos_beta);
  my_sincos(gamma,&sin_gamma,&cos_gamma);
  
  r[FSTRACK_X][FSTRACK_X] = cos_alpha*cos_beta*cos_gamma - sin_alpha*sin_gamma; 
  r[FSTRACK_X][FSTRACK_Y] = sin_alpha*cos_beta*cos_gamma + cos_alpha*sin_gamma;
  r[FSTRACK_X][FSTRACK_Z] = -sin_beta*cos_gamma;
  r[FSTRACK_Y][FSTRACK_X] = -cos_alpha*cos_beta*sin_gamma -sin_alpha*cos_gamma;
  r[FSTRACK_Y][FSTRACK_Y] = -sin_alpha*cos_beta*sin_gamma +cos_alpha*cos_gamma;
  r[FSTRACK_Y][FSTRACK_Z] = sin_beta*sin_gamma;
  r[FSTRACK_Z][FSTRACK_X] = cos_alpha*sin_beta;
  r[FSTRACK_Z][FSTRACK_Y] = sin_alpha*sin_beta;
  r[FSTRACK_Z][FSTRACK_Z] = cos_beta;
}
/* 

determine the polar basis vectors in cartesian space at position
r(r,theta,phi). this matrix has the r, theta, phi vectors in cart. components 
as rows R = ( e_r e_theta e_phi )

*/
void calc_polar_basis_at_r(COMP_PRECISION polar_basis[3][3],
			   COMP_PRECISION *r)
{
  COMP_PRECISION ct,cp,st,sp;

  my_sincos(r[FSTRACK_PHI],&sp,&cp);
  my_sincos(r[FSTRACK_THETA],&st,&ct);

  polar_basis[FSTRACK_X][FSTRACK_R]= st * cp;
  polar_basis[FSTRACK_X][FSTRACK_THETA]= ct * cp;
  polar_basis[FSTRACK_X][FSTRACK_PHI]= -sp;
  polar_basis[FSTRACK_Y][FSTRACK_R]= st * sp;
  polar_basis[FSTRACK_Y][FSTRACK_THETA]= ct * sp;
  polar_basis[FSTRACK_Y][FSTRACK_PHI]= cp;
  polar_basis[FSTRACK_Z][FSTRACK_R]= ct;
  polar_basis[FSTRACK_Z][FSTRACK_THETA]= -st;
  polar_basis[FSTRACK_Z][FSTRACK_PHI]= 0.0;
}
/* 

determine the cartesian basis vectors in polar space at position
r(r,theta,phi). 

*/
void calc_cart_basis_at_r(COMP_PRECISION cart_basis[3][3],
			  COMP_PRECISION *r)
{
  COMP_PRECISION ct,cp,st,sp;

  my_sincos(r[FSTRACK_PHI],&sp,&cp);
  my_sincos(r[FSTRACK_THETA],&st,&ct);

  cart_basis[FSTRACK_R][FSTRACK_X]= st * cp;
  cart_basis[FSTRACK_R][FSTRACK_Y]= st * sp;
  cart_basis[FSTRACK_R][FSTRACK_Z]= ct;
  cart_basis[FSTRACK_THETA][FSTRACK_X]= ct * cp;
  cart_basis[FSTRACK_THETA][FSTRACK_Y]= ct * sp;
  cart_basis[FSTRACK_THETA][FSTRACK_Z]= -st;
  cart_basis[FSTRACK_PHI][FSTRACK_X]= -sp;
  cart_basis[FSTRACK_PHI][FSTRACK_Y]=  cp;
  cart_basis[FSTRACK_PHI][FSTRACK_Z]= 0.0;
}
/* 
   
compute sines and cosines of x, given in radians

*/
void my_sincos(COMP_PRECISION x, COMP_PRECISION *sinv, COMP_PRECISION *cosv)
{
  *sinv = sin(x);
  *cosv = cos(x);
}
/* 
   
compute sines and cosines of x, given in degrees

*/
void my_sincosd(COMP_PRECISION x, 
		COMP_PRECISION *sinv, COMP_PRECISION *cosv)
{
  COMP_PRECISION xr;
  xr = DEG2RAD(x);
  *sinv = sin(xr);
  *cosv = cos(xr);
}
/* 
   rotate a 3x3 tensor using a general rotation matrix r 
   xout = r . xin . r^T
*/
void rotate_mat(COMP_PRECISION xin[3][3],
		COMP_PRECISION xout[3][3],
		COMP_PRECISION r[3][3])
{
  COMP_PRECISION tmp[3][3];
  int i,j,k;
  // calculate xin . r^T
  for(i=0;i<3;i++)
    for(j=0;j<3;j++){
      tmp[i][j] = 0.0;
      for(k=0;k<3;k++)
	tmp[i][j] += xin[i][k] * r[j][k];
    }
  // calculate r . tmp
  for(i=0;i<3;i++)
    for(j=0;j<3;j++){
      xout[i][j] = 0.0;
      for(k=0;k<3;k++)
	xout[i][j] += r[i][k] * tmp[k][j];
    }
}

/* 
   convert spherical coordinate system vector into 
   cartesian vector given a basis polar_basis, 

   abase_vec2bbase_vec(polar_vec,cart_vec,polar_basis)

   or

   convert Caretsian system vector into 
   spherical system vector given a basis cart_basis, i.e.

   abase_vec2bbase_vec(cart_vec,polar_vec,cart_basis)


*/
void abase_vec2bbase_vec(COMP_PRECISION *abase_vec,
			 COMP_PRECISION *bbase_vec,
			 COMP_PRECISION a_basis[3][3])
{
  calc_Ax3x3(a_basis,abase_vec,bbase_vec);
}

/* calculate matrix times vector
   y[3] = A[3][3] . x[3]  
*/
void calc_Ax3x3(COMP_PRECISION a[3][3], COMP_PRECISION *x,
		COMP_PRECISION *y)
{
  y[FSTRACK_X] = a[FSTRACK_X][FSTRACK_X] * x[FSTRACK_X] + a[FSTRACK_X][FSTRACK_Y] * x[FSTRACK_Y] + a[FSTRACK_X][FSTRACK_Z] * x[FSTRACK_Z];
  y[FSTRACK_Y] = a[FSTRACK_Y][FSTRACK_X] * x[FSTRACK_X] + a[FSTRACK_Y][FSTRACK_Y] * x[FSTRACK_Y] + a[FSTRACK_Y][FSTRACK_Z] * x[FSTRACK_Z];
  y[FSTRACK_Z] = a[FSTRACK_Z][FSTRACK_X] * x[FSTRACK_X] + a[FSTRACK_Z][FSTRACK_Y] * x[FSTRACK_Y] + a[FSTRACK_Z][FSTRACK_Z] * x[FSTRACK_Z];

}
// calculate x and y components from a unity vector with azimuth azi in degrees
void calc_xy_from_azi_deg(COMP_PRECISION azi, COMP_PRECISION *ax,COMP_PRECISION *ay)
{
  COMP_PRECISION tmp;
  tmp = azi * PIOVERONEEIGHTY;
  my_sincos(tmp,ax,ay);
}
// calculate azimuth in degrees from ax and ay
void calc_azi_from_axay(COMP_PRECISION ax, COMP_PRECISION ay, COMP_PRECISION *azi)
{
  *azi = atan2(ax,ay)*ONEEIGHTYOVERPI;
}
//
// similar to above but as function and with e_theta and e_phi
COMP_PRECISION azi_from_etep(COMP_PRECISION et, COMP_PRECISION ep)
{
  COMP_PRECISION azi;
  azi=atan2(ep,-et)*ONEEIGHTYOVERPI;
  fix_deg_angle(&azi);
  return azi;
}
// assuming unity vector:
// calculate dip (angle between vector and t-p plane, downward is negative
COMP_PRECISION dip_from_er_unityvec(COMP_PRECISION er)
{
  return  asin(er)*ONEEIGHTYOVERPI;
}

// adjust an angle that is given in degrees so that it is 0<=deg<=360
void fix_deg_angle(COMP_PRECISION *azi)
{
  while(*azi < 0)		
    (*azi) += 360.0;
  while(*azi > 360)
    (*azi) -= 360.0;
}
/*
  
  distance of two points on sphere in km
  x and y are vectors in r, theta, phi convention

*/
COMP_PRECISION dist_on_sphere(COMP_PRECISION *x,
			      COMP_PRECISION *y)
{
  COMP_PRECISION tmp1,tmp2,tmp3;
  tmp1 = (y[FSTRACK_THETA]-x[FSTRACK_THETA])/2.0;
  tmp1 = sin(tmp1);
  tmp1 = tmp1 * tmp1;
  tmp2 = (x[FSTRACK_PHI]-y[FSTRACK_PHI])/2.0;
  tmp2 = sin(tmp2);
  tmp2 = tmp2 * tmp2;
  tmp3 = tmp1;
  tmp3 += sin(x[FSTRACK_THETA]) * sin(y[FSTRACK_THETA]) * tmp2;
  if(fabs(tmp3) > EPS_PREC)
    tmp3 = sqrt(tmp3);
  else
    tmp3 = 0.0;
  return 2.0*asin(tmp3)*R_E;	/* i used asinl before, but that
				   didn't work on some compilers !? */
}
/* 
   compute lon (xpp[FSTRACK_X]), lat(xpp[FSTRACK_Y]) and depth z (xpp[FSTRACK_Z]) in km 
   from spherical coordinate vector xp 
*/
void lonlatz_from_xp(COMP_PRECISION *xp,COMP_PRECISION *xpp)
{
  xpp[FSTRACK_X] = PHI2LONGITUDE(xp[FSTRACK_PHI]);
  fix_deg_angle((xpp+FSTRACK_X));
  xpp[FSTRACK_Y] = THETA2LATITUDE(xp[FSTRACK_THETA]);
  xpp[FSTRACK_Z] = ZDEPTH(xp[FSTRACK_R]);
}

void xp_from_lonlatz(COMP_PRECISION lon, COMP_PRECISION lat, COMP_PRECISION z,
		     COMP_PRECISION *xp)
{
  xp[FSTRACK_PHI] =   LONGITUDE2PHI(lon);
  xp[FSTRACK_THETA] = LATITUDE2THETA(lat);
  xp[FSTRACK_R] =  ND_RADIUS(z);
}

/* 

substract the mean orientation from an oriented 
(180 degree periodic) vector whose whose azimuth [in degrees] is 
given in azi[]

returns the signed difference diff_azi (-90..90)
from the mean azimuthm mean_azi

*/
void diff_from_orient_mean(COMP_PRECISION *azi,int n,
			   COMP_PRECISION *diff_azi, 
			   COMP_PRECISION *mean_azi,
			   COMP_PRECISION *dt_weights,
			   int use_dt_for_avg,
			   COMP_PRECISION *mean_vec_dt)
{
  int i;
  COMP_PRECISION *x,*y,wsum =0.0,xm,ym;
  my_vecalloc(&x,n," diff_from_orient_mean");
  my_vecalloc(&y,n," diff_from_orient_mean");
  for(i=0;i < n;i++){		/* get the mean of the 2\theta */
    my_sincosd(2.0*azi[i],(x+i),(y+i));
  }
  if(use_dt_for_avg){		/* weighted by delay time */
    xm = wmean(x,dt_weights,n);
    ym = wmean(y,dt_weights,n);
  }else{
    xm = mean(x,n);
    ym = mean(y,n);
  }
  *mean_azi = atan2(xm,ym); 
  *mean_azi = RAD2DEG((*mean_azi)/2.0);	/* convert to degrees */
  fix_deg_angle(mean_azi);		/* move to 0...180 range */
  if(*mean_azi > 180)
    *mean_azi -= 180;
  /* amplitude */
  *mean_vec_dt = hypot(xm,ym);

  for(i=0;i<n;i++){		/* 
				   convert azimuths into deviations
				   from mean azimuth, accounting for
				   the 180 periodicity but allowing
				   for -/+ deviations
				*/
    diff_azi[i] = diff_orient_angle_azid(azi[i],*mean_azi);
  }
  free(x);free(y);
}



/* 

given the azimuth azi and mean azimuth mazi [both in degrees], compute
the difference in degrees, assuming that this is an oriented (ie. 180
periodic) quantity, allowing for positive and negative deviations

there's another version of this based on vectors (better)
*/
COMP_PRECISION diff_orient_angle_azid(COMP_PRECISION azi, 
				      COMP_PRECISION mazi)
{
  COMP_PRECISION da;
  /* make sure angles are in 0...360 */
  fix_deg_angle(&mazi);
  fix_deg_angle(&azi);
  da = azi - mazi;
  if(da > 0.0){
    if(da > 180.)
      da -= 180.;
    if(da > 90.)
      da -= 180.0;
  }else{
    if(da < -180.0)
      da += 180.0;
    if(da < -90.0)
      da += 180.0;
  }
  /* output is -90 ... 90  */
  return da;
}
COMP_PRECISION diff_orient_vector_azid(COMP_PRECISION *a, 
				       COMP_PRECISION *b,
				       my_boolean unity_vectors)
{
  double al,bl,cp;
  cp = vec2d_crossp(a,b);		
  if(!unity_vectors){
    al = hypot(a[0],a[1]);
    bl = hypot(b[0],b[1]);
    return (COMP_PRECISION)(RAD2DEG(asin(cp/(al*bl))));
  }else{
    return (COMP_PRECISION)(RAD2DEG(asin(cp)));
  }
}
