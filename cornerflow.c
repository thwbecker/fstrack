#include "fstrack.h"

/*
  refer to McKenzie (GJRAS, 58, 689, 1979) Table 1 for the 
  definitions of the parameters
  we refer to vr or vx as vel, and modes 0,1,2,3 refer 
  to BC types (a), (b), (c), (d)
  by default, mode=0, vel=1, alpha=90
  
  typos have been corrected
  


  $Id: cornerflow.c,v 1.3 2010/12/31 18:59:17 becker Exp $
  

*/
//
// obtain constants
//
//
void  calc_cornerflow_constants(int mode, COMP_PRECISION vel, 
				COMP_PRECISION alpha, 
				COMP_PRECISION *a, COMP_PRECISION *b, 
				COMP_PRECISION *c, COMP_PRECISION *d,
				char **argv)
{
  COMP_PRECISION cosa,sina;
  cosa = cos(alpha);
  sina = sin(alpha);
  switch(mode){
    //
    // listed are the BCs for boundary one and two as in Table 1
    // NOTE: the D coefficient for (c) should be -vr..., there's a typo
    // in the original Table which has been corrected below
    //
  case 0:// s_rt = 0; v = vx
    *a = (2.0 * vel * cosa*cosa)/   (2.0 * alpha - sin(2.0*alpha));
    *b = 0.0;
    *c = 0.0;
    *d= -(2.0 * vel)/               (2.0 * alpha - sin(2.0*alpha));
    break;
  case 1:// s_rt = 0; v = vr;
    fprintf(stderr,"%s: cornerflow: mode 1 (b) not working yet\n",argv[0]);exit(-1);
    *a = (2.0 * vel * alpha * cosa)/(2.0 * alpha - sin(2.0*alpha));
    *b = 0.0;
    *c = 0.0;
    *d= -(2.0 * vel )/              (2.0 * alpha - sin(2.0*alpha));
    break;
  case 2:// v = 0; v = vr;
    *a = (vel * alpha * sina)/(alpha*alpha - sina*sina);
    *b = 0.0;
    *c = (vel*(alpha * cosa - sina))/(alpha * alpha - sina * sina);
    *d = (vel * alpha * sina) /      (alpha * alpha - sina * sina);
    break;
  case 3:// v = vr; v = -vr;
    *a = (vel * alpha)/(alpha + sina);
    *b = 0.0;
    *c = (-vel*(1.0+ cosa))/(alpha + sina);
    *d = (vel * sina)      /(alpha + sina);
    break;
  default:
    fprintf(stderr,"%s: cornerflow: error: mode %i undefined\n",argv[0],mode);
    exit(-1);
    break;
  }
  if(fabs(*a)<EPS_COMP_PREC) *a = 0.0;
  if(fabs(*b)<EPS_COMP_PREC) *b = 0.0;
  if(fabs(*c)<EPS_COMP_PREC) *c = 0.0;
  if(fabs(*d)<EPS_COMP_PREC) *d = 0.0;
}
//
// obtain stream function at xcyl given a,b,c,d
COMP_PRECISION stream(COMP_PRECISION *xcyl, 
		      COMP_PRECISION a, COMP_PRECISION b, COMP_PRECISION c, 
		      COMP_PRECISION d)
{
  COMP_PRECISION cost,sint;
  cost = cos(xcyl[FSTRACK_THETA]);
  sint = sin(xcyl[FSTRACK_THETA]);
  return (xcyl[FSTRACK_R] * (a * sint + b * cost + 
		     c * xcyl[FSTRACK_THETA] * sint + 
		     d * xcyl[FSTRACK_THETA] * cost));
}
//
// obtain velocity at cylidrical coordinates xcyl given constants a,b,c,d
// returns veccyl
//
void cylvel(COMP_PRECISION *xcyl, COMP_PRECISION *veccyl,
	    COMP_PRECISION a,COMP_PRECISION b, COMP_PRECISION c, 
	    COMP_PRECISION d)
{
  COMP_PRECISION cost,sint;
  cost = cos(xcyl[FSTRACK_THETA]);sint = sin(xcyl[FSTRACK_THETA]);
  veccyl[FSTRACK_R] = a * cost - b * sint + c * (sint + xcyl[FSTRACK_THETA] * cost) + 
    d * (cost - xcyl[FSTRACK_THETA] * sint);
  veccyl[FSTRACK_THETA] = 
    - a * sint - b * cost - c * xcyl[FSTRACK_THETA] * sint - d * xcyl[FSTRACK_THETA] * cost;
  veccyl[FSTRACK_R]     *= xcyl[FSTRACK_R];
  veccyl[FSTRACK_THETA] *= xcyl[FSTRACK_R];
}

// convert from cylindrical to cartesian coordiates
void cyl2cart(COMP_PRECISION *xcyl, COMP_PRECISION *xcart)
{
  xcart[FSTRACK_X] =  sin(xcyl[FSTRACK_THETA]) * xcyl[FSTRACK_R];
  xcart[FSTRACK_Y] = -cos(xcyl[FSTRACK_THETA]) * xcyl[FSTRACK_R];
}
// convert from cartesian to cylindrical coordinates
void cart2cyl(COMP_PRECISION *xcart, COMP_PRECISION *xcyl)
{
  xcyl[FSTRACK_THETA] = atan2(xcart[FSTRACK_X],-xcart[FSTRACK_Y]);
  xcyl[FSTRACK_R] =     hypot(xcart[FSTRACK_X],xcart[FSTRACK_Y]);
}
// convert from cylindrical vector vcyl at cylindrical coords xcyl to
// cartesian vector
void cylvec2cartvec(COMP_PRECISION *xcyl, COMP_PRECISION *vcyl, 
		    COMP_PRECISION *vcart)
{
  COMP_PRECISION sint,cost;
  cost = cos(xcyl[FSTRACK_THETA]);sint = sin(xcyl[FSTRACK_THETA]);
  vcart[FSTRACK_X] =  sint * vcyl[FSTRACK_R] + cost * vcyl[FSTRACK_THETA];
  vcart[FSTRACK_Y] = -cost * vcyl[FSTRACK_R] + sint * vcyl[FSTRACK_THETA];
}
