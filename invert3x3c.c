#include "fstrack.h"
/*

  invert a 3 by 3 matrix that is given in vector form
  
  0 1 2      [0][0] [0][1] [0][2]
  3 4 5  or  [1][0] [1][1] [1][2]
  6 7 8      [2][0] [2][1] [2][2]

  this routine has the equivalent functionality to invert3x3 which 
  is in FORTRAN. However, it is five times faster

  $Id: invert3x3c.c,v 1.3 2004/04/19 18:41:26 becker Exp $

*/
void invert3x3c(COMP_PRECISION *a, COMP_PRECISION *ainv)
{
  COMP_PRECISION          m00 = a[4] * a[8] - a[7] * a[5];
  COMP_PRECISION          m01 = a[5] * a[6] - a[8] * a[3];
  COMP_PRECISION          m02 = a[3] * a[7] - a[6] * a[4];
  // compute determinant
  register COMP_PRECISION d = a[0] * m00 + a[1] * m01 + a[2] * m02;
  if(fabs(d) < EPS_PREC){
    fprintf(stderr,"invert3x3c: WARNING: matrix (near) singular: det: %g\n",
	    d);
    print_vector(a,9,stderr);
  }
  ainv[0] = m00 / d;
  ainv[1] = (a[7] * a[2] - a[1] * a[8]) / d;
  ainv[2] = (a[1] * a[5] - a[4] * a[2]) / d;
  ainv[3] = m01 / d;
  ainv[4] = (a[8] * a[0] - a[2] * a[6]) / d;
  ainv[5] = (a[2] * a[3] - a[5] * a[0]) / d;
  ainv[6] = m02 / d;
  ainv[7] = (a[6] * a[1] - a[0] * a[7]) / d;
  ainv[8] = (a[0] * a[4] - a[3] * a[1]) / d;
}
// FORTRAN wrapper
void invert3x3c_(COMP_PRECISION *a, COMP_PRECISION *ainv)
{
  invert3x3c(a,ainv);
}
