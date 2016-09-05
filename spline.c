#include "fstrack.h"

/* 

cubic spline interpolation routines, modified from numerical recipes


all routines changed to [0...n-1] convention

$Id: spline.c,v 1.2 2004/03/03 20:47:23 becker Exp twb $

*/

/* 

get the second derivatives y2[] for tabulated values y[n] and
locations x[n] with first derivatives yp1 and ypn and the endpoints. if those
are passed as >1e30, will pick natural conditions

(changed to 0 .. n-1 convention
*/
void nr_spline(COMP_PRECISION *x,COMP_PRECISION *y,int n,
	       COMP_PRECISION yp1,COMP_PRECISION ypn, 
	       COMP_PRECISION *y2)
{
  int i,k,nm1;
  COMP_PRECISION p,qn,sig,un,*u;
  nm1=n-1;
  my_vecalloc(&u,n,"spline: 1");
  if (yp1 > 0.99e30)		/* natural spline */
    y2[0]=u[0]=0.0;
  else {
    y2[0] = -0.5;
    u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
  }
  for (i=1;i < nm1;i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - 
      (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  if (ypn > 0.99e30)
    qn=un=0.0;
  else {
    qn=0.5;
    un=(3.0/(x[nm1]-x[nm1-1]))*
      (ypn-(y[nm1]-y[nm1-1])/(x[nm1]-x[nm1-1]));
  }
  y2[nm1]=(un-qn*u[nm1-1])/(qn*y2[nm1-1]+1.0);
  for (k=nm1-1;k>=0;k--)
    y2[k]=y2[k]*y2[k+1]+u[k];
  free(u);
}
/* 

given the xa[n] locations at which ya[n] was specified and the second 
derivatives as given by y2a and nr_spline, return the interpolated 
value at x, y

if calc_first_der is set, also return the first derivative

(chnaged to 0 ... n-1 convention)
*/
void nr_splint(COMP_PRECISION *xa,COMP_PRECISION *ya,
	       COMP_PRECISION *y2a,
	       int n,COMP_PRECISION x,COMP_PRECISION *y,
	       my_boolean calc_first_der, COMP_PRECISION *y1)
     
{
  int klo,khi,k,nm1;
  COMP_PRECISION h,b,a;
  nm1=n-1;
  
  klo=0;
  khi=nm1;
  /* 
     find indices by bisection, good for random x
  */
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    if (xa[k] > x) 
      khi = k;
    else 
      klo = k;
  }
  h = xa[khi] - xa[klo];
  if (fabs(h) < EPS_PREC) {
    PE("nr_splint: bad xa input ");
    exit(-1);
  }
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
  if(calc_first_der){		
    /* compute the first derivative */
    *y1  = (ya[khi]-ya[klo])/h;
    *y1 -= (3.0*a*a-1.0)/6.0 * h * y2a[klo];
    *y1 += (3.0*b*b-1.0)/6.0 * h * y2a[khi];
  }
}


