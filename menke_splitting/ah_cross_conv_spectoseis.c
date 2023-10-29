#define TRUE (1)
#define FALSE (0)
#define DTOR (3.141592654/180.0)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <rpc/rpc.h>
#include "SHahhead.h"
/* 


this code is based on bill menke's ah_cross_conv_1

as of Oct 2005


Comments by TWB (twb@usc.edu)

- I changed the I/O so that spectoseis output is read for radial and
transverse components, and to make things interoperable Vera' routines

$Id: ah_cross_conv_spectoseis.c,v 1.4 2005/10/22 01:37:43 becker Exp becker $ 


this is the working routine, driver is in ah_cross_conv_spectoseis_driver

*/
void ah_cross_conv_spectoseis(float *,float *,float *,
			      int *,float *,float *,float *);
void ah_cross_conv_spectoseis_(float *,float *,float *,
			       int *,float *,float *,float *);

void ah_cross_conv_spectoseis_(float *time,float *data1,float *data2,
			       int *nout,float *fastaz,float *delta,float *misfit)
{
  ah_cross_conv_spectoseis(time,data1,data2,nout,fastaz,delta,misfit);
}

void ah_cross_conv_spectoseis(float *time,float *data1,float *data2,
			      int *nout,float *fastaz,float *delta,float *misfit)
{
  float *data3,*data4,*data5,*data6,*data7;
  double  E, e, Emin, tlag;
  int i, j, ii, jj, ierror, lag1, lag1p, lagmin, maxlag, itheta, first,n;
  double ea, eb, Ea, Eb, Eamin, Ebmin, R, Rmin, theta, thetamin;
  double a0, a1, b0, b1, a0p, a1p, b0p, b1p;
  double a0min, a1min, b0min, b1min;
  float azi,rayp;
  double phi, ctheta, stheta;
 

  n = *nout;			/* for fortran interop */

  if( (data3=(float*)malloc(n*sizeof(float))) == NULL ) {
    fprintf(stderr,"error: %s: out of memory!\n", "ah_cross_conv_spectoseis");
    exit(-1);
  }
  if( (data4=(float*)malloc(n*sizeof(float))) == NULL ) {
    fprintf(stderr,"error: %s: out of memory!\n", "ah_cross_conv_spectoseis");
    exit(-1);
  }
  if( (data5=(float*)malloc(n*sizeof(float))) == NULL ) {
    fprintf(stderr,"error: %s: out of memory!\n", "ah_cross_conv_spectoseis");
    exit(-1);
  }
  if( (data6=(float*)malloc(n*sizeof(float))) == NULL ) {
    fprintf(stderr,"error: %s: out of memory!\n", "ah_cross_conv_spectoseis");
    exit(-1);
  }
  if( (data7=(float*)malloc(n*sizeof(float))) == NULL ) {
    fprintf(stderr,"error: %s: out of memory!\n", "ah_cross_conv_spectoseis");
    exit(-1);
  }
  /* 
       max time lag 
  */
  
  tlag = 3.5;

  /* find sampling rate */
  for(i=0;i<n;i++)
   if(i==1){
      *delta = time[1] - time[0];
   }else if(i > 1){
     if(fabs((time[i] - time[i-1])- (*delta))>1e-4){
       fprintf(stderr,"%s: error: give even time sampling: %g %g\n",
	       "ah_cross_conv_spectoseis",*delta,time[i]-time[i-1]);
       exit(-1);
     }
   }
  

  maxlag = (int) (0.5+tlag/(*delta));
  
  first=TRUE;
  
  for( itheta=0; itheta<180; itheta++ ) { 
    
    theta = (double) itheta;
    
    ctheta=cos(DTOR*theta); stheta=sin(DTOR*theta);
    a0p = ctheta*ctheta;
    a1p = stheta*stheta;
    b0p = stheta*ctheta;
    b1p = (-b0p);
    /* R/T system, don't need bazi */
    a0=a0p; a1=a1p; b0=b0p; b1=b1p;
    
    
    for( lag1p=1; lag1p<=maxlag; lag1p++ ) {
      
      lag1=(lag1p);
      
      E = 0.0; Ea=0.0; Eb=0.0;
      for( ii=0; ii<n; ii++ ) {
	jj = ii-lag1;
	if( jj<0 ) jj+=n; /* circular */
	e =  b0 * data1[ii];
	e += a0 * (-data2[ii]);
	ea = b0 * data1[ii];
	eb = a0 * data2[ii];
	e  += b1 * data1[jj] + a1 * (-data2[jj]);
	ea += b1 * data1[jj];
	eb += a1 * data2[jj];
	data3[ii]=(float)ea;
	data4[ii]=(float)eb;
	data5[ii]=(float)e;
	E += e*e; Ea += ea*ea; Eb += eb*eb;
      }
	
      R = E/(Ea+Eb);
      
      if( first ) {
	a0min=a0; a1min=a1; b0min=b0; b1min=b1;
	Emin=E; Eamin=Ea; Ebmin=Eb; Rmin=R;
	lagmin=lag1; thetamin=theta; first=FALSE;
	for( ii=0; ii<n; ii++ ) {
	  data6[ii]=data3[ii];
	  data7[ii]=data4[ii];
	}
      }
      else if( R<Rmin ) {
	a0min=a0; a1min=a1; b0min=b0; b1min=b1;
	Emin=E; Eamin=Ea; Ebmin=Eb; Rmin=R;
	lagmin=lag1; thetamin=theta;
	for( ii=0; ii<n; ii++ ) {
	  data6[ii]=data3[ii];
	  data7[ii]=data4[ii];
	}
      }
      
      phi = -theta;
      /* azimuth to 0-360 range */
      if( phi<0.0 ) phi+=360.0;
      else if( phi>=360.0 ) phi-=360.0;
      /* fast axis has 2-theta variation */
      if( phi>180.0 ) phi-=180.0;
      // this would be the lag scan output
      //fprintf(f4,"%f\t%f\t%f\n", phi, *delta*(double)lag1, R );
      
    } /* end lag loop */
  } /* end theta loop */
    
  
  phi = -thetamin;
  /* azimuth to 0-360 range */
  if( phi<0.0 ) phi+=360.0;
  else if( phi>=360.0 ) phi-=360.0;
  /* fast axis has 2-theta variation */
  if( phi>180.0 ) phi-=180.0;
  

  *misfit = Rmin;
  *fastaz = phi;
  *delta *= (double)lagmin;
  
    
 
  free(data3); free(data4);
  free(data5); free(data6); free(data7);
}


  
