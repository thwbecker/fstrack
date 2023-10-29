#define TRUE (1)
#define FALSE (0)
#define DTOR (3.141592654/180.0)

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <rpc/rpc.h>
#include "SHahhead.h"
/* 


this code is based on bill menke's ah_cross_conv_1

as of Oct 2005


Comments by TWB (twb@usc.edu)

- I changed the I/O so that spectoseis output is read for radial and
transverse components, and to make things interoperable Vera' routines

$Id: ah_cross_conv_spectoseis.c,v 1.4 2005/10/22 01:37:43 becker Exp becker $ 


*/
void ah_cross_conv_spectoseis(float *,float *,float *,
			      int *,float *,float *,float *);

char *progname;

main(argc,argv)
     int argc;
     char *argv[];
{
  FILE *in;
  float *data1, *data2,*time;
  int i, n;
  float azi,rayp,delta,phi,Rmin;

  progname = argv[0];
  
  if( argc < 3 ) {
    fprintf(stderr,"usage: spectoseis_stream | %s rayp azi > results_text_file\n",argv[0]);
    fprintf(stderr,"\n");
		
    fprintf(stderr,"purpose: estimate simple one-layer anisotropy parameters\n");
    fprintf(stderr,"    from a single radially-polarized (e.g. SKS) phase\n");
    fprintf(stderr,"    using the Menke & Levin (2002) cross-convolution Method\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"this works in analogy to ah_cross_conv but reads spectoseis_stream output and takes\n");
    fprintf(stderr,"columns 1, 5, 6 as time, radial, and transverse component.\n\n");
    fprintf(stderr,"rayp and azi are just passed through, for consistency\n\n");
    
    fprintf(stderr,"Synopysis of method:\n");
    fprintf(stderr,"    When the rt flag is given and the radial horizontal component\n");
    fprintf(stderr,"    (ie. pointing back to the earthquake, A) and tangential horizontal\n");
    fprintf(stderr,"    (ie. 90 deg CCW from radial, B) component data are supplied, then we assume\n");
    fprintf(stderr,"    radial A(t) = s(t)*a(t) and tangential B(t) = s(t)*b(t) where\n");
    fprintf(stderr,"    where   s(t) = unknown source wavelet, and a(t), b(t) are the\n");
    fprintf(stderr,"    respose fuctions for a single anisotropic layer at normal incidence\n");
    fprintf(stderr,"            a(t) = a0 delta(t) + a1 delta(t-t1) and\n");
    fprintf(stderr,"            b(t) = b0 delta(t) + b1 delta(t-t1)\n");
    fprintf(stderr,"            a0 = cos(theta)**2; a1 = sin(theta)**2 and\n");
    fprintf(stderr,"            b0 = cos(theta)*sin(theta); b1 = -b0 and\n");
    fprintf(stderr,"            theta = CCW angle from radial to fast axis\n");
    fprintf(stderr,"    This program determines  t1 and theta by a grid search by\n");
    fprintf(stderr,"    minimizing the amplitude-normalized squared L2 norm, R, of the cross-convolution,\n");
    fprintf(stderr,"            e(t) = A(t)*b(t) - B(t)*a(t).\n");
    fprintf(stderr,"    Here R = ||e||**2 / [ ||A(t)*b(t)||**2 + ||B(t)*a(t)||**2 ]\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"    When the ne flag is given and the north and east horizontal\n");
    fprintf(stderr,"    components are supplied, then all calculateion are performed\n");
    fprintf(stderr,"    in the north/east coordinate system and e(t) is redefined in terms\n");
    fprintf(stderr,"    of north and east components\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"output_ah_file, ah file with best-fitting\n");
    fprintf(stderr,"    A(t)*b(t)\n");
    fprintf(stderr,"    B(t)*a(t)\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"error_surface_text_file, the backazimuth followed by a table of\n");
    fprintf(stderr,"    phi, azimuth of fast axis in deg east of north\n");
    fprintf(stderr,"    t1, anisotropic delay time in seconds\n");
    fprintf(stderr,"    R, value of amplitude-normalized squared L2 norm\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"results_text_file, a table of:\n");
    fprintf(stderr,"    backazimuth of earthquake in deg east of north\n");
    fprintf(stderr,"    phi, best-fitting azimuth of fast axis in deg east of north\n");
    fprintf(stderr,"    t1, best-fitting anisotropic delay time in seconds\n");
    fprintf(stderr,"    R, minimum value of amplitude-normalized squared L2 norm\n");
    fprintf(stderr,"\n");
    exit(-1);
  }
  
  
  sscanf(argv[1],"%f",&rayp);
  sscanf(argv[2],"%f",&azi);
  in = stdin;


  n = 1;
  if( (time=(float*)malloc(n*sizeof(float))) == NULL ) {
    fprintf(stderr,"error: %s: out of memory!\n", argv[0]);
    exit(-1);
  }
  if( (data1=(float*)malloc(n*sizeof(float))) == NULL ) {
    fprintf(stderr,"error: %s: out of memory!\n", argv[0]);
    exit(-1);
  }
  if( (data2=(float*)malloc(n*sizeof(float))) == NULL ) {
    fprintf(stderr,"error: %s: out of memory!\n", argv[0]);
    exit(-1);
  }
  /* read in radial and transverse component */
  n=0;
  while( fscanf(in,"%f %*f %*f %*f %f %f %*f %*f %*f %*f\n",
		(time+n),(data1+n),(data2+n))==3){
    if( (time=(float*)realloc(time,(n+2)*sizeof(float))) == NULL ) {
      fprintf(stderr,"error: %s: out of memory!\n", argv[0]);
      exit(-1);
    }
    if( (data1=(float*)realloc(data1,(n+2)*sizeof(float))) == NULL ) {
      fprintf(stderr,"error: %s: out of memory!\n", argv[0]);
      exit(-1);
    }
    if( (data2=(float*)realloc(data2,(n+2)*sizeof(float))) == NULL ) {
      fprintf(stderr,"error: %s: out of memory!\n", argv[0]);
      exit(-1);
    }
   
    n++;
  } /* done data readind */


  ah_cross_conv_spectoseis(time,data1,data2,&n,&phi,&delta,&Rmin);
  
  if(Rmin > 0.5){		/* no splitting */
    fprintf(stdout,"%f %f nan 0.0 %f\n", rayp, azi, Rmin);
  }else{
    fprintf(stdout,"%f %f %f %f %f\n", rayp, azi, phi, delta,Rmin);
  }

  fflush(stdout);
  
    
 
  free(data1); free(data2); 
}


  
