#include "fstrack.h"
#include <math.h>
/* 

subroutines analyzing timeseries and such

$Id: series_analyze.c,v 1.4 2010/12/31 18:59:17 becker Exp $

*/

/*

  calculate mean and standard deviation of x_i
  if hypoth is set, calc mean and stddev of sqrt(x_i^2 + y_i^2)
  if weighted is set, uses weights[n] for weighting the mean and so on

  this way of calculating the stddev is inaccurate but fast

*/
void calc_mean_and_stddev(float *x, float *y,
			  int n, COMP_PRECISION *mean,
			  COMP_PRECISION *stddev,
			  COMP_PRECISION *rms, 
			  my_boolean hypoth, my_boolean weighted,
			  float *weight)
{
  COMP_PRECISION sum1=0.0,sum2=0.0,tmp,ws;
  int i;
  if(n <= 1){
    fprintf(stderr,"calc_mean_and_stddev: error: n: %i\n",n);exit(-1);
  }
  ws=0.0;
  if(hypoth){// sqrt(x^2+y+2)
    for(i=0;i<n;i++){
      if(weighted){
	tmp = hypot(x[i],y[i]) * weight[i];
	ws += weight[i];
      }else{
	tmp = hypot(x[i],y[i]);
	ws += 1.0;
      }
      sum1 += tmp;sum2 += tmp * tmp;
    }
  }else{
    for(i=0;i<n;i++){
      if(weighted){
	tmp  = x[i] * weight[i];
	ws += weight[i];
      }else{
	tmp = x[i];
	ws += 1.0;
      }
      sum1 += tmp;sum2 += tmp*tmp;
    }
  }
  // standard deviation
  tmp = (ws * sum2 - sum1 * sum1) / (ws*(ws-1.0));
  if(tmp > 0)
    *stddev = sqrt(tmp);
  else
    *stddev = 0.0;
  *rms  = sqrt(sum2 / ws);// RMS 
  *mean = sum1 / ws;      // mean 
}


/* arithmetic mean */
COMP_PRECISION mean(COMP_PRECISION *x, int n)
{
  COMP_PRECISION sum=0.0;
  int i;
  for(i=0;i<n;i++)
    sum += x[i];
  if(n>0)
    return sum/(COMP_PRECISION)n;
  else
    return 0.0;
}
/* weighted arithmetic mean */
COMP_PRECISION wmean(COMP_PRECISION *x, COMP_PRECISION *w, int n)
{
  COMP_PRECISION sum = 0.0, sumw = 0.0;
  int i;
  for(i=0;i<n;i++){
    sum += x[i] * w[i];
    sumw += w[i];
  }
  if(n>0){
    return sum/(COMP_PRECISION)sumw;
  }else
    return 0.0;
}


/* standard deviation, quick and dirty */
COMP_PRECISION std(COMP_PRECISION *x, int n)
{
  int i;
  COMP_PRECISION sum1=0.0,sum2=0.0,ws,tmp;
  ws = (COMP_PRECISION)n;
  for(i=0;i<n;i++){
    sum1 += x[i];
    sum2 += x[i] * x[i];
  }
  // standard deviation
  tmp = (ws * sum2 - sum1 * sum1) / (ws*(ws-1.0));
  return sqrt(tmp);
}


/* 

statistic moments from numerical recipes, changed to C 0...n-1 style 
vectors

*/

void stat_moments(COMP_PRECISION *data,int n,
		  COMP_PRECISION *ave, /*  0: mean */
		  COMP_PRECISION *adev, /* 1: mean absolute deviation */
		  COMP_PRECISION *sdev, /* 2: standard deviation */
		  COMP_PRECISION *var, /*  3: variance */
		  COMP_PRECISION *skew, /* 4: skewness */
		  COMP_PRECISION *curt,/*  5: curtosis */
		  my_boolean report) 
     
{
  int j,nn;
  COMP_PRECISION ep=0.0,s,p;

  if (n <= 1) {
    fprintf(stderr,"stat_moments: n must be at least 1 in moments, n: %i\n",n);
    exit(-1);
  }
  s=0.0;
  nn=0;
  for (j=0;j<n;j++){
    if(finite(data[j])){
      s += data[j];
      nn++;
    }
  }
  /*  if(nn<=1){ */
  /*     fprintf(stderr,"stat_moments: non finite must be at least 2 in moments, nn: %i\n",nn); */
  /*     exit(-1); */
  /*   } */
  *ave=s/nn;
  *adev=(*var)=(*skew)=(*curt)=0.0;
  for (j=0;j<n;j++) {
    if(finite(data[j])){
      *adev += fabs(s=data[j]-(*ave));
      *var += (p=s*s);
      *skew += (p *= s);
      *curt += (p *= s);
    }
  }
  *adev /= nn;
  *var=(*var-ep*ep/nn)/(nn-1);
  *sdev=sqrt(*var);
  if (fabs(*var)>EPS_COMP_PREC) {
    *skew /= (n*(*var)*(*sdev));
    *curt=(*curt)/(n*(*var)*(*var))-3.0;
  } else{
    fprintf(stderr,"stat_moments: no skew/kurtosis when variance = %g (n: %i)\n",
	    *var,nn);
    *skew = *curt = 0.0;
  }
  if(report)
    fprintf(stderr,"moments: ave: %11g adev: %11g sdev: %11g skew: %11g curt: %11g nan: %i\n",
	    *ave,*adev,*sdev,*skew,*curt,n-nn);
}

/* 


fit harmonic terms assuming x has angles in radians

uses harm_fit_func as a fitting function and SVD LSQR routines to fit

npara = 2 + 2 * nharm



input: xrad[ndata] y[ndata]: x and y values of the data
       sigy[ndata]:         uncertainties in y[], only used if use_data_sigma
                             is TRUE, else pass as NULL
       use_svd:              if TRUE, use SVD; else use normal equations and Gauss - Jordan
       wmax_thres:           threshold for SVD singular value cut off
       nharm:                degree of harmonic function: 0: linear a + b x 1: a + b x + c sin(2x) + d cos(2x) ...

output: 
        chi2:  misfit chi square
        a[npara]: parameters
	siga[npara]: uncertainties 

we assume that we don't know the std of y

*/
void fit_harmonic(COMP_PRECISION *xrad, COMP_PRECISION *y, 
		  COMP_PRECISION *sigy,int ndata_orig, 
		  int nharm,COMP_PRECISION wmax_thres,
		  COMP_PRECISION *a, COMP_PRECISION *siga, 
		  COMP_PRECISION *chi2, 
		  my_boolean use_data_sigma,
		  my_boolean use_svd, 
		  my_boolean verbose)
{
  COMP_PRECISION *cov,sfac,*xloc=NULL,*yloc=NULL;
  int i,j,np,k,dof,ndata,nan;
  /* 
     
     take care of Nans

  */
  for(i=ndata=0;i < ndata_orig;i++){
    if(finite(xrad[i]) && finite(y[i])){
      ndata++;
      my_vecrealloc(&xloc,ndata,"fit_harmonic: xloc");
      my_vecrealloc(&yloc,ndata,"fit_harmonic: xloc");
      xloc[ndata-1] = xrad[i];
      yloc[ndata-1] = y[i];
    }
  }
  nan = ndata_orig - ndata;
  if(verbose && nan)
    fprintf(stderr,"fit_harmonic: excluding %i nans out of %i\n",
	    nan,ndata_orig);
  if((ndata < 3)||(nharm < 0)){
    fprintf(stderr,"fit_harmonic: parameter error: ndata: %i nharm: %i\n",
	    ndata,nharm);
    exit(-1);
  }
 
  np = 2 + 2 * nharm;
  /* 
     resize the parameter and uncertainty of parameter arrays 
  */
  if(!use_data_sigma){
    /* 
       we don't have actual uncertainties, need to fake them 
    */
    for(i=0;i<ndata;i++)
      sigy[i] = 1.0;		/* fake uncertainties */
  }
  /* covariance remains local */
  my_vecalloc(&cov,np * np,"fit_harmonic");	/* covariance */
  /* 
     perform the fit 
  */
  if(use_svd){
    svdfit_driver(xloc,yloc,sigy,ndata,np,wmax_thres,
		  (void (*)(void))(harm_fit_func),nharm,a,
		  cov,chi2);
  }else{
    lfit_driver(xloc,yloc,sigy,ndata,np,
		(void (*)(void))(harm_fit_func),nharm,a,
		cov,chi2);
  }
  if(!use_data_sigma){
    /* 
       rescale covariance matrix to obtain sigma in parameters
    */
    sfac = *chi2/(COMP_PRECISION)((ndata-2));
    for(i=0;i< np;i++)
      for(j=0;j< np;j++)
	cov[i*(np)+j] *= sfac;
  }
  /* 
     take sqrt of covariance matrix and obtain parameter uncertainties
  */
  for(i=0;i< np;i++){
    for(j=0;j< np;j++){
      k = i*np+j;
      cov[k] += MY_SIGN(cov[k]) * sqrt(fabs(cov[k])); 
    }
    siga[i] = cov[i*np+i]; 
  }
  dof = ndata - np - 2;

  if(verbose){
    /* 
       output of results 
    */
    if(use_svd)
      fprintf(stderr,"fit_harmonic: SVD fit %i harmonics to %i points: chi2: %11g rchi2: %8.3f (thres: %.2e)\n",
	      nharm,ndata,*chi2,sqrt(*chi2/dof),wmax_thres);
    else
      fprintf(stderr,"fit_harmonic: NE  fit %i harmonics to %i points: chi2: %11g rchi2: %8.3f\n",
	      nharm,ndata,*chi2,sqrt(*chi2/dof));
    
    fprintf(stderr,"fit_harmonic: a_00: %12g +/- %12g (%8.1f%%) (const.)\n",
	    a[0],siga[0],fabs(siga[0]/(a[0])*100));
    fprintf(stderr,"fit_harmonic: a_01: %12g +/- %12g (%8.1f%%) (linear)\n",
	    a[1],siga[1],fabs(siga[1]/(a[1])*100));
    for(i=1;i<=nharm;i++){
      fprintf(stderr,"fit_harmonic: a_%02i: %12g +/- %12g (%8.1f%%) (sin %2it)\n",
	      i*2,a[i*2],siga[i*2],fabs(siga[i*2]/(a[i*2])*100),
	      i);
      fprintf(stderr,"fit_harmonic: a_%02i: %12g +/- %12g (%8.1f%%) (cos %2it)\n",
	      i*2+1,a[i*2+1],siga[i*2+1],fabs(siga[i*2+1]/(a[i*2+1])*100),
	      i);
    }
    for(i=0;i<np;i++){
      fprintf(stderr,"fit_harmonic: |Cij| para %3i: ",i);
      for(j=0;j<=i;j++)
	fprintf(stderr,"%6.3f ",cov[i*np+j]/(siga[i]));
      fprintf(stderr,"\n");
    }
  }

  free(cov);
  free(xloc);
  free(yloc);
}

/* 
   subroutine to generate value of fitting function at x 

   f(x) = a_0 + a_1 *x + a_2 sin(x) + a_3 cos(x) + a_4 sin(2x) + a_5 cos(2x) .... to harmonic nharm*x

*/
void harm_fit_func(COMP_PRECISION x,COMP_PRECISION *y_base, 
		   int npara, int nharm)
{
  int i;
  COMP_PRECISION arg;
  if(npara != 2 + 2 * nharm){
    fprintf(stderr,"harm_fit_function: error: nharm: %i npara: %i\n",
	    nharm,npara);
    exit(-1);
  }
  y_base[0] = 1.0;		/* a_0 */
  y_base[1] = x;		/* a_1 */
  for(i=1;i<= nharm;i++){
    arg = (COMP_PRECISION)i * x;
    my_sincos(arg,(y_base+i*2),(y_base+i*2+1));
  }
  //fprintf(stderr,"x: %11g ",x);for(i=0;i<npara;i++)fprintf(stderr,"a%i: %12.4e ",i,y_base[i]);fprintf(stderr,"\n");
}
