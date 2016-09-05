/* 

numerical recipes implementation of SVD and Gauss-Jordan direct matrix
method fitting routines


$Id: datafit_util.c,v 1.4 2011/04/12 06:19:09 becker Exp $

*/
#include "precision.h"

#include "nr_defines.h"
/* 

main SVD fitting routine

input:

data: x[ndata] locations, y[ndata] values with sigy[ndata] uncertainties

fit function: of form
    (*funcs)(x[i],afunc,ma,func_para1): which will evaluate the ma basis functions
                                        at location x given the void parameter func_para1

output:

parameters: a[npara] estimates and cov[npara][npara] covariance

performs SVD fit using the supplied function and computes the best fit parameters a[npara]


THIS ROUTINE WILL WORK WHEN CALLED THE REGULAR C WAY, USING
X[0....N-1] ADDRESSING. LIKEWISE, THE FITTING FUNCTION WILL BE
EXPECTED TO BE IN THE NORMAL C FASHION. THIS MEANS THAT CALLS BELOW
HAVE BEEN MODIFIED SO AS TO TAKE THE NUMERICAL RECIPES [1...N] CALLING
WAY INTO ACCOUNT


*/
void svdfit_driver(COMP_PRECISION *x, /* input: x values of data [ndata] */
		   COMP_PRECISION *y, /* y values of data  [ndata] */
		   COMP_PRECISION *sigy, /* standard deviations [ndata] */
		   int ndata,	/* number of data points */
		   int npara, 	/* number  of parameters */
		   COMP_PRECISION wmax_thres, /* threshold for SVD fitting */
		   void (*funcs)(),
		   int func_para1, /* parameter for the fit function */
		   COMP_PRECISION *a, /* output: parameters a[npara] */
		   COMP_PRECISION *cov, /* covariance matrix [npara][npara] */
		   COMP_PRECISION *chisq /* misfit */)
{
  COMP_PRECISION **u,**v,**nr_cov,*w;
  int i,j;

  /* allocate work space */
  u = nr_matrix(1,(long)ndata,1,(long)npara);
  v = nr_matrix(1,(long)npara,1,(long)npara);
  w = nr_vector(1,(long)npara);
  nr_cov = nr_matrix(1,(long)npara,1,(long)npara);
  /* 
     call the numerical recipes routines, with the necessary shift for
     the external vectors 
  */
  nr_svdfit((x-1),(y-1),(sigy-1),ndata,(a-1),npara,u,v,w,chisq,
	    wmax_thres,funcs,func_para1);
  /* 
     get covariance matrix 
  */
  nr_svdvar(v,npara,w,nr_cov);
  for(i=0;i<npara;i++)
    for(j=0;j<npara;j++)
      cov[i*npara+j] = nr_cov[i+1][j+1];
  /* free temporary arrays */
  free(w);
  nr_free_matrix(nr_cov,1,npara,1,npara);
  nr_free_matrix(v,1,npara,1,npara);
  nr_free_matrix(u,1,ndata,1,npara);
}

/* 

same as above: driver routine for numerical recipes lfit Gauss-Jordan directe
LSQR solver

THIS ROUTINE WILL WORK WHEN CALLED THE REGULAR C WAY, USING
X[0....N-1] ADDRESSING. LIKEWISE, THE FITTING FUNCTION WILL BE
EXPECTED TO BE IN THE NORMAL C FASHION. THIS MEANS THAT CALLS BELOW
HAVE BEEN MODIFIED SO AS TO TAKE THE NUMERICAL RECIPES [1...N] CALLING
WAY INTO ACCOUNT


*/
int lfit_driver(COMP_PRECISION *x, /* input: x values of data [ndata] */
		COMP_PRECISION *y, /* y values of data  [ndata] */
		COMP_PRECISION *sigy, /* standard deviations [ndata] */
		int ndata,	/* number of data points */
		int npara, 	/* number  of parameters */
		void (*funcs)(),
		int func_para1, /* parameter for the fit function */
		COMP_PRECISION *a, /* output: parameters a[npara] */
		COMP_PRECISION *cov, /* covariance matrix [npara][npara] */
		COMP_PRECISION *chisq /* misfit */)
{
  COMP_PRECISION **nr_cov;
  int *ia;
  int i,j,err;
  /* allocate work space */
  nr_cov = nr_matrix(1,npara,1,npara);
  my_ivecalloc(&ia,npara,"lfit_driver");
  /* 
     the ia array holds the inversion flags: 1: solve for parameter 0: hold fixed
     by default, we invert for all parameters
  */
  for(i=0;i<npara;i++)
    ia[i] = 1;		
  /* 
     call the numerical recipes routines, with the necessary shift for
     the external vectors
  */
  err = nr_lfit((x-1),(y-1),(sigy-1),ndata,(a-1),(ia-1),
		npara,nr_cov,
		chisq,funcs,func_para1);
  if(err != 0){
    /* 
       resort covariance matrix 
    */
    for(i=0;i<npara;i++)
      for(j=0;j<npara;j++){
	cov[i*npara+j] = nr_cov[i+1][j+1];
      }
  }
  /* free temporary arrays */
  free(ia);
  nr_free_matrix(nr_cov,1,npara,1,npara);
  return err;
}


/* 


fit polynomials to f(x)

uses fit_func as a fitting function and SVD LSQR routines to fit

input: x[ndata] y[ndata]: x and y values of the data
       wmax_thres:           threshold for SVD singular value cut off
       npoly:                degree of polynomial function:  a + a1 x + a2 x^2 + ...

       npara = npoly + 1

output: 
        chi2:  misfit chi square
        a[npara]: parameters
	siga[npara]: uncertainties 

we assume that we don't know the std of y

*/
void fit_poly(struct nr_datas *d,int ndata_orig, 
	      int npoly,COMP_PRECISION wmax_thres,
	      COMP_PRECISION *a, COMP_PRECISION *siga, 
	      COMP_PRECISION *chi2, 
	      my_boolean use_data_sigma,
	      my_boolean verbose,
	      COMP_PRECISION *sfac) /* sigma-like factor
				       determined from chi^2/n-m
				    */
{
  COMP_PRECISION *cov,*xloc,*yloc,*sigyloc;
  int i,j,npara,k,dof,ndata,nan;
  /* 
     
     take care of Nans
     
  */
  my_vecalloc(&xloc,ndata_orig,"fit_poly: xloc");
  my_vecalloc(&yloc,ndata_orig,"fit_poly: xloc");
  my_vecalloc(&sigyloc,ndata_orig,"fit_poly: xloc");
  for(i=ndata=0;i < ndata_orig;i++){
    if(finite(d[i].x) && finite(d[i].y)){
      xloc[ndata] = d[i].x;
      yloc[ndata] = d[i].y;
      sigyloc[ndata] = d[i].sigy;
      ndata++;
    }
  }
  nan = ndata_orig - ndata;
  if(verbose && nan)
    fprintf(stderr,"fit_poly: excluding %i nans out of %i\n",
	    nan,ndata_orig);
  if((ndata < 3)||(npoly<0)){
    fprintf(stderr,"fit_poly: parameter error: ndata: %i npoly: %i\n",
	    ndata,npoly);
    exit(-1);
  }
  /* 
     resize the parameter and uncertainty of parameter arrays 
  */
  npara = npoly + 1;
  if(!use_data_sigma){
    /* 
       we don't have actual uncertainties, need to fake them 
    */
    for(i=0;i<ndata;i++)
      sigyloc[i] = 1.0;		/* fake uncertainties */
  }
  /* covariance remains local */
  my_vecalloc(&cov,npara * npara,"fit_poly");	/* covariance */
  /* 
     perform the fit 
  */
  if(1)
    svdfit_driver(xloc,yloc,sigyloc,ndata,npara,wmax_thres,
		  (void (*)(void))(poly_fit_func),
		  npoly,a,cov,chi2);
  else
    lfit_driver(xloc,yloc,sigyloc,ndata,npara,
		(void (*)(void))(poly_fit_func),
		npoly,a,cov,chi2);

   if(!use_data_sigma){
    /* 
       rescale covariance matrix to obtain sigma in parameters
    */
    *sfac = *chi2/(COMP_PRECISION)((ndata-2));
    for(i=0;i< npara;i++)
      for(j=0;j< npara;j++)
	cov[i*npara+j] *= (*sfac);
  }
  /* 
     take sqrt of covariance matrix and obtain parameter uncertainties
  */
  for(i=0;i< npara;i++){
    for(j=0;j< npara;j++){
      k = i*npara+j;
      cov[k] += MY_SIGN(cov[k]) * sqrt(fabs(cov[k])); 
    }
    siga[i] = cov[i*npara+i]; 
  }
  dof = ndata - npara - 2;
  if(verbose){
    /* 
       output of results 
    */
    fprintf(stderr,"fit_poly: SVD fit %i polynomial to %i points: chi2: %11g rchi2: %8.3f (thres: %.2e)\n",
	    npoly,ndata,*chi2,sqrt(*chi2/dof),wmax_thres);
    for(i=0;i<=npoly;i++)
      fprintf(stderr,"fit_poly: a_%2i: %12g +/- %12g (%8.1f%%)\n",
	      i,a[i],siga[i],fabs(siga[i]/a[i])*100);
    for(i=0;i<npara;i++){
      fprintf(stderr,"fit_poly: |Cij| para %3i: ",i);
      for(j=0;j<=i;j++)
	fprintf(stderr,"%6.3f ",cov[i*npara+j]/siga[i]);
      fprintf(stderr,"\n");
    }
  }

  free(cov);
  free(xloc);
  free(yloc);
  free(sigyloc);
}



/* 

evaluate a full set of basis functions at location x

input: x:        location
       funcs:    basis function evlauation function
       npara:    number of parameters
       a[npara]: parameters

return value of model at location x

this has been programmed to be interoperable with the 
fitting function (0..n-1) style

*/
COMP_PRECISION evaluate_model(COMP_PRECISION x,int npara,
			      COMP_PRECISION *a,
			      void (*funcs)(),int func_para1)
{
  int i;
  COMP_PRECISION ymod,*base;
  my_vecalloc(&base,npara,"evaluate_model");
  /* get basis functions */
  (*funcs)(x,base,npara,func_para1);
  /* add up */
  for(ymod=0.0,i=0;i < npara;i++)
    ymod += base[i] * a[i];
  free(base);
  return ymod;
}

/* 
   subroutine to generate value of polynomial 
   fitting function at x
*/
void poly_fit_func(COMP_PRECISION x,COMP_PRECISION *y_base, 
		   int npara, int npoly)
{
  int i;
  y_base[0] = 1.0;		/* a_0 */
  for(i=1;i<= npoly;i++)
    y_base[i] = y_base[i-1] * x;
}

/* 
   evaluate polynomial of order n  (2 = power 2) at location 
   x
*/
COMP_PRECISION poly_val(COMP_PRECISION x, COMP_PRECISION *a, 
			int n)
{
  int i,np1;
  COMP_PRECISION tmp,res;
  np1 = n+1;
  res = a[0];
  tmp = x;
  for(i=1;i < np1;i++){
    res += a[i] * tmp;
    tmp *= x;
  }
  return res;
}
/* 
   
compute the reduced degree of freedom of ndata data values fitted by 
a npoly order polynomial

 */
int red_dof(int ndata, int npoly)
{
  return ndata - npoly -1;
}


/* 

numerical recipes Gauss-Jordan normal equation fitting procedure

modified such that funcs work in the regular C way (0...n-1)


*/

int nr_lfit(x,y,sig,ndat,a,ia,ma,covar,chisq,funcs,func_para1)
     COMP_PRECISION **covar,*chisq,a[],sig[],x[],y[];
     int ia[],ma,ndat,func_para1;
     void (*funcs)();
{
  
  int i,j,k,l,m,mfit=0,err;
  COMP_PRECISION ym,wt,sum,sig2i,**beta,*afunc;

  beta=nr_matrix(1,ma,1,1);
  afunc=nr_vector(1,ma);
  /* count the actual fitting parameters */
  for (j=1;j<=ma;j++)
    if (ia[j]) mfit++;
  if (mfit == 0) 
    nr_error("nr_lfit: no parameters to be fitted");
  for (j=1;j<=mfit;j++) {	/* init */
    for (k=1;k<=mfit;k++) 
      covar[j][k]=0.0;
    beta[j][1]=0.0;
  }
  /* assemble AT A matrix */
  for (i=1;i<=ndat;i++) {
    (*funcs)(x[i],(afunc+1),ma,func_para1);
    ym=y[i];
    if (mfit < ma) {
      for (j=1;j<=ma;j++)
	if (!ia[j]) 
	  ym -= a[j]*afunc[j];
    }
    sig2i=1.0/(sig[i]*sig[i]);
    for (j=0,l=1;l<=ma;l++) {
      if (ia[l]) {
	wt=afunc[l]*sig2i;
	for (j++,k=0,m=1;m<=l;m++)
	  if (ia[m]) 
	    covar[j][++k] += wt*afunc[m];
	beta[j][1] += ym*wt;
      }
    }
  }
  for (j=2;j<=mfit;j++)
    for (k=1;k<j;k++)
      covar[k][j]=covar[j][k];
  /* Gauss Jordan */
  if((err=nr_gaussj(covar,mfit,beta,1)) != 0){
    
    for (j=0,l=1;l<=ma;l++)
      if (ia[l]) 
	a[l]=beta[++j][1];
    /* evaluate misfit */
    *chisq=0.0;
    for (i=1;i<=ndat;i++) {
      (*funcs)(x[i],(afunc+1),ma,func_para1);
      for (sum=0.0,j=1;j<=ma;j++) 
	sum += a[j]*afunc[j];
      *chisq += NR_SQR((y[i]-sum)/sig[i]);
    }
    /* get covariance */
    nr_covsrt(covar,ma,ia,mfit);
  }
  free(afunc);
  nr_free_matrix(beta,1,ma,1,1);
  return err;
}



/* 

numerical recipes svdfit procedure

modified such that funcs work in the regular C way (0...n-1)

*/
void nr_svdfit(x,y,sig,ndata,a,ma,u,v,w,chisq,wmax_thres,funcs,func_para1)
     COMP_PRECISION **u,**v,*chisq,a[],sig[],w[],x[],y[],wmax_thres;
     int ma,ndata;
     int func_para1;
     void (*funcs)();
{
  int j,i;
  COMP_PRECISION wmax,tmp,thresh,sum,*b,*afunc;
  
  b=nr_vector(1,ndata);
  afunc=nr_vector(1,ma);
  for (i=1;i<=ndata;i++) {
    (*funcs)(x[i],(afunc+1),ma,func_para1);
    tmp=1.0/sig[i];
    for (j=1;j<=ma;j++) u[i][j]=afunc[j]*tmp;
    b[i]=y[i]*tmp;
  }
  nr_svdcmp(u,ndata,ma,w,v);
  wmax=0.0;
  for (j=1;j<=ma;j++)
    if (w[j] > wmax) wmax=w[j];
  thresh = wmax_thres * wmax;
  for (j=1;j<=ma;j++)
    if (w[j] < thresh) w[j]=0.0;
  nr_svbksb(u,w,v,ndata,ma,b,a);
  *chisq=0.0;
  for (i=1;i<=ndata;i++) {
    (*funcs)(x[i],(afunc+1),ma,func_para1);
    for (sum=0.0,j=1;j<=ma;j++) sum += a[j]*afunc[j];
    *chisq += (tmp=(y[i]-sum)/sig[i],tmp*tmp);
  }
  free(afunc);
  free(b);
}


