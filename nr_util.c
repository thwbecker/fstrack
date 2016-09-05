#include "precision.h"
#include "nr_defines.h"
/* other stuff */
#define NR_END 1
#define NR_FREE_ARG char*


/* 

general numerical recipes routines


*/
int *nr_ivector(nl,nh)
     long nh,nl;
     /* allocate an int vector with subscript range v[nl..nh] */
{
  int *v;
  
  v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
  if (!v) nr_error("nr_ivector: allocation failure in ivector()");
  return v-nl+NR_END;
}

void nr_free_matrix(m,nrl,nrh,ncl,nch)
COMP_PRECISION **m;
long nch,ncl,nrh,nrl;
/* free a COMP_PRECISION matrix allocated by matrix() */
{
  free((NR_FREE_ARG) (m[nrl]+ncl-NR_END));
  free((NR_FREE_ARG) (m+nrl-NR_END));
}

COMP_PRECISION **nr_matrix(nrl,nrh,ncl,nch)
long nch,ncl,nrh,nrl;
/* allocate a COMP_PRECISION matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  COMP_PRECISION **m;
  
  /* allocate pointers to rows */
  m=(COMP_PRECISION **) malloc((size_t)((nrow+NR_END)*sizeof(COMP_PRECISION*)));
  if (!m) nr_error("nr_matrix: allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;
  
  /* allocate rows and set pointers to them */
  m[nrl]=(COMP_PRECISION *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(COMP_PRECISION)));
  if (!m[nrl]) nr_error("nr_matrix: allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;
  
  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
  
  /* return pointer to array of pointers to rows */
  return m;
}
COMP_PRECISION *nr_vector(nl,nh)
long nh,nl;
/* allocate a vector with subscript range v[nl..nh] */
{
  COMP_PRECISION *v;
  
  v=(COMP_PRECISION *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(COMP_PRECISION)));
  if (!v) nr_error("nr_vector: allocation failure in vector()");
  return v-nl+NR_END;
}

void nr_error(char *error_text)
/* Numerical Recipes standard error handler */
{
  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(-1);
}
void nr_free_vector(v,nl,nh)
COMP_PRECISION *v;
long nh,nl;
/* free a float vector allocated with vector() */
{
  free((NR_FREE_ARG) (v+nl-NR_END));
}
/* 
   Gauss Jordan matrix elimation 
   returns non-zero, if error
*/
int nr_gaussj(a,n,b,m)
     COMP_PRECISION **a,**b;
     int m,n;
{
  int *indxc,*indxr,*ipiv;
  int i,icol,irow,j,k,l,ll;
  COMP_PRECISION big,dum,pivinv;
  //for(i=1;i<=n;i++)for(j=1;j<=n;j++)fprintf(stderr,"%3i %3i %11g\n",i,j,a[i][j]);
  icol=irow=0;
  indxc=nr_ivector(1,n);
  indxr=nr_ivector(1,n);
  ipiv=nr_ivector(1,n);
  for (j=1;j<=n;j++) ipiv[j]=0;
  for (i=1;i<=n;i++) {
    big=0.0;
    for (j=1;j<=n;j++)
      if (ipiv[j] != 1)
	for (k=1;k<=n;k++) {
	  if (ipiv[k] == 0) {
	    if (fabs(a[j][k]) >= big) {
	      big=fabs(a[j][k]);
	      irow=j;
	      icol=k;
	    }
	  } else if (ipiv[k] > 1) {
	    fprintf(stderr,"nr_gaussj: Singular Matrix-1");
	    return -1;
	  }
	}
    ++(ipiv[icol]);
    if (irow != icol) {
      for (l=1;l<=n;l++) NR_SWAP(a[irow][l],a[icol][l]);
      for (l=1;l<=m;l++) NR_SWAP(b[irow][l],b[icol][l]);
    }
    indxr[i]=irow;
    indxc[i]=icol;
    if (a[icol][icol] == 0.0){
      fprintf(stderr,"gaussj: Singular Matrix-2");
      return -2;
    }
    pivinv=1.0/a[icol][icol];
    a[icol][icol]=1.0;
    for (l=1;l<=n;l++) a[icol][l] *= pivinv;
    for (l=1;l<=m;l++) b[icol][l] *= pivinv;
    for (ll=1;ll<=n;ll++)
      if (ll != icol) {
	dum=a[ll][icol];
	a[ll][icol]=0.0;
	for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
	for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
      }
  }
  for (l=n;l>=1;l--) {
    if (indxr[l] != indxc[l])
      for (k=1;k<=n;k++)
	NR_SWAP(a[k][indxr[l]],a[k][indxc[l]]);
  }
  free(ipiv);
  free(indxr);
  free(indxc);
  return 0;
}
