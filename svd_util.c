#include "precision.h"
#include "nr_defines.h"
/* 

SVD utility routines from numerical recipes



*/


/* 
   
SVD decomposition

 */
void nr_svdcmp(a,m,n,w,v)
     COMP_PRECISION **a,**v,w[];
     int m,n;
{
  int flag,i,its,j,jj,k,l,nm;
  COMP_PRECISION anorm,c,f,g,h,s,scale,x,y,z,*rv1;
  l=nm=0;
  rv1=nr_vector(1,n);
  g=scale=anorm=0.0;
  for (i=1;i<=n;i++) {
    l=i+1;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i <= m) {
      for (k=i;k<=m;k++) scale += fabs(a[k][i]);
      if (fabs(scale)>EPS_COMP_PREC) {
	for (k=i;k<=m;k++) {
	  a[k][i] /= scale;
	  s += a[k][i]*a[k][i];
	}
	f=a[i][i];
	g = -NR_SIGN(sqrt(s),f);
	h=f*g-s;
	a[i][i]=f-g;
	for (j=l;j<=n;j++) {
	  for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
	  f=s/h;
	  for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
	}
	for (k=i;k<=m;k++) a[k][i] *= scale;
      }
    }
    w[i]=scale *g;
    g=s=scale=0.0;
    if (i <= m && i != n) {
      for (k=l;k<=n;k++) scale += fabs(a[i][k]);
      if (fabs(scale)>EPS_COMP_PREC) {
	for (k=l;k<=n;k++) {
	  a[i][k] /= scale;
	  s += a[i][k]*a[i][k];
	}
	f=a[i][l];
	g = -NR_SIGN(sqrt(s),f);
	h=f*g-s;
	a[i][l]=f-g;
	for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
	for (j=l;j<=m;j++) {
	  for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
	  for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
	}
	for (k=l;k<=n;k++) a[i][k] *= scale;
      }
    }
    anorm=NR_FMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
  }
  for (i=n;i>=1;i--) {
    if (i < n) {
      if (fabs(g)>EPS_COMP_PREC) {
	for (j=l;j<=n;j++)
	  v[j][i]=(a[i][j]/a[i][l])/g;
	for (j=l;j<=n;j++) {
	  for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
	  for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
	}
      }
      for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
    }
    v[i][i]=1.0;
    g=rv1[i];
    l=i;
  }
  for (i=NR_IMIN(m,n);i>=1;i--) {
    l=i+1;
    g=w[i];
    for (j=l;j<=n;j++) a[i][j]=0.0;
    if (fabs(g)>EPS_COMP_PREC) {
      g=1.0/g;
      for (j=l;j<=n;j++) {
	for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
	f=(s/a[i][i])*g;
	for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
      }
      for (j=i;j<=m;j++) a[j][i] *= g;
    } else for (j=i;j<=m;j++) a[j][i]=0.0;
    ++a[i][i];
  }
  for (k=n;k>=1;k--) {
    for (its=1;its<=30;its++) {
      flag=1;
      for (l=k;l>=1;l--) {
	nm=l-1;
	if ((COMP_PRECISION)(fabs(rv1[l])+anorm) == anorm) {
	  flag=0;
	  break;
	}
	if ((COMP_PRECISION)(fabs(w[nm])+anorm) == anorm) break;
      }
      if (flag) {
	c=0.0;
	s=1.0;
	for (i=l;i<=k;i++) {
	  f=s*rv1[i];
	  rv1[i]=c*rv1[i];
	  if ((COMP_PRECISION)(fabs(f)+anorm) == anorm) break;
	  g=w[i];
	  h=nr_pythag(f,g);
	  w[i]=h;
	  h=1.0/h;
	  c=g*h;
	  s = -f*h;
	  for (j=1;j<=m;j++) {
	    y=a[j][nm];
	    z=a[j][i];
	    a[j][nm]=y*c+z*s;
	    a[j][i]=z*c-y*s;
	  }
	}
      }
      z=w[k];
      if (l == k) {
	if (z < 0.0) {
	  w[k] = -z;
	  for (j=1;j<=n;j++) v[j][k] = -v[j][k];
	}
	break;
      }
      if (its == 30) 
	nr_error("svdcmp: no convergence in 30 svdcmp iterations");
      x=w[l];
      nm=k-1;
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=nr_pythag(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+NR_SIGN(g,f)))-h))/x;
      c=s=1.0;
      for (j=l;j<=nm;j++) {
	i=j+1;
	g=rv1[i];
	y=w[i];
	h=s*g;
	g=c*g;
	z=nr_pythag(f,h);
	rv1[j]=z;
	c=f/z;
	s=h/z;
	f=x*c+g*s;
	g = g*c-x*s;
	h=y*s;
	y *= c;
	for (jj=1;jj<=n;jj++) {
	  x=v[jj][j];
	  z=v[jj][i];
	  v[jj][j]=x*c+z*s;
	  v[jj][i]=z*c-x*s;
	}
	z=nr_pythag(f,h);
	w[j]=z;
	if (z) {
	  z=1.0/z;
	  c=f*z;
	  s=h*z;
	}
	f=c*g+s*y;
	x=c*y-s*g;
	for (jj=1;jj<=m;jj++) {
	  y=a[jj][j];
	  z=a[jj][i];
	  a[jj][j]=y*c+z*s;
	  a[jj][i]=z*c-y*s;
	}
      }
      rv1[l]=0.0;
      rv1[k]=f;
      w[k]=x;
    }
  }
  free(rv1);
}


/* 
   covariance resorting for nr_lfit 
*/
void nr_covsrt(covar,ma,ia,mfit)
     COMP_PRECISION **covar;
     int ia[],ma,mfit;
{
  int i,j,k;
  for (i=mfit+1;i<=ma;i++)
    for (j=1;j<=i;j++) covar[i][j]=covar[j][i]=0.0;
  k=mfit;
  for (j=ma;j>=1;j--) {
    if (ia[j]) {
      for (i=1;i<=ma;i++) 
	NR_SWAP(covar[i][k],covar[i][j]);
      for (i=1;i<=ma;i++) 
	NR_SWAP(covar[k][i],covar[j][i]);
      k--;
    }
  }
}

COMP_PRECISION nr_pythag(a,b)
     COMP_PRECISION a,b;
{
  COMP_PRECISION absa,absb;
  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) return absa*sqrt(1.0+NR_SQR(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+NR_SQR(absa/absb)));
}

/* 
   compute covariance matrix 
*/
void nr_svdvar(v,ma,w,cvm)
     COMP_PRECISION **cvm,**v,*w;
     int ma;
{
  int k,j,i;
  COMP_PRECISION sum,*wti;
  wti=nr_vector(1,ma);
  for (i=1;i<=ma;i++) {
    wti[i]=0.0;
    if (w[i]) 
      wti[i]=1.0/(w[i]*w[i]);
  }
  for (i=1;i<=ma;i++) {
    for (j=1;j<=i;j++) {
      for (sum=0.0,k=1;k<=ma;k++) 
	sum += v[i][k]*v[j][k]*wti[k];
      cvm[j][i]=cvm[i][j]=sum;
    }
  }
  free(wti);
}
/* SVD back substitution */
void nr_svbksb(u,w,v,m,n,b,x)
     COMP_PRECISION **u,**v,b[],w[],x[];
     int m,n;
{
  int jj,j,i;
  COMP_PRECISION s,*tmp;

  tmp=nr_vector(1,n);
  for (j=1;j<=n;j++) {
    s=0.0;
    if (w[j]) {
      for (i=1;i<=m;i++) s += u[i][j]*b[i];
      s /= w[j];
    }
    tmp[j]=s;
  }
  for (j=1;j<=n;j++) {
    s=0.0;
    for (jj=1;jj<=n;jj++) s += v[j][jj]*tmp[jj];
    x[j]=s;
  }
  free(tmp);
}

