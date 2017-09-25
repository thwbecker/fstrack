#include "fstrack_flow.h"


int main(int argc, char **argv)
{
  int n,i,j,nharm,npara;
  long int seed=-1;
  COMP_PRECISION *x,*y,*a=NULL,*siga=NULL,xmod,z,s[36],
    av[2]={1,1},bv[2]={.9,.9},t,p,
    ave,adev,sdev,var,skew,curt,dx,arg,chi2,*sigy;

  myrand(&seed);

  //exit(-1);
  n=1000;  

  my_vecalloc(&x,n,"");
  my_vecalloc(&y,n,"");
  my_vecalloc(&sigy,n,"");

  dx = TWO_PI/(n-1);

  for(i=0;i<n;i++){
    x[i] = dx * i;
    y[i] = (float)i/(float)n + 
      sin(x[i] +0.1) + cos(2*x[i] +0.1) + myrandnr(0.1,&seed);
    //x[i] = mygauss_randnr(1,&seed);
    //fprintf(stdout,"%g\n",x[i]);
  }
  
  stat_moments(y,n,&ave,&adev,&sdev,&var,&skew,&curt,TRUE);

  //for(nharm=0;nharm<6;nharm++){
  for(nharm=0;nharm<3;nharm++){
    npara = 2 + 2 *nharm;

    my_vecrealloc(&a,npara,"");
    my_vecrealloc(&siga,npara,"");
    fit_harmonic(x,y,sigy,n,nharm,1e-5,a,siga,
		 &chi2,FALSE,TRUE,TRUE);
  }
  nharm--;
  for(i=0;i<n;i++){
    
    printf("%g %g %g\n",x[i],y[i],
	   evaluate_model(x[i],npara,a, 
			  (void (*)(void))(harm_fit_func),
			  nharm));
  }
  
  free(x);free(y);free(a);free(siga);

  return 0;
}
