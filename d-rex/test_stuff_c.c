#include "drex.h"

int main(void)
{
  COMP_PRECISION sav[36],alpha,beta=0,gamma=0,
    savr[36],savr21[21];

  drex_get_sav_constants(sav,DREX_C_OL_KR, 0,0);
  for(alpha=0;alpha<360;alpha+=30){
    for(beta=0;beta<180;beta+=20){
      for(gamma=0;gamma<180;gamma+=30){
	drex_rotate_6x6_ftrn(sav,savr,&alpha,&beta,&gamma);
	drex_fullsym6_ftrn(savr);
	drex_v21d_ftrn(savr,savr21);
 
	fprintf(stderr,"%g %g %g %g\n",
		alpha,beta,gamma,
		sqrt(vecdotp(savr21,savr21,21)));
      }
    }
  }

  return 0;
}
/* dot product */
COMP_PRECISION vecdotp(COMP_PRECISION *a,COMP_PRECISION *b,
		       int n)
{
  int i;
  COMP_PRECISION x;
  for(x=0.0,i=0;i < n;i++)
    x += a[i] * b[i];
  return x;
}
void transpose_3by3_inplace(COMP_PRECISION *a)
{
  COMP_PRECISION c[9];
  c[0] = a[0];c[1] = a[3];c[2] = a[6];
  c[3] = a[1];c[4] = a[4];c[5] = a[7];
  c[6] = a[2];c[7] = a[5];c[8] = a[8];
  a_equals_b_vector(a,c,9);	/* to allow calling as (x,x) */
}
void transpose_nbyn_inplace(COMP_PRECISION *a,int n)
{
  COMP_PRECISION *c;
  int nn,j,i;
  nn = n * n;
  my_vecalloc(&c,nn,"transpose_nbyn_inplace");
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      c[j*n+i] = a[i*n+j];
  a_equals_b_vector(a,c,nn);	
  free(c);
}
