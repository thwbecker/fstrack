#include "drex.h"

/*
  
allocate floating point vector

*/
void my_vecalloc(COMP_PRECISION **x,int n,char *message)
{
  *x = (COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*(size_t)n);
  if(! (*x))
    MEMERROR(message);
}
/* single prevision version */
void my_svecalloc(float **x,int n,char *message)
{
  *x = (float *)malloc(sizeof(float)*(size_t)n);
  if(! (*x))
    MEMERROR(message);
}

/*
  
re-allocate floating point vector

*/
void my_vecrealloc(COMP_PRECISION **x,int n,char *message)
{
  *x = (COMP_PRECISION *)realloc(*x,sizeof(COMP_PRECISION)*(size_t)n);
  if(! (*x))
    MEMERROR(message);
}
//  \vec{a} = \vec{b}
void a_equals_b_vector(COMP_PRECISION *a,COMP_PRECISION *b,
		       int n)
{
  int i;
  for(i=0;i<n;i++)
    a[i] = b[i];
}
/* 
   b = b - a
 */
void sub_a_from_b_vector(COMP_PRECISION *a, COMP_PRECISION *b,
			 int n)
{
  int i;
  for(i=0;i<n;i++)
    b[i] -= a[i];
}
/* 
   compute the transpose of B and assign it to A

   A(3,3) = B(3,3)^T

   a and b are passed as 9 element vectors
*/
void a_is_b_transpose_3by3(COMP_PRECISION *a,COMP_PRECISION *b)
{
  COMP_PRECISION c[9];		/* don't delete this! */
  c[0] = b[0];c[1] = b[3];c[2] = b[6];
  c[3] = b[1];c[4] = b[4];c[5] = b[7];
  c[6] = b[2];c[7] = b[5];c[8] = b[8];
  a_equals_b_vector(a,c,9);	/* to allow calling as (x,x) */
}
void transpose_3by3_inplace(COMP_PRECISION *a)
{
  COMP_PRECISION c[9];
  c[0] = a[0];c[1] = a[3];c[2] = a[6];
  c[3] = a[1];c[4] = a[4];c[5] = a[7];
  c[6] = a[2];c[7] = a[5];c[8] = a[8];
  a_equals_b_vector(a,c,9);	/* to allow calling as (x,x) */
}
void a_is_b_transpose_nbyn(COMP_PRECISION *a,COMP_PRECISION *b,
			   int n)
{
  COMP_PRECISION *c;
  int nn,j,i;
  nn = n * n;
  my_vecalloc(&c,nn,"a_is_b_transposenbyn");
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      c[j*n+i] = b[i*n+j];
  a_equals_b_vector(a,c,nn);	/* to allow calling as (x,x) */
  free(c);
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
/* norm */
COMP_PRECISION norm(COMP_PRECISION *x, int n)
{
  return sqrt(vecdotp(x,x,n));
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
void print_vector(COMP_PRECISION *t, int n, FILE *out)
{
  int i;
  for(i=0;i<n;i++)
    fprintf(out,"%11g ",t[i]);
  fprintf(out,"\n");
}
