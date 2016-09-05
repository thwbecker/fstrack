#include "fstrack.h"
#include "precision.h"
#include "misc.h"
//#include <malloc.h>
/*
  
  miscellaneous routines

  $Id: misc.c,v 1.42 2010/12/31 18:59:17 becker Exp $

*/


/*
  
open file safely

*/
FILE *myopen(char *filename, char *rwmode, char *program)
{
  FILE *stream;
  stream = fopen(filename,rwmode);
  if(!stream){
    fprintf(stderr,"%s: can not open file \"%s\" for mode %s, exiting\n",
	    program,filename,rwmode);
    exit(-1);
  }
  return stream;
}

/* calculate a sqrt function which returns zero and not underflow 
   is it's argument is small 
*/
COMP_PRECISION save_sqrt(COMP_PRECISION x)
{
  double xd;
  xd = (double)x;
  if(xd < 1.0e-16)
    return 0.0;
  else
    return ((COMP_PRECISION)sqrt(xd));
}
/*
  
allocate floating point vector

*/
void my_vecalloc(COMP_PRECISION **x,int n,char *message)
{
  *x = (COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*(size_t)n);
  if(! (*x))
    MEMERROR(message);
#ifdef MEM_ALLOC_FSTRACK_DEBUG
  fprintf(stderr,"my_vecalloc: newly allocating %.4f MB of memory (calling routine: %s)\n",
	  (float)(sizeof(COMP_PRECISION)*(size_t)n)/(float)ONE_MEGABYTE,message);
#endif
}
/* single prevision version */
void my_svecalloc(float **x,int n,char *message)
{
  *x = (float *)malloc(sizeof(float)*(size_t)n);
  if(! (*x))
    MEMERROR(message);
#ifdef MEM_ALLOC_FSTRACK_DEBUG
  fprintf(stderr,"my_svecalloc: newly allocating %.4f MB of memory (calling routine: %s)\n",
	  (float)(sizeof(float)*(size_t)n)/(float)ONE_MEGABYTE,message);
#endif
}

/* for velocity storage */
void my_vecvalloc(VPREC **x,int n,char *message)
{
  *x = (VPREC *)malloc(sizeof(VPREC)*(size_t)n);
  if(! (*x))
    MEMERROR(message);
#ifdef MEM_ALLOC_FSTRACK_DEBUG
  fprintf(stderr,"my_vecalloc: newly allocating %.4f MB of memory (calling routine: %s)\n",
	  (float)(sizeof(VPREC)*(size_t)n)/(float)ONE_MEGABYTE,message);
#endif
}



/*

  
allocate integer vector

*/
void my_ivecalloc(int **ix,int n,char *message)
{
  *ix = (int *)malloc(sizeof(int)*(size_t)n);
  if(! (*ix))
    MEMERROR(message);
#ifdef MEM_ALLOC_FSTRACK_DEBUG
  fprintf(stderr,"my_vecalloc: newly allocating %.4f MB of memory (calling routine: %s)\n",
	  (float)(sizeof(int)*(size_t)n)/(float)ONE_MEGABYTE,message);
#endif
}
/* reallocate floating point vector */
void my_vecrealloc(COMP_PRECISION **x,int n,char *message)
{
  *x = (COMP_PRECISION *)realloc(*x,
				 sizeof(COMP_PRECISION)*(size_t)n);
  if(!(*x))
    MEMERROR(message);
#ifdef MEM_ALLOC_FSTRACK_DEBUG
  fprintf(stderr,"my_vecrealloc: re-allocating %.4f MB of memory (calling routine: %s)\n",
	  (float)(sizeof(COMP_PRECISION)*(size_t)n)/(float)ONE_MEGABYTE,
	  message);
#endif
}


/* random number */
//
// returns uniformly distributed random number 
// 0 < x < 1
//
COMP_PRECISION myrand(long *seed)
{
  return((COMP_PRECISION)RAND_GEN(seed));
}

/*


  ran1 random number generator as of numerical recipes in C
  page 280


*/

#define RAN1_IA 16807
#define RAN1_IM 2147483647
#define RAN1_AM (1.0/RAN1_IM)
#define RAN1_IQ 127773
#define RAN1_IR 2836
#define RAN1_NTAB 32
#define RAN1_NDIV (1+(RAN1_IM-1)/RAN1_NTAB)
#define RAN1_EPS 5.0e-15
#define RAN1_RNMX (1.0-RAN1_EPS)

double ran1(long *idum)
{
  int j;
  long k;
  static long iy=0;
  static long iv[RAN1_NTAB];
  double temp;
  
  if (*idum <= 0 || !iy) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    for (j=RAN1_NTAB+7;j>=0;j--) {
      k=(*idum)/RAN1_IQ;
      *idum=RAN1_IA*(*idum-k*RAN1_IQ)-RAN1_IR*k;
      if (*idum < 0) *idum += RAN1_IM;
      if (j < RAN1_NTAB) iv[j] = *idum;
		}
    iy=iv[0];
  }
  k=(*idum)/RAN1_IQ;
  *idum=RAN1_IA*(*idum-k*RAN1_IQ)-RAN1_IR*k;
  if (*idum < 0) *idum += RAN1_IM;
  j=iy/RAN1_NDIV;
  iy=iv[j];
  iv[j] = *idum;
  if ((temp=RAN1_AM*iy) > RAN1_RNMX) return RAN1_RNMX;
  else return temp;
}
#undef RAN1_IA
#undef RAN1_IM
#undef RAN1_AM
#undef RAN1_IQ
#undef RAN1_IR
#undef RAN1_NTAB
#undef RAN1_NDIV
#undef RAN1_EPS
#undef RAN1_RNMX

/*


  ran2 number generator from numerical recipes, page 282


 */
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define RNMX (1.0-EPS_PREC)

double ran2(long *idum)
{
  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  double temp;
  static int ntabp7 = NTAB + 7;
  
  if (*idum <= 0) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=ntabp7;j>=0;j--) {
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;
  *idum=IA1*(*idum-k*IQ1)-k*IR1;
  if (*idum < 0) *idum += IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV;
  iy=iv[j]-idum2;
  iv[j] = *idum;
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef RNMX



//
// get Gaussian distribution with unity variance (or standard deviation)
//
double gasdev(long *idum)
{
  static int iset=0;
  static double gset;
  double fac,rsq,v1,v2;
  
  if  (iset == 0) {
    do {
      v1=2.0*RAND_GEN(idum)-1.0;
      v2=2.0*RAND_GEN(idum)-1.0;
      rsq=v1*v1+v2*v2;
    } while ((rsq >= 1.0) || (fabs(rsq)<EPS_COMP_PREC));
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  } else {
    iset=0;
    return gset;
  }
}
//
// this returns a random number between 0 < y < x
// uniformly distributed
//
COMP_PRECISION myrandnr(COMP_PRECISION x,long *seed)
{
  return((COMP_PRECISION)((double)x*RAND_GEN(seed)));
}
//
// Gauss distributed number with zero mean and x 
// standard deviation
//
COMP_PRECISION mygauss_randnr(COMP_PRECISION x,long *seed)
{
  return((COMP_PRECISION)((double)x*gasdev(seed)));
}

/* 

get linear interpolation factors f1 and f2=1-f1 at location xloc for x[n]
where x[] increases monotonousely, i is the left index, j the right

y = f1 * y[i] + f2 * y[j];

*/
void get_lin_int_weights(COMP_PRECISION xloc, COMP_PRECISION *x, int n,
			 int *i, int *j, COMP_PRECISION *f1, COMP_PRECISION *f2)
{
  *j = 0;
  while((*j < n) && (x[*j] < xloc))
    *j += 1;
  if(*j == 0)
    *j = 1;
  *i = *j - 1;
  *f2 = (xloc - x[*i])/(x[*j] - x[*i]);
  *f1 = 1.0 -  *f2;
}

COMP_PRECISION lin_inter(COMP_PRECISION xloc, COMP_PRECISION *x, COMP_PRECISION *y,
			 int n)
{
  int i, j;
  COMP_PRECISION f1,f2;
  get_lin_int_weights(xloc,x,n,&i,&j,&f1,&f2);
  return f1 * y[i] + f2 * y[j];
}
/* 

read symmetric sav tensor

*/
int read_sym_6by6(COMP_PRECISION *sav,FILE *in)
{
  int i,j,nr;
  nr = 0;
  for(i=0;i<6;i++)
    for(j=i;j<6;j++){
      /* read in Sav in upper right hand side format and fill in
	 symmetric pieces */
      nr += fscanf(in,FLT_FORMAT,(sav+i*6+j));
      sav[j*6+i] = sav[i*6+j];
    }
  if(nr != 21){
    fprintf(stderr,"%s: read error, not enough fields (%i out of 21) for 6x6 tensor in upper right hand side format\n",
	    "read_sym_6by6",nr);
    return 0;
  }else{
    return 1;
  }
}

COMP_PRECISION restrict_ani_factor(int restrict_ani, COMP_PRECISION depth)
{
  if(depth < 0){
    fprintf(stderr,"%s: error: expecting positive depth, %g\n",
	    "restrict_ani_factor",depth);
    exit(-1);
  }
  /* check if we need to restrict the anisotropy */
  switch(restrict_ani){
  case 0:
    return 1.0;
    break;
  case 1:			/* top zeroed out */
    if(depth <= 75)
      return 0.0;
    else
      return 1.0;
    break;
  case 2:			/* bottom zeroed out */
    if(depth >= 300)
      return 0.0;
    else
      return 1.0;
    break;
  default:
    fprintf(stderr,"%s: error: restrict_ani mode %i undefined\n",
	    "restrict_ani_factor",restrict_ani);
    exit(-1);
    break;
  }
  return 1.0;
}
