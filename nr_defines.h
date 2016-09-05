#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "precision.h"

#ifndef LOADED_MY_NR_DEFINES
/* 

abbreviated version of the numerical recipes headers

*/
#define NR_SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define MY_SIGN(x) (((x)>=0.0)?(1.0):(-1.0))

static int __tmp__iminarg1,__tmp__iminarg2;
#define NR_IMIN(a,b) (__tmp__iminarg1=(a),__tmp__iminarg2=(b),(__tmp__iminarg1) < (__tmp__iminarg2) ?\
        (__tmp__iminarg1) : (__tmp__iminarg2))


static COMP_PRECISION __tmp__sqrarg;
#define NR_SQR(a) ((fabs((__tmp__sqrarg=(a)))<EPS_COMP_PREC) ? (0.0) : (__tmp__sqrarg*__tmp__sqrarg))


static COMP_PRECISION __tmp__dmaxarg1,__tmp__dmaxarg2;
#define NR_FMAX(a,b) (__tmp__dmaxarg1=(a),__tmp__dmaxarg2=(b),(__tmp__dmaxarg1) > (__tmp__dmaxarg2) ? \
         (__tmp__dmaxarg1) : (__tmp__dmaxarg2))

#define NR_FREE_ARG char*
#define NR_END 1


static COMP_PRECISION __tmp_swap__;
#define NR_SWAP(a,b) {__tmp_swap__=(a);(a)=(b);(b)=__tmp_swap__;}
/* 

data structure

 */
struct nr_datas{
  COMP_PRECISION x,y,sigx,sigy;
};

/* nr_util.c */
int *nr_ivector(long, long);
void nr_free_matrix(COMP_PRECISION **, long, long, long, long);
COMP_PRECISION **nr_matrix(long, long, long, long);
COMP_PRECISION *nr_vector(long, long);
void nr_error(char *);
void nr_free_vector(COMP_PRECISION *, long, long);
int nr_gaussj(COMP_PRECISION **, int, COMP_PRECISION **, int);
/* datafit_util.c */
void svdfit_driver(COMP_PRECISION *, COMP_PRECISION *, COMP_PRECISION *, int, int, COMP_PRECISION, void (*)(void), int, COMP_PRECISION *, COMP_PRECISION *, COMP_PRECISION *);
int lfit_driver(COMP_PRECISION *, COMP_PRECISION *, COMP_PRECISION *, int, int, void (*)(void), int, COMP_PRECISION *, COMP_PRECISION *, COMP_PRECISION *);
void fit_poly(struct nr_datas *, int, int, COMP_PRECISION, COMP_PRECISION *, COMP_PRECISION *, COMP_PRECISION *, unsigned short, unsigned short, COMP_PRECISION *);
COMP_PRECISION evaluate_model(COMP_PRECISION, int, COMP_PRECISION *, void (*)(void), int);
void poly_fit_func(COMP_PRECISION, COMP_PRECISION *, int, int);
COMP_PRECISION poly_val(COMP_PRECISION, COMP_PRECISION *, int);
int red_dof(int, int);
int nr_lfit(COMP_PRECISION [], COMP_PRECISION [], COMP_PRECISION [], int, COMP_PRECISION [], int [], int, COMP_PRECISION **, COMP_PRECISION *, void (*)(void), int);
void nr_svdfit(COMP_PRECISION [], COMP_PRECISION [], COMP_PRECISION [], int, COMP_PRECISION [], int, COMP_PRECISION **, COMP_PRECISION **, COMP_PRECISION [], COMP_PRECISION *, COMP_PRECISION, void (*)(void), int);
/* svd_util.c */
void nr_svdcmp(COMP_PRECISION **, int, int, COMP_PRECISION [], COMP_PRECISION **);
void nr_covsrt(COMP_PRECISION **, int, int [], int);
COMP_PRECISION nr_pythag(COMP_PRECISION, COMP_PRECISION);
void nr_svdvar(COMP_PRECISION **, int, COMP_PRECISION *, COMP_PRECISION **);
void nr_svbksb(COMP_PRECISION **, COMP_PRECISION [], COMP_PRECISION **, int, int, COMP_PRECISION [], COMP_PRECISION []);


/* other stuff */


void my_vecalloc(COMP_PRECISION **,int ,char *);
void my_vecrealloc(COMP_PRECISION **,int ,char *);

void my_ivecalloc(int **,int ,char *);


#define LOADED_MY_NR_DEFINES
#endif

