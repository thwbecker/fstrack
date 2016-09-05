#include <stdio.h>
#include <stdlib.h>
#ifndef MEMERROR
#define MEMERROR(x) {fprintf(stderr,"memory allocation error: %s\n",x);exit(-1);}
#endif
COMP_PRECISION save_sqrt(COMP_PRECISION);
void my_vecalloc(COMP_PRECISION **, int, char *);
void my_ivecalloc(int **, int, char *);
void my_vecrealloc(COMP_PRECISION **, int, char *);
COMP_PRECISION myrand(long *);
COMP_PRECISION ran1(long *);
COMP_PRECISION ran2(long *);
COMP_PRECISION gasdev(long *);
COMP_PRECISION myrandnr(COMP_PRECISION, long *);
COMP_PRECISION mygauss_randnr(COMP_PRECISION, long *);
void get_lin_int_weights(COMP_PRECISION, COMP_PRECISION *, int, int *, int *, COMP_PRECISION *, COMP_PRECISION *);
COMP_PRECISION lin_inter(COMP_PRECISION, COMP_PRECISION *, COMP_PRECISION *, int);
