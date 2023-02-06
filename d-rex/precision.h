#include <stdio.h>
#include <math.h>
#include <stdlib.h>
//
// precision related, if changed, rerun 'make proto' for the auto_proto.h file
//
#ifndef MY_PRECISION_HAS_BEEN_SOURCED

#ifndef SINGLE_PRECISION // double precision is the default

#define COMP_PRECISION double
#define EXT_FTRN_PREC double
#define FLT_FORMAT "%lf"
#define TWO_FLT_FORMAT "%lf %lf"
#define THREE_FLT_FORMAT "%lf %lf %lf"
#define FOUR_FLT_FORMAT "%lf %lf %lf %lf"
#define EPS_PREC 8.0e-15
// eispack symmetric real double precision eigensystem routine
#define SYMREAL_ES_ROUTINE rs_
#define REAL_ES_ROUTINE rg_

#else // single precision

#define COMP_PRECISION float
#define EXT_FTRN_PREC float
#define FLT_FORMAT "%f"
#define TWO_FLT_FORMAT "%f %f"
#define THREE_FLT_FORMAT "%f %f %f"
#define FOUR_FLT_FORMAT "%f %f %f %f"
#define EPS_PREC 7.0e-8
// eispack symmetric real single precision eigensystem routine
#define SYMREAL_ES_ROUTINE s_rs_
#define REAL_ES_ROUTINE s_rg_

#endif

// needed for compatibility
#ifndef EPS_COMP_PREC
#define EPS_COMP_PREC EPS_PREC
#endif

#define my_boolean unsigned short

#ifndef TRUE 
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

#define MY_PRECISION_HAS_BEEN_SOURCED
#endif
