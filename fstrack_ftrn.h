/*

 fortran routines
*/
#ifdef GFORTRAN_STDOUT_UNIT
#define FTRN_STDOUT GFORTRAN_STDOUT_UNIT
#else
#define FTRN_STDOUT 6
#endif
#ifdef GFORTRAN_STDIN_UNIT
#define FTRN_STDIN GFORTRAN_STDIN_UNIT
#else
#define FTRN_STDIN 5
#endif
