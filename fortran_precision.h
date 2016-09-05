c
c fortran related definitions, only read by fortran code
c
c default is to use double precision
#ifndef SINGLE_PRECISION
#define COMP_PRECISION double precision
      implicit real*8 (a-h,o-z)
      parameter (pi = 3.14159265358979d0, twopi = 6.28318530717959d0)
      parameter (eps_prec = 5.0d-15)
#else
#define COMP_PRECISION real*4
      implicit real*4 (a-h,o-z)
      parameter (pi = 3.1415927, twopi = 6.2831853)
      parameter (eps_prec = 5.0e-7)
#endif

#define VPRECF real*4
c
c define more meaningful abbreviations for work array
c which holds the x vector (3-D) and the deformation matrix (3 x 3)
c x(j)
#define XTRACER(j) work(j)  
c def(1,j)
#define DEF1(j) work(3+j) 
c def(2,j)
#define DEF2(j) work(6+j) 
c def(3,j)
#define DEF3(j) work(9+j) 


c some array helpers
#define _R_ 1
#define _T_ 2
#define _P_ 3
