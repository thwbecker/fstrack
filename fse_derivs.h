#ifndef NO_UNDERSCORE
/* those for GCC with two underscore */

#ifdef GCC_USCR

#define fse_derivs fse_derivs__

#else

#define fse_derivs fse_derivs_

#endif



#endif

extern void fse_derivs(EXT_FTRN_PREC *,EXT_FTRN_PREC *,
		       EXT_FTRN_PREC *,int *,
		       VPREC *,VPREC *,VPREC *,
		       int *,int *,int *,EXT_FTRN_PREC *,
		       EXT_FTRN_PREC *,EXT_FTRN_PREC *,
		       EXT_FTRN_PREC *,int *,
		       int *,EXT_FTRN_PREC *,EXT_FTRN_PREC *,
		       EXT_FTRN_PREC *,int *,int *,
		       EXT_FTRN_PREC *,int *,int *, 
		       EXT_FTRN_PREC *);



