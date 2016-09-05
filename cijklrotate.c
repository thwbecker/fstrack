#include "fstrack.h"

/* 

read vera style (3,3,3,3) elasticity tensor,
rotate, and write Cijkl format


*/

int main(int argc, char **argv)
{
  my_boolean rotate = FALSE;
  COMP_PRECISION alpha=0,beta=0,gamma=0,c[81],rot[9],cr[81];
  /* fortran I/O */
  int fp = 5;			/* stdin */
  int op = 6;			/* stdout  */
  if(argc < 2){
    fprintf(stderr,"%s alpha [beta, %g] [gamma, %g]\n",
	    argv[0],beta,gamma);
    fprintf(stderr,"\t reads Cijkl tensor from stdin, rotates, and writes to stdout\n");
    exit(-1);
  }
  sscanf(argv[1],"%lf",&alpha);
  if(argc > 2)
    sscanf(argv[2],"%lf",&beta);
  if(argc > 3)
    sscanf(argv[3],"%lf",&gamma);
  /* convert angles */
  alpha = DEG2RAD(alpha);
  beta  = DEG2RAD(beta);
  gamma = DEG2RAD(gamma);
  if(fabs(beta) > EPS_COMP_PREC){
    fprintf(stderr,"%s: WARNING: rotating tensor by %g degrees down from horizontal (beta angle)\n",
	    argv[0],beta);
    rotate = TRUE;
  }
  if(fabs(alpha) > EPS_COMP_PREC){
    fprintf(stderr,"%s: WARNING: rotating tensor by %g degrees around vertical (alpha angle)\n",
	    argv[0],alpha);
    rotate = TRUE;
  }
  if(fabs(gamma) > EPS_COMP_PREC){
    fprintf(stderr,"%s: WARNING: rotating tensor by %g degrees around x (gamma angle)\n",
	    argv[0],gamma);
    rotate = TRUE;
  }
  /* read in tensor from stdin */
  fprintf(stderr,"%s: reading Cijkl from stdin, rotating by %g, %g, %g\n",
	  argv[0],alpha,beta,gamma);
  vera_read_cijkl(c,&fp);	/* read C in fortran format */
  if(rotate){
    /* get rotation matrix */
    drex_calc_rotmat_cart(rot,&alpha,&beta,&gamma);
    /* rotate a 3,3,3,3 tensor Fortran style */
    drex_rot4(c,rot,cr);
    /* copy over */
    a_equals_b_vector(c,cr,81);
  }
  /* print cijkl tensor fortran style */
  vera_print_cijkl_ftrn(c,&op);

}
