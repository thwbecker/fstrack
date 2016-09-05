#include "fstrack.h"
/*

  read in tracer.f.s formay, ie.

  lon lat r lrr lrt lrp ltt ltp lpp time dz dx

  and prints cartesian eigenvectors to stdout in
		       
  time x y z e1 e2 e3 e1x e1y e1z e2x e2y e2z e3x e3y e3z

  format
		      
 */
int main(void)
{
  COMP_PRECISION evec[9],ecvec[9],eval[3],l[9],xc[3],xp[3],t,det;
  while(fscanf(stdin,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %*f %*f",
	       (xp+2),(xp+1),(xp),(l+RR),(l+RP),(l+RT),(l+TT),(l+TP),(l+PP),&t)==10){
    xp[2] = LONGITUDE2PHI(xp[2]);//phi 
    xp[1] = LATITUDE2THETA(xp[1]);//theta
    xp[0] = ND_RADIUS(xp[0]);//r
    calc_eigensystem_sym(l,eval,evec,TRUE);
    det=((eval[2]*eval[1]*eval[0])-1.0)/3.;// remove possible spurious det != 1
    l[RR] -= det;l[TT] -= det;l[PP] -= det;
    calc_eigensystem_sym(l,eval,evec,TRUE);
    rtp2xyz(xp,xc);// location in cartesian 
    polar_vec2cart_vec_at_xp((evec+6),(ecvec+6),xp);// e1 in cart
    polar_vec2cart_vec_at_xp((evec+3),(ecvec+3),xp);// e2 in cart
    polar_vec2cart_vec_at_xp((evec  ),(ecvec  ),xp);// e3 in cart
    fprintf(stdout,"%11g %11g %11g %11g %13.6e %13.6e %13.6e %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",
	    t,xc[FSTRACK_X],xc[FSTRACK_Y],xc[FSTRACK_Z],
	    sqrt(eval[2]),sqrt(eval[1]),sqrt(eval[0]),
	    ecvec[6+FSTRACK_X],ecvec[6+FSTRACK_Y],ecvec[6+FSTRACK_Z],
	    ecvec[3+FSTRACK_X],ecvec[3+FSTRACK_Y],ecvec[3+FSTRACK_Z],
	    ecvec[  FSTRACK_X],ecvec[  FSTRACK_Y],ecvec[  FSTRACK_Z]);
  }
  return 0;
}
