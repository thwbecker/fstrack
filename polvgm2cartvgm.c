/*     
       routine reads r theta phi TIME d_r(vr) d_t(vr) d_p(vr) ...
                                      d_r(vt) d_t(vt) d_p(vt) ...
		                      d_r(vp) d_t(vp) d_p(vp) 

       and converts the velocity gradient matrix to cartesian
       
       output is

       x y z TIME d_x vx d_y vx d_z vx ... 
                  d_x vy d_y vy d_z vy ...
                  d_x vz d_y vz d_z vz
       
  
      $Id: polvgm2cartvgm.c,v 1.4 2004/04/17 19:42:52 becker Exp $
*/
  
#include "fstrack.h"

int main(int argc, char **argv)
{
  COMP_PRECISION pvgm[3][3],cvgm[3][3],time,px[3],cx[3];
  int n=0;
  /* read in VGM in spherical coordinates */
  while(fscanf(stdin,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
	       (px+FSTRACK_R),(px+FSTRACK_THETA),(px+FSTRACK_PHI),&time,
	       &pvgm[FSTRACK_R][FSTRACK_R],    &pvgm[FSTRACK_R][FSTRACK_THETA],    &pvgm[FSTRACK_R][FSTRACK_PHI],
	       &pvgm[FSTRACK_THETA][FSTRACK_R],&pvgm[FSTRACK_THETA][FSTRACK_THETA],&pvgm[FSTRACK_THETA][FSTRACK_PHI],
	       &pvgm[FSTRACK_PHI][FSTRACK_R],  &pvgm[FSTRACK_PHI][FSTRACK_THETA],  &pvgm[FSTRACK_PHI][FSTRACK_PHI])==13){
    rtp2xyz(px,cx);		/* convert location from polar to cartesian*/
    /* convert matrix from polar to cartesian system */
    polar_to_cart_mat_at_r3x3(pvgm,cvgm,px);
    fprintf(stdout,"%16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n",
	    cx[FSTRACK_X],cx[FSTRACK_Y],cx[FSTRACK_Z],time,
	    cvgm[FSTRACK_X][FSTRACK_X],cvgm[FSTRACK_X][FSTRACK_Y],cvgm[FSTRACK_X][FSTRACK_Z],
	    cvgm[FSTRACK_Y][FSTRACK_X],cvgm[FSTRACK_Y][FSTRACK_Y],cvgm[FSTRACK_Y][FSTRACK_Z],
	    cvgm[FSTRACK_Z][FSTRACK_X],cvgm[FSTRACK_Z][FSTRACK_Y],cvgm[FSTRACK_Z][FSTRACK_Z]);
    if(fabs(pvgm[FSTRACK_X][FSTRACK_X]+pvgm[FSTRACK_Y][FSTRACK_Y]+pvgm[FSTRACK_Z][FSTRACK_Z]) > 1e-7){ /* check trace, should be zero */
      fprintf(stderr,"trace: %11g %11g\n",
	      pvgm[FSTRACK_X][FSTRACK_X]+pvgm[FSTRACK_Y][FSTRACK_Y]+pvgm[FSTRACK_Z][FSTRACK_Z],
	      cvgm[FSTRACK_X][FSTRACK_X]+cvgm[FSTRACK_Y][FSTRACK_Y]+cvgm[FSTRACK_Z][FSTRACK_Z]);
    }
    n++;
  }
  if(n)
    fprintf(stderr,"%s: converted %i entries\n",argv[0],n);
  else{
    fprintf(stderr,"%s: error, n is zero\n",argv[0]);
    exit(-1);
  }
  return 0;
}
