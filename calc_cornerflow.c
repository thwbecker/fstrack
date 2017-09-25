#include "fstrack.h"
/* 

*/

int main(int argc,char **argv)
{
  COMP_PRECISION alpha,vel,a,b,c,d,xcyl[3]={0,0,0},phi,
    xcart[3]={0,0,0},veccyl[3],veccart[3];
  int mode;
  if(argc != 4){
    fprintf(stderr,"%s mode(1-4) alpha[deg] vel\n",argv[0]);
    exit(-1);
  }
  
  sscanf(argv[1],"%i",&mode);
  sscanf(argv[2],"%lf",&alpha);
  sscanf(argv[3],"%lf",&vel);
  fprintf(stderr,"%s: corner flow, type %i: alpha: %g vel: %g\n",argv[0],mode,alpha,vel);
  alpha = DEG2RAD(alpha);
  /*  */
  calc_cornerflow_constants(mode-1, vel, alpha, &a, &b, &c, &d,argv);
  
  for(xcart[FSTRACK_X]=-10;xcart[FSTRACK_X]<=10+EPS_COMP_PREC;xcart[FSTRACK_X]+=0.05)
    for(xcart[FSTRACK_Y]=0;xcart[FSTRACK_Y]>=-10+EPS_COMP_PREC;xcart[FSTRACK_Y]-=0.05){
      cart2cyl(xcart,xcyl);	/* convert to cylindrical */
      cylvel(xcyl,veccyl,a,b,c,d);/* obtain cylindrical velocity components */
      cylvec2cartvec(xcyl, veccyl, veccart);// convert cly vel to cartesian vel
      
      phi = stream(xcyl,a,b,c,d);/* obtain cylindrical velocity components */

      
      fprintf(stdout,"%11g %11g %11g %11g %11g\n",xcart[FSTRACK_X],xcart[FSTRACK_Y],
	      veccart[FSTRACK_Y],veccart[FSTRACK_Y],phi);
    }
}

