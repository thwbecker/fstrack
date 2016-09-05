#include "fstrack.h"
/* 

temperature model

input is depth in [km]


temperature returned is in Kelvin

T(0) is defined to be 298K, room temperature


*/

#define TMODEL_S_N 14
COMP_PRECISION temperature_model(COMP_PRECISION z,int mode)
{
  COMP_PRECISION t;
  COMP_PRECISION zts[TMODEL_S_N]={0,  11,  70, 119.5, 120.5, 170, 219.5, 220.5, 270, 320, 370, 419.5, 420.5, 470},
    ts[TMODEL_S_N]={            298, 550,1460,  1670,  1670,1760,1838,  1838,  1910,1975,2036,2085,  2085,  2130};
  switch(mode){
  case 0:			/* constant */
    t = 1400;
    break;
  case 1:			/* Stacey 75 oceanic up to 410 */
    if((z<0)||(z>zts[TMODEL_S_N-1])){
      fprintf(stderr,"temperature_model: error: z (%g) out of bounds (0 .. %g) for Stacey (1975)\n",
	      z,zts[TMODEL_S_N-1]);
      exit(-1);
    }
    t = lin_inter(z,zts,ts,TMODEL_S_N);
    break;
  default:
    fprintf(stderr,"temperature_model: error: mode %i undefined\n",mode);
    exit(-1);
    break;
  }
  return t;
}
