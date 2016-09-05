#include "fstrack.h"
/* 

pressure model

input is depth in [km]

pressure returned is in GPa, linearly interpolated from PREM

*/

#define PMODEL_P_N 25

/*  */
COMP_PRECISION pressure_model(COMP_PRECISION z)
{
  COMP_PRECISION p;
  COMP_PRECISION zpp[PMODEL_P_N]={0,2.999,3.001,14.999,15.0001,24.399,24.401,40,60,79.999,80.001,115,150,185,219.999,220.001,265,310,355,399.999,400.001,450,500,550,600};
  COMP_PRECISION pp[PMODEL_P_N]={0.0001,0.0299,0.0303,0.3364,0.3370,0.6040,0.6043,1.1239,1.7891,2.4539,2.4546,3.6183,4.7824,5.9466,7.1108,7.1115,8.6497,10.2027,11.7702,13.3520,13.3527,15.2251,17.1311,19.0703,21.0425};
  if((z < 0)||(z > zpp[PMODEL_P_N-1])){
      fprintf(stderr,"temperature_model: error: z (%g) out of bounds (0 .. %g) for pressure_model\n",
	      z,zpp[PMODEL_P_N-1]);
      exit(-1);
  }
  return lin_inter(z,zpp,pp,PMODEL_P_N);
}
