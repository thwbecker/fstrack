#include "fstrack_flow.h"
/*

  reads in velocities and extracts the strain-rate matrix G

  prints the latter to stdout


 */

int main(int argc, char **argv)
{
  COMP_PRECISION x[12],vp[3],e[9],lon,lat,z;
  int n,preject;
  struct mod *model;
  model=(struct mod *)calloc(1,sizeof(struct mod));
  set_defaults(model);

  read_vel_grids(model,model->verbose,TRUE  /* zero out radial flow at boundary? */);
  

  fprintf(stderr,"%s: reading lon lat z locations from stdin, using velocities at time: %g\n",
	  argv[0],model->itime);
  fprintf(stderr,"%s: output format: lon lat z e_rr e_rt e_rp e_tt e_tp e_pp sqrt(0.5 E:E) v_r v_theta v_phi\n",
	   argv[0]);
  n=preject=0;
  //
  // read in locations at which to determine strain-rates
  while(fscanf(stdin,THREE_FLT_FORMAT,&lon,&lat,&z)==3){
    if((z < 0 )||(z>R_E)){
      fprintf(stderr,"%s: error: z should be >= 0 (and in km)\n",
	      argv[0]);
      exit(-1);
    }
    x[FSTRACK_PHI] = LONGITUDE2PHI(lon);
    x[FSTRACK_THETA] = LATITUDE2THETA(lat);
    x[FSTRACK_R] = ND_RADIUS(z);
    /* compute velocities and strain-rates in spherical system */
    calc_vel_and_strain_rate_units(model->itime,x,model->dp,vp,e,
				   FALSE);
    //
    // output in: 
    // 
    // lon lat z e_rr e_rt e_rp e_tt e_tp e_pp 
    //   ... sqrt(0.5 E:E) v_r v_theta v_phi
    //
    // format. velocities in cm/yr
    //
    fprintf(stdout,"%11.3e %11.3e %11.3e   %11.3e %11.3e %11.3e %11.3e %11.3e %11.3e   %11.3e  %11.3e %11.3e %11.3e \n",
	    lon,lat,z,
	    e[RR],e[RT],e[RP],e[TT],e[TP],e[PP],
	    sqrt(0.5 * double_dot(e,e)),
	    vp[FSTRACK_R],vp[FSTRACK_THETA],vp[FSTRACK_PHI]);
    n++;
  }
  fprintf(stderr,"%s: processed %i locations, rejected %i polar locations\n",
	  argv[0],n,preject);
  return 0;
}
