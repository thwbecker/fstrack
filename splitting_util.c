/* 
   
routines to perform analysis of elastic tensors in terms of splitting
(see also single_layer)


$Id: splitting_util.c,v 1.8 2016/09/05 04:44:47 becker Exp $

*/
#include "fstrack.h"
  
/* 

do some statistics on the azimuthal dependence of splitting

input:

n_azi:       number of samples with azimuth
nharm:       number of harmonics to fit
xazi[]:      azimuth samples on x axis in degrees
fazi[]:      fast azimuth of splitting in degrees
dt[]:        splitting times as functions of xazi azimuth [s]

output:

npara:     number of parameters = 
mean_fazi: mean fast azimuth
std_fazi:  standard deviation
mean_dt:   mean splitting time
std_dt:    standard deviation

output of fitting parameters: ALL ARE NORMALIZED TO THE STD

npara = 2 + nharm * 2

a_fazi[npara]:     pass as NULL: a_0 + a_1 *x + a_2 sin(x) + a_3 cos(x) ... terms
a_dt[npara]:       same for dt
ma_fazi[nharm],sma_fazi[nharm]: magnitude of the harmonic terms and their 
                                uncertainties, normalized by 1/sqrt(2)

ma_dt,sma_dt: same for dt



note that to assemble the best-fit functions, you will need to 
do

fazi_model = 	 mean_fazi + std_fazi * 
		 evaluate_model(razi,npara,a_fazi, 
				(void (*)(void))(harm_fit_func),
				nharm),

and


dt_model=        mean_dt + std_dt * evaluate_model(razi,npara,a_dt, 
				(void (*)(void))(harm_fit_func),
				nharm));			

				(note that dt is not reduced by the mean)

*/
void analyze_splitting_azi_dep(int n_azi,int nharm,
			       COMP_PRECISION *xazi,
			       COMP_PRECISION *fazi,
			       COMP_PRECISION *dt,
			       COMP_PRECISION *mean_fazi,
			       COMP_PRECISION *std_fazi,
			       COMP_PRECISION *mean_dt,
			       COMP_PRECISION *std_dt,
			       COMP_PRECISION *a_fazi,
			       COMP_PRECISION *a_dt,
			       COMP_PRECISION *ma_fazi,
			       COMP_PRECISION *sma_fazi,
			       COMP_PRECISION *ma_dt,
			       COMP_PRECISION *sma_dt,
			       my_boolean verbose,
			       int use_dt_for_avg)
{
  COMP_PRECISION chi2,*dfazi,dt_stats[6],dfazi_stats[6],
    *dfazi_sig,*dt_copy,*rxazi,*harm_siga,mean_vec_dt;
  static COMP_PRECISION fit_svd_cutoff = 1e-4; /* for harmonic fit */
  static my_boolean fit_use_svd = FALSE, fit_be_verbose = FALSE;
  int i,npara;
  npara = 2 + nharm * 2;
  if(nharm > MAX_NR_HARM){
    fprintf(stderr,"analyze_splitting_azi_dep: error, too many harmonics (%i vs %i)\n",
	    nharm,MAX_NR_HARM);
    exit(-1);
  }
  /* local */
  my_vecalloc(&dfazi,n_azi,"sav2splitting");
  my_vecalloc(&dfazi_sig,n_azi,"sav2splitting");
  my_vecalloc(&dt_copy,n_azi,"sav2splitting");
  my_vecalloc(&rxazi,n_azi,"sav2splitting");
  my_vecalloc(&harm_siga,npara,"sav2splitting");
  for(i=0;i < n_azi;i++)		/* get x axis values in radians */
    rxazi[i] = DEG2RAD(xazi[i]);
  /* make a copy of the dt vector */
  a_equals_b_vector(dt_copy,dt,n_azi);
  /* 
     compute mean azimuth and difference from mean azimuth,
     differences will be in dfazi
  */
  diff_from_orient_mean(fazi,n_azi,dfazi,mean_fazi,dt,
			use_dt_for_avg,&mean_vec_dt);
  /* 
     compute  stats for dt and deviations 
  */
  /* 
     deviation from mean fast azimuth (note that the aithmetic mean
     need not be exactly zero, since this is an orientational quantity
  */
  stat_moments(dfazi,n_azi,dfazi_stats,(dfazi_stats+1),(dfazi_stats+2),
	       (dfazi_stats+3),(dfazi_stats+4),(dfazi_stats+5),FALSE);
  /* variations in delay time  */
  stat_moments(dt_copy,n_azi,dt_stats,(dt_stats+1),(dt_stats+2),
	       (dt_stats+3),(dt_stats+4),(dt_stats+5),FALSE);
  /* assign to those quantities that we actually want */
  *std_fazi = dfazi_stats[2];
  if(use_dt_for_avg){
    /* vector average */
    *mean_dt = mean_vec_dt;
  }else{
    *mean_dt = dt_stats[0]; 
  }
  *std_dt = dt_stats[2];
  /* remove the mean dt */
  for(i=0;i < n_azi;i++)
    dt_copy[i] -= (*mean_dt);
  /* 
     scale fluctuations in fast azimuth by std 
  */
  if(dfazi_stats[2] > EPS_COMP_PREC){
    scale_vector(dfazi,1/dfazi_stats[2],n_azi);
  }else{
    fprintf(stderr,"analyze_splitting_azi_dep: fast azimuth constant!\n");
  }
  /* 
     scale reduce fluctuations in dt by std 
  */
  if(dt_stats[2] > EPS_COMP_PREC){
    scale_vector(dt_copy,1/dt_stats[2],n_azi);
  }else{
    fprintf(stderr,"analyze_splitting_azi_dep: splitting time constant!\n");
  }
  /* 
     fit harmonics to the normalized deviations from the mean
     fast azimuth
  */
  if(n_azi < 3){
    fprintf(stderr,"analyze_splitting_azi_dep: n_azi: %i < 3, no harmonic analysis\n",
	    n_azi);
    nharm = 0;
  }

  if(nharm){
    fit_harmonic(rxazi,dfazi,dfazi_sig,n_azi,nharm,
		 fit_svd_cutoff,
		 a_fazi,harm_siga,
		 &chi2,FALSE,fit_use_svd,fit_be_verbose);
    /* assign the strength of harmonic terms of fazi variation */
    for(i=0;i < nharm;i++){	
      ma_fazi[i]  = ONE_OVER_SQRT_TWO * hypot(a_fazi[(i+1)*2],a_fazi[(i+1)*2+1]);
      sma_fazi[i] = ONE_OVER_SQRT_TWO * hypot(harm_siga[(i+1)*2],harm_siga[(i+1)*2+1]);
    }
    /* 
     fit to normalized dt variations (note that harm_siga gets reused)
    */
    fit_harmonic(rxazi,dt_copy,dfazi_sig,n_azi,nharm,
		 fit_svd_cutoff,
		 a_dt,harm_siga,
		 &chi2,FALSE,fit_use_svd,fit_be_verbose);
    /* strength of harmonic terms */
    for(i=0;i < nharm;i++){	  
      ma_dt[i] = ONE_OVER_SQRT_TWO *  hypot(a_dt[(i+1)*2],a_dt[(i+1)*2+1]);
      sma_dt[i] = ONE_OVER_SQRT_TWO * hypot(harm_siga[(i+1)*2],harm_siga[(i+1)*2+1]);
    }
  }
  free(dfazi_sig);free(dfazi);free(rxazi);free(dt_copy);free(harm_siga);
}
