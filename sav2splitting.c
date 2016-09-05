/* 

convert Sav stiffness tensor format to splitting estimates

sav2splitting sav.file output_format beta_dip layer_d density alpha_rot gamma_rot hexsys_rotation ssc_old_mode inc_def

input format:

lon lat depth upper_triangl_C_ij (j fast), i.e.
			  
SAV_11 SAV_12 SAV_13 ... SAV_16 \
       SAV_22 SAV_23 ... SAV_26 \
       ... 
                         SAV_66


this version uses a single layer

output: 

... depends


parameters:

densities: 3.35  for olivine, approx
           3.353 for olivine/enstatite(30%)

$Id: sav2splitting.c,v 1.11 2010/12/31 18:59:17 becker Exp $


*/
#include "fstrack.h"

#include "prem.h"

#define SPLITTING_VERA_FORMAT 0
#define SPLITTING_STATS 1
#define SPLITTING_STATS_SCAN 2
#define SPLITTING_BEST_FIT 3
#define SPLITTING_SCAN 4
#define SPLITTING_FULL_SCAN 5
/* 

usage:  
 

          sav2splitting sav.dat [mode, SPLITTING_STATS] [beta, 0]...\
	             [layerd, 50] [rho, 3.353] [alpha,0] [gamma, 0] ...\
		     [rotate_into_hex, 0] [scca_old_mode, 1] [inc_def, 5]
		     
		     if rho is negative, will look for prem density from z value,
		     within reason


		     beta = dip
		     (alpha, beta, gamma are Euler angles as in Dahlen and Tromp)

 */
int main(int argc, char **argv)
{
  COMP_PRECISION sav[36],savr[36],cijkl[81],x[3],incidence,azi,
    *fazi,*vsphase,
    ciso[36],chex[36],ctet[36],cort[36],cmon[36],ctri[36],k,g,sav_scc[36],
    scc_irot[9],*dt,*sfastx,*sfasty,dinc,dazi,mean_fazi,*xazi,std_fazi,ti_azi,vel[9],
    mean_dt,std_dt,
    razi,symm_frac[6],ti_vec[6],tiamp,vs1,vs2,vsmean,
    alpha,beta,gamma,
    /* fitting parameters */
    *a_fazi,*a_dt,
    /* magnitude of the fitting terms */
    *ma_fazi,*msa_fazi,
    *ma_dt,*msa_dt;
  int i,nr,n,j,iop=FTRN_STDOUT,output_format,npara,tioff;
  int n_azi,n_inc,nharm;	
  static COMP_PRECISION azi_min,azi_max,inc_min,inc_max,inc_def;
  static COMP_PRECISION rho, layerd;
  const my_boolean verbose = FALSE;
  const int use_dt_for_avg = FALSE; /* use dt for weighting averages? */
  struct prem_model prem[1];
  my_boolean rotate, hexsys_rotate, prem_dens;
  FILE *in;
  int scca_old_mode;
  /* 
     defaults
  */
  alpha = beta = gamma = 0.0;	/* rotation angles */
  output_format = SPLITTING_STATS;
  inc_min = 10.0;
  layerd = 50.0;		/* layer thickness in [km] */
  rho = 3.353;			/* density in [g/cm^3] */
  hexsys_rotate = FALSE;	/* do not rotate into hex system */
  rotate = FALSE;
  prem_dens = FALSE;
  nharm = 6;			/* harmonic degree */
  scca_old_mode = 1;	/* default is old mode */
  inc_def=5;		/* default incidence angle */

  /* command line arguments */
  if((argc > 1)&&(strcmp(argv[1],"stdin")!=0)){
    in = myopen(argv[1],"r",argv[0]);
    fprintf(stderr,"%s: reading tensor from %s\n",argv[0],argv[1]);
  }else{
    in = stdin;
    fprintf(stderr,"%s: reading tensor from %s\n",argv[0],"stdin");
  }
  if(argc > 2)
    sscanf(argv[2],"%i",&output_format);
  if(argc > 3)
    sscanf(argv[3],"%lf",&beta);
  if(argc > 4)
    sscanf(argv[4],"%lf",&layerd);
  if(argc > 5)
    sscanf(argv[5],"%lf",&rho);
  if(argc > 6){
    sscanf(argv[6],"%lf",&alpha);
  }
  if(argc > 7){
    sscanf(argv[7],"%lf",&gamma);
  }
  if(argc > 8){
    sscanf(argv[8],"%i",&i);	/* rotate into hex symmetry system */
    if(i)
      hexsys_rotate = TRUE;
    else
      hexsys_rotate = FALSE;
  }
  if(argc>9)			
    sscanf(argv[9],"%i",&scca_old_mode);
  if(argc>10)			
    sscanf(argv[10],"%lf",&inc_def);
 
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
  if(rotate && hexsys_rotate){
    fprintf(stderr,"%s: error: both regular rotation and hex system rotation chosen\n",
	    argv[0]);
    exit(-1);
  }
  if(rho < 0){			/* get density from PREM */
    prem_dens = TRUE;
    fprintf(stderr,"%s: WARNING: varying density with z depth and PREM\n",argv[0]);
    prem_read_model(PREM_MODEL_FILE,prem,FALSE);
  }
  if(scca_old_mode)
    fprintf(stderr,"%s: using old mode of SCCA alignment\n",argv[0]);
  else
    fprintf(stderr,"%s: using new mode of SCCA alignment\n",argv[0]);
  
  /* 
     parameters
  */

  /* 
     min/max ranges and number of samples for azimuth sweep 
  */
  azi_min = 0.0;azi_max = 360.0;n_azi = 181;
  switch(output_format){
  case SPLITTING_VERA_FORMAT:	/* 
				   scan through azimuths 
				   and indices and print output
				   as in Vera's old 
				*/
    inc_min = 3.0;inc_max = 15.0;n_inc = 13;break;
  case SPLITTING_SCAN:
    /* scan through azimuths at constant incidence */
    inc_min = inc_def;inc_max = inc_min;n_inc = 1;break;
  case SPLITTING_FULL_SCAN:
    /* fine scan through azimuths and incidence */
    n_azi = 721;		/*  */
    //inc_min = 0.25;inc_max = 89.75;n_inc = 91;
    inc_min = 0;inc_max = 90;n_inc = 181;
    break;
   case SPLITTING_STATS:
    /* statistics for azimuthal dependence */
    inc_min = inc_def;inc_max = inc_min;n_inc = 1;break;
  case SPLITTING_STATS_SCAN:
    /* scan through indicence  */
    inc_min = 3;inc_max = 15.0;n_inc = 5;break;
  case SPLITTING_BEST_FIT:
    break;
  default:
    inc_min = inc_def;inc_max = inc_min;n_inc = 1;break;
  }
  /* 
     prepare sweep variables
  */
  dazi = (n_azi>2)?((azi_max-azi_min)/(n_azi-1)):(0.0);
  dinc = (n_inc>2)?((inc_max-inc_min)/(n_inc-1)):(0.0);


  npara = 2 + nharm * 2;

  fprintf(stderr,"%s: dip: %g alpha: %g gamma: %g rho: %g d: %g nharm: %i npara: %i nazi: %i inc_def: %g\n",
	  argv[0],beta,alpha,gamma,rho,layerd,nharm,npara,n_azi,inc_def);
  /* 
     allocate for fitting 
  */

  my_vecalloc(&a_fazi,  npara,"sav2splitting");
  my_vecalloc(&a_dt,    npara,"sav2splitting");
  my_vecalloc(&ma_fazi, nharm,"sav2splitting");
  my_vecalloc(&msa_fazi,nharm,"sav2splitting");
  my_vecalloc(&ma_dt,   nharm,"sav2splitting");
  my_vecalloc(&msa_dt,  nharm,"sav2splitting");
  
  /* allocate */
  my_vecalloc(&xazi,n_azi,"sav2splitting");
  my_vecalloc(&fazi,n_azi,"sav2splitting");
  my_vecalloc(&sfastx,n_azi,"sav2splitting");
  my_vecalloc(&sfasty,n_azi,"sav2splitting");
  my_vecalloc(&dt,n_azi,"sav2splitting");
  my_vecalloc(&vsphase,n_azi*2,"sav2splitting");
  /* 
     azimuth values
  */
  for(i=0,azi=azi_min;i < n_azi;i++,azi += dazi)
    xazi[i] = azi;
  /* 
     input 
  */
  n=0;
  while(fscanf(in,THREE_FLT_FORMAT,x,(x+1),(x+2))==3){
    /* 

    read in 6,6 matrix

    */
    if(!read_sym_6by6(sav,in)){
      fprintf(stderr,"%s: tensor read error\n",argv[0]);
      exit(-1);
    }

    if(rotate){		/* rotate tensor out of horizontal */
      fprintf(stderr,"%s: WARNING: rotating tensor by %g %g %g (deg)\n",argv[0],alpha,beta,gamma);
      /* angles are in degress */
      drex_rotate_6x6_deg_ftrn(sav,savr,&alpha,&beta,&gamma);
      a_equals_b_vector(sav,savr,36);
      zero_small_entries(sav,36);
    }

    /* compute the best-fit transverse isotropy axis from the Sav tensor */
    drex_decsym(sav,&k,&g,vel,symm_frac,ti_vec,
		ciso,chex,ctet,cort,cmon,ctri,sav_scc,scc_irot,
		&scca_old_mode);
    if(hexsys_rotate){
      /* 
	 
      rotate into hex sccs system before analysis 

      */
      fprintf(stderr,"%s: WARNING: rotating into best-fit hex system\n",argv[0]);
      /* rotate the SCC tensor such that the fast axes is
	 "North-South" and not up */
      alpha = 90.0;beta = 90.0;gamma = 0.0;
      drex_rotate_6x6_deg_ftrn(sav_scc,savr,&alpha,&beta,&gamma);
      a_equals_b_vector(sav_scc,savr,36);
      zero_small_entries(sav_scc,36);
      /* assign rotated to original for later */
      a_equals_b_vector(sav,sav_scc,36);
      /* redo decomposition for TI axes */
      drex_decsym(sav_scc,&k,&g,vel,symm_frac,ti_vec,
		  ciso,chex,ctet,cort,cmon,ctri,sav_scc,scc_irot,
		  &scca_old_mode);
    }
    /* slow and fast S from hexagonal fit */
    vs1 = vel[7]/sqrt(rho);
    vs2 = vel[8]/sqrt(rho);
    vsmean = (vs1+vs2)/2.0;
    /* anisotropy in percent */
    tiamp = (vs1-vs2)/vsmean*100.0;
    /* 
       select the fast axis, best-fitting might be slow
    */
    if(tiamp > 0)
      tioff = 0;
    else{
      fprintf(stderr,"%s: WARNING: slow axes is best-fit hex axes\n",argv[0]);
      tiamp = -tiamp;
      tioff = 3;
    }
    
    /* this returns the cartesian vector, i.e. r-t-p is z-x-y where
       x=South,y=East,z=up */
    ti_azi = RAD2DEG(atan2(ti_vec[tioff+FSTRACK_Y],-ti_vec[tioff+FSTRACK_X]));
    fix_deg_angle(&ti_azi);
    if(ti_azi>180)
      ti_azi -= 180;
    /* 
       convert tensor to normalized cijkl format 
    */
    if(prem_dens){
      prem_get_rho(&rho,1.0-x[2]/6371.0,prem);
      rho /= 1000.0;
      if(rho < 2.6)
	rho=2.6;
      fprintf(stderr,"%s: using density %g at depth %g\n",argv[0],rho,x[2]);
    }
    vera_sav_to_cijkl_ftrn(sav,&rho,cijkl);
    /* 
       sweep through incidence
    */

    for(i=0,incidence=inc_min;i<n_inc;i++,incidence += dinc){ /* begin inc loop */
      /* 
	 sweep through azimuths  
      */
      for(j=0;j < n_azi;j++){
	/* compute the splitting parameters for a single layer,
	   azimuth, and incidence */
	vera_layer_split_from_tensor_ftrn(cijkl,&incidence,(xazi+j),
					  &layerd,(sfastx+j),(sfasty+j),
					  (fazi+j),(dt+j),(vsphase+j*2));
      }
      
      /* 
	 process azimuthal dependence
      */
      /* 
	 compute statistics 
      */
      analyze_splitting_azi_dep(n_azi,nharm,xazi,fazi,dt,&mean_fazi,&std_fazi,
				&mean_dt,&std_dt,a_fazi,a_dt,ma_fazi,
				msa_fazi,ma_dt,msa_dt,TRUE,use_dt_for_avg);
      if(verbose){
	/* 
	   print mean fast azimuth and splitting time as well as std. 
	*/
	fprintf(stderr,"%s: inc: %5.1f fast_azi: %6.1f (+/-%5.1f%% of 180) dt: %6.2fs (+/-%5.1f%%)\n",
		argv[0],incidence,mean_fazi,std_fazi/180.*100.,
		mean_dt,std_dt/mean_dt*100.);
	/* harmonic terms */
	fprintf(stderr,"%s: fit to dazi: ",argv[0]);
	for(j=0;j < nharm;j++)	/* normalized strength of terms and
				   relative to uncertainties */
	  fprintf(stderr,"%1it: %6.2f (%6.2f) ",
		  j+1,ma_fazi[j],log(ma_fazi[j]/msa_fazi[j]));
	fprintf(stderr,"\n");
	fprintf(stderr,"%s: fit to   dt: ",argv[0]);
	for(j=0;j < nharm;j++)
	  fprintf(stderr,"%1it: %6.2f (%6.2f) ",
		  j+1,ma_dt[j],log(ma_dt[j]/msa_dt[j]));
	fprintf(stderr,"\n");
      }
      switch(output_format){
      case SPLITTING_FULL_SCAN:
	/* print in 
	   
	   back-azimuth incidence fast_azi split_time vs1 vs2 

	   format */
	for(j=0;j < n_azi;j++){
	  fprintf(stdout,"%6.1f %6.3f\t%8.3f %8.5f\t%9.5f %9.5f\n",
		  xazi[j],incidence,fazi[j],dt[j],vsphase[j*2],vsphase[j*2+1]);
	}
	break;
      case SPLITTING_SCAN:
      case SPLITTING_VERA_FORMAT:
	/* print */
	for(j=0;j < n_azi;j++){
	  vera_print_splitting_ftrn(&incidence,(xazi+j),(sfastx+j),(sfasty+j),(fazi+j),(dt+j),&iop);    
	}
	break;
      case SPLITTING_STATS_SCAN:
      case SPLITTING_STATS:
	/* standard output:

	lon lat incidence mean_fazi std_fazi mean_dt std_dt af_t af_t2 ... ti^h_azi ti^h_amp
	
	*/
	printf("%6.2f %6.2f %6.2f\t%6.2f %6.2f\t %6.2f %6.4f\t",
	       x[0],x[1],incidence,mean_fazi,std_fazi,mean_dt,std_dt);
	for(j=0;j < nharm;j++)printf("%8.4f ",ma_fazi[j]);
	printf("\t%6.2f %6.4f\n",ti_azi,
	       sqrt(1.0-ti_vec[tioff+FSTRACK_Z]*ti_vec[tioff+FSTRACK_Z])*tiamp);
	//printf("\n");
	break;
      case SPLITTING_BEST_FIT:
	/* 
	   print out the variations along side the best-fit
	*/
	
	for(j=0;j < n_azi;j++){
	  razi = DEG2RAD(xazi[j]); /* azimuth in radians */
	  printf("%g\t%g\t\t%g %g\t\t%g %g\n",
		 incidence,xazi[j],fazi[j],
		 mean_fazi + std_fazi * 
		 evaluate_model(razi,npara,a_fazi, 
				(void (*)(void))(harm_fit_func),
				nharm),
		 dt[j],mean_dt + std_dt * evaluate_model(razi,npara,a_dt, 
				(void (*)(void))(harm_fit_func),
				nharm));
	}
	break;
      default:
	fprintf(stderr,"%s: output mode %i undefined\n",argv[0],output_format);
	exit(-1);
      }
    } /* end incidence loop */
    n++;
  }
  free(dt);free(fazi);free(sfastx);free(sfasty);free(vsphase);
  free(xazi);free(a_fazi);free(a_dt);
  free(ma_fazi);free(msa_fazi);
  free(ma_dt);free(msa_dt);
  if(verbose)
    fprintf(stderr,"%s: converted %i tensors using rho: %g g/cm^3 layerd: %g km\n",
	    argv[0],n,rho,layerd);
}
  
