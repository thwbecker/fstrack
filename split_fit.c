/* 

reads in an sav file and tries to match a splitting observation


*/
#include "fstrack.h"
void calc_skfit_misfit(COMP_PRECISION, COMP_PRECISION ,COMP_PRECISION ,COMP_PRECISION ,
		       COMP_PRECISION, COMP_PRECISION, COMP_PRECISION *);



/* 

 
*/
int main(int argc, char **argv)
{
  COMP_PRECISION sav[36],savr[36],cijkl[81],x[3],incidence,azi,
    *fazi,*vsphase,tmp,misfit[3],
    ciso[36],chex[36],ctet[36],cort[36],cmon[36],ctri[36],k,g,sav_scc[36],sav_ast[36],sav_lith[36],
    scc_irot[9],*dt,*sfastx,*sfasty,dinc,dazi,mean_fazi,*xazi,std_fazi,ti_azi,vel[9],
    mean_dt,std_dt,best_la,best_alpha,best_beta,best_la_b0,best_alpha_b0,min_misfit,min_misfit_b0,
    symm_frac[6],ti_vec[6],tiamp,vs1,vs2,vsmean,best_fazi,best_fazi_b0,best_dt, best_dt_b0,
    alpha,beta,gamma,alpha_max,orig_fazi,orig_dt,*a_fazi,*a_dt,*ma_fazi,*msa_fazi,
    *ma_dt,*msa_dt;
  int i,j,tioff,ibc,ilc,iac;
  int n_azi,n_inc,nharm=0,npara,debug=0;	
  COMP_PRECISION azi_min,azi_max,inc_min,inc_max,inc_def,fit_fazi,fit_dt;
  COMP_PRECISION rho, layerd,layerd_total,layer_add;
  const int use_dt_for_avg = FALSE; /* use dt for weighting averages? */

  FILE *in;
  int scca_old_mode,hexsys_rotate;
  /* 
     defaults
  */
  alpha = beta = gamma = 0.0;	/* rotation angles */
  inc_min = 10.0;
  /* initial */
  layerd = 190.0;		/* layer thickness in [km] */
  rho = 3.40434;			/* density in [g/cm^3] */

  nharm = 0;			/* harmonic degree */
  hexsys_rotate = 0;		/* rotate lithospheric tensor into hex system */
  scca_old_mode = 0;		/* SCC alignment mode */
  inc_def=5;		/* default incidence angle */

  /* command line arguments */
  if(argc<5){
    fprintf(stderr,"%s: sav_asth sav_lith fazi[deg] dt [debug, %i]\n",argv[0],debug);
    exit(-1);
  }

  sscanf(argv[3],FLT_FORMAT,&fit_fazi);
  sscanf(argv[4],FLT_FORMAT,&fit_dt);
  if(argc>5)
    sscanf(argv[5],"%i",&debug);
  azi_min = 0.0;azi_max = 360.0;n_azi = 181;
  /* scan through azimuths at constant incidence */
  inc_min = inc_def;inc_max = inc_min;n_inc = 1;

  /* 
     prepare sweep variables
  */
  dazi = (n_azi>2)?((azi_max-azi_min)/(n_azi-1)):(0.0);
  dinc = (n_inc>2)?((inc_max-inc_min)/(n_inc-1)):(0.0);

  npara = 2 + nharm * 2;

  fprintf(stderr,"%s: rho: %g layerd: %g - fit: fazi %g dt %g\n",
	  argv[0],rho,layerd,fit_fazi,fit_dt);
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
  fprintf(stderr,"%s: reading %s for asthenospheric SAV\n",argv[0],argv[1]);
  in = myopen(argv[1],"r",argv[0]);
  fscanf(in,THREE_FLT_FORMAT,x,(x+1),(x+2));
  if(!read_sym_6by6(sav_ast,in)){
    fprintf(stderr,"%s: tensor read error\n",argv[0]);
    exit(-1);
  }
  fclose(in);

  fprintf(stderr,"%s: reading %s for lithospheric SAV\n",argv[0],argv[2]);
  in = myopen(argv[2],"r",argv[0]);fscanf(in,THREE_FLT_FORMAT,&tmp,&tmp,&tmp);
  if(!read_sym_6by6(sav_lith,in)){
    fprintf(stderr,"%s: tensor read error\n",argv[0]);
    exit(-1);
  }
  fclose(in);
  /*  */
  if(hexsys_rotate){
    drex_decsym(sav_lith,&k,&g,vel,symm_frac,ti_vec,
		ciso,chex,ctet,cort,cmon,ctri,sav_scc,scc_irot,
		&scca_old_mode);
    /* 
       
       rotate into hex sccs system before analysis 

    */
    /* rotate the SCC tensor such that the fast axes is
       "North-South" and not up */
    alpha = 90.0;beta = 90.0;gamma = 0.0;
    drex_rotate_6x6_deg_ftrn(sav_scc,savr,&alpha,&beta,&gamma);a_equals_b_vector(sav_scc,savr,36);
    /* assign rotated to original for later */
    a_equals_b_vector(sav_lith,sav_scc,36);
  }


  min_misfit = min_misfit_b0 = 1e20;
  for(ibc=0,beta = 0;beta <= 90.0001; beta += 5,ibc++){ /* dip scan, leave out here */
    for(ilc=0,layer_add = 0;layer_add <= 250.01;layer_add+=5,ilc++){
      /*  */
      layerd_total = layer_add + layerd;
      if(ilc==0)
	alpha_max =1;
      else
	alpha_max = 180;
      for(iac=0,alpha=0;alpha < alpha_max;alpha+=5,iac++){
      
	//fprintf(stderr,"%s: d_asth: %g d_add: %g d_total: %g alpha: %g\n",argv[0],layerd,layer_add,layerd_total,alpha);
	
	a_equals_b_vector(sav,sav_ast,36);
	scale_vector(sav,layerd/layerd_total,36);
	/* lithospheric part */
	/* angles are in degress */
	drex_rotate_6x6_deg_ftrn(sav_lith,savr,&alpha,&beta,&gamma);
	scale_vector(savr,layer_add/layerd_total,36);
	add_a_to_b_vector(savr, sav,36);
	
	/* compute the best-fit transverse isotropy axis from the Sav
	   tensor and simple anisotropy metrics */
	drex_decsym(sav,&k,&g,vel,symm_frac,ti_vec,
		    ciso,chex,ctet,cort,cmon,ctri,sav_scc,scc_irot,
		    &scca_old_mode);
	
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
	  //fprintf(stderr,"%s: WARNING: slow axes is best-fit hex axes\n",argv[0]);
	  tiamp = -tiamp;
	  tioff = 3;
	}
	/* this returns the cartesian vector, i.e. r-t-p is z-x-y where
	   x=South,y=East,z=up */
	ti_azi = RAD2DEG(atan2(ti_vec[tioff+FSTRACK_Y],-ti_vec[tioff+FSTRACK_X]));
	fix_deg_angle(&ti_azi);
	if(ti_azi>180)
	  ti_azi -= 180;
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
					      &layerd_total,(sfastx+j),(sfasty+j),
					      (fazi+j),(dt+j),(vsphase+j*2));
	  }
	  
	  /* 
	     compute statistics 
	  */
	  analyze_splitting_azi_dep(n_azi,nharm,xazi,fazi,dt,&mean_fazi,&std_fazi,
				    &mean_dt,&std_dt,a_fazi,a_dt,ma_fazi,
				    msa_fazi,ma_dt,msa_dt,FALSE,use_dt_for_avg);
	  if((iac == 0)&&(ilc ==0)){
	    /* zero addition, original split */
	    orig_fazi = mean_fazi;
	    orig_dt = mean_dt;
	  }
	  
	  /* standard output:
	     
	     mean_fazi std_fazi mean_dt std_dt ti^h_azi ti^h_amp
	     
	  */
	  //printf("%6.2f %6.2f\t %6.2f %6.4f\t%6.2f %6.4f\n",mean_fazi,std_fazi,mean_dt,std_dt,ti_azi,sqrt(1.0-ti_vec[tioff+FSTRACK_Z]*ti_vec[tioff+FSTRACK_Z])*tiamp);
	  calc_skfit_misfit(fit_fazi,fit_dt,mean_fazi,std_fazi,mean_dt,std_dt,misfit);
	  if(debug){
	    printf("%g %g %g \t%e %e %e\n",layer_add,alpha,beta,misfit[0],misfit[1],misfit[2]);
	  }else{
	    if(ibc==0){
	      if(min_misfit_b0 > misfit[2]){
		min_misfit_b0 = misfit[2];
		best_alpha_b0 = alpha;
		best_la_b0 = layer_add;
		best_fazi_b0 = fit_fazi;
		best_dt_b0 = fit_dt;
	      }
	    }
	    if(min_misfit > misfit[2]){
	      min_misfit = misfit[2];
	      best_alpha = alpha;
	      best_la = layer_add;
	      best_beta = beta;
	      best_fazi = fit_fazi;
	      best_dt = fit_dt;
	    }
	  }
	} /* end incidence loop */
      }	/* end beta */
    }	/* end alpha */
  }	/* end layer */
  printf("best: %8.3e at a %03.1f b %02.1f ld %03.1f with %5.1f %4.2f\tbest_b0: %8.3e at a %03.1f ld %02.1f with %5.1f %4.2f\tobserved: %5.1f %4.2f asthenosphere:  %5.1f %4.2f\n",
	 min_misfit, best_alpha, best_beta, best_la,best_fazi,best_dt,
	 min_misfit_b0, best_alpha_b0, best_la_b0,best_fazi_b0,best_dt_b0,
	 fit_fazi,fit_dt,
	 orig_fazi,orig_dt);
  
  free(dt);free(fazi);free(sfastx);free(sfasty);free(vsphase);
  free(xazi);free(a_fazi);free(a_dt);
  free(ma_fazi);free(msa_fazi);
  free(ma_dt);free(msa_dt);
}


void calc_skfit_misfit(COMP_PRECISION fit_fazi, COMP_PRECISION fit_dt,
		       COMP_PRECISION mean_fazi,COMP_PRECISION std_fazi,
		       COMP_PRECISION mean_dt,COMP_PRECISION std_dt,
		       COMP_PRECISION *misfit)
{
  COMP_PRECISION azi_scale,dt_scale;

  azi_scale = (std_fazi > 10)?(std_fazi):(10);
  dt_scale = (std_dt > 0.2)?(std_dt):(0.2);

  misfit[0] = fabs(fit_fazi - mean_fazi);
  if(misfit[0] > 180)
    misfit[0] -= 180;
  if(misfit[0] > 90)
    misfit[0] = 180 - misfit[0];
  misfit[0] /= azi_scale;
  
  misfit[1]  = fabs(fit_dt - mean_dt)/dt_scale;

  misfit[2] = misfit[0] + misfit[1];
}
