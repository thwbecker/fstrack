/* 
   
Runge Kutta driver routine for forward integration of 
y(nvar) state from x1 to x2 using accuracy eps,
starting step h1, and minimum step him

all derivative input parameters are passed in the dp structure


return the new state y[] and 
the number of good and bad steps (nok and nbad) 


modified from numerical recipes

WARNING: call all these routines have been changed to C-style
y(0...n-1)

returns 0 on sucecss, else error codes

have to be initialized before:
dp->kmax   is the number of intermediate steps to be saved, should be >=0
dp->dxsave is the approximate time interval at which to save intermediate steps


this routine assumes that x[0..2] are locations, x[3...11] are strains,
and the rest is dealt with by DREX

$Id: odeint.c,v 1.9 2010/12/31 18:59:17 becker Exp $

*/
#include "fstrack_flow.h"
#include "nr_defines.h"

#define RK_MAXSTP 50000


int nr_odeint(COMP_PRECISION *ystart, /* start vector */
	      int nvar,	/* number of variables */
	      COMP_PRECISION x1, /* start time */
	      COMP_PRECISION x2, /* stop time */
	      COMP_PRECISION eps, /* accuracy */
	      COMP_PRECISION h1, /* starting timestep */
	      COMP_PRECISION hmin, /* min timestep */
	      COMP_PRECISION hmax, /* max timstep, gets passed 
				      to nr_rkqs */
	      int *nok,int *nbad, /* nr of good and bad step */
	      struct der_par *dp /* parameters for derivatives */
	      )
{
  /* 
     
  regular RK stuff 
  
  */
  static my_boolean init = FALSE;
  static int oldn = 0;
  int nstp,i,o1,error_code,ilim;
  COMP_PRECISION xsav=0.0,x,hnext,hdid,h;
  COMP_PRECISION *yscal,*y,*dydx;
  my_boolean lya_exit;
  /* 
     
  additional stuff 

  */
  /* 

  Lyapunov 

  */
  static COMP_PRECISION tlalphaold[3];
  int icntbailout = 0;
  my_boolean ifrenorm = FALSE;
  COMP_PRECISION alpha[3],slalpha,tlalpha;
  /* 

  start code 

  */
  if((!init)||(nvar > oldn)){	
    /* 

    first call, or larger temporary storage needed:
    
    initialize the helper arrays for intermediate storage

    */
    if(nvar > 12)
      PE("nr_odeint: attempting texture computation");
    if(dp->kmax > 1000)
      PEE("nr_odeint: really more than 1000 intermediate results to be stored? exiting...");
    if(dp->kmax < 0)
      PEE("nr_oderint: dp->kmax should be >= 0");
    if(dp->kmax){
      /* 
	 we do want intermediate results to be stored
      */
      if(init){
	/* 
	   already allocated, free yp matrix, 
	   nvar might have changed 
	*/
	free(dp->yp);
      }else{
	/* only allocate xp once */
	my_vecalloc(&dp->xp,dp->kmax,"nr_odeint: 1");
      }
      my_vecalloc(&dp->yp,nvar*dp->kmax,"nr_odeint: 2");
    }
    oldn = nvar;
    /* 
       Lyapunov calculation? 
    */
    if(dp->calc_lyapunov){
      if(dp->rlbailout <= 0.0){
	fprintf(stderr,"nr_odeint: error, calc_lyapnuov is on, but rlbailout is %g\n",
		dp->rlbailout);
	exit(-1);
      }
      for(i=0;i < 3;i++)
	tlalphaold[i] = -1.0;
    }
    /* 
       DREX
    */
    if((nvar > 12)&&(dp->texture_mode != NO_TEXTURE)){
      drex_test_if_initialized(dp->drex,"nr_odeint");
      if(dp->drex_epsfac < 0){
	fprintf(stderr,"nr_odeint: error: drex eps factor is negative: %g\n",
		dp->drex_epsfac);
	exit(-1);
      }
      if(dp->strain_rate_control){
	fprintf(stderr,"nr_odeint: using strain rate control. %g times 1/strain_rate is max timestep\n",
		dp->eps_strain);
      }else{
	fprintf(stderr,"nr_odeint: texture has modified RK eps, factor: %.3e (%s)\n",
		dp->drex_epsfac,
		((determine_n_error(1000,dp)==1000)?("but still checking"):("no error checking")));
      }
    }
    init = TRUE;
  }
  /* 

  allocate local vectors 

  */
  my_vecalloc(&yscal,nvar,"nr_odeint: 3");
  my_vecalloc(&y,nvar,"nr_odeint: 4");
  my_vecalloc(&dydx,nvar,"nr_odeint: 5");
  /* 

  initialize local and global arrays 

  */
  lya_exit = FALSE;
  x = x1;
  h = NR_SIGN(h1, x2 - x1);
  *nok = (*nbad) = dp->kount = 0;
  a_equals_b_vector(y,ystart,nvar); /* save start vector */
  if (dp->kmax > 0) 		/* make sure first step is stored */
    xsav = x - dp->dxsav * 2.0;
  for(nstp=0;nstp < RK_MAXSTP;nstp++) { 
    /* take at most MAXSTP steps */
    /* 
       
    start main loop

    */
    /* 
       compute y derivatives with respect to x at time x, those are needed for 
       the rkqs step
    */
    rk_derivs(x,y,dydx,nvar,dp);
    /* 
       determine the scale for each component of y[nvar] which
       will be used to determine the relative error
    */
    determine_error_scale(x,y,dydx,h,nvar,dp,yscal);
    /* 
       should we store intermediate results? 
    */
    if ((dp->kmax > 0) && 
	(dp->kount < dp->kmax - 2) && 
	(fabs(x - xsav) > fabs(dp->dxsav))) {
      /* 
	 assign intermediate solution 
      */
      dp->xp[dp->kount] = x;	/* time of storage */
      o1 = dp->kount*nvar;
      /* solution vector storage */
      a_equals_b_vector((dp->yp+o1),y,nvar);
      dp->kount++;
      xsav = x;
    }
    if ((x+h-x2)*(x+h-x1) > 0.0) /* if stepsize overshot, decrease */
      h = x2-x;
    /* 
       advance by one RK step, this routine returns the actual
       step size taken, hdid, and the next suggested one, hnext
    */
    error_code = 
      nr_rkqs(y,dydx,nvar,&x,h,hmax,eps,yscal,&hdid,&hnext,dp);
    if(error_code != 0)		/* return with error code */
      return odeint_exit(&y,&dydx,&yscal,x,dp,error_code);
    /* 

       Lyapunov part 

    */
    if(dp->calc_lyapunov){
      /* 
	 
      check if one of the abs(components) of the strain tensor are
      larger than renorm
      
      */
      ifrenorm = FALSE;
      ilim = MIN(nvar,12);
      for(i=3;i < ilim;i++){
	if(fabs(y[i]) > dp->renorm){
	  ifrenorm = TRUE;
	  break;
	}
      }
      if(ifrenorm){
	//
	//     renormalize the deformation matrix, call with the whole y[12]
	//     vector
	//
	gramschmidt(y,alpha);
	//
	//     and increment lyapunov exponents
	//
	icntbailout = 0;
	for(i=0;i < 3;i++){	/* check all components */
	  slalpha = log(alpha[i]);
	  /* 
	     increment the exponent
	  */
	  dp->rlyap[i] += slalpha;
	  /*

	    check for bailout

	  */
	  if(fabs(dp->rlyap[i])>EPS_PREC){
	    // 
	    // compute the fractional change in the Lyapunov sum
	    //
	    tlalpha = fabs(slalpha/dp->rlyap[i]);
	    //
	    //     need at least one renorm step (tlalphaold() init as -1)
	    //     and a decrease in the Delta addition
	    //
	    if((tlalphaold[i] > tlalpha) && /* increment to exponent needs to 
					       decrease */
	       (fabs(x) > LYA_TMIN_DEF) && /* have to advect for more than tmin */
	       (tlalpha < dp->rlbailout)) /* fractional has to be smaller than 
					     the bailout criterion  */
	      icntbailout++;
	    //     save increment for next check
	    tlalphaold[i] = tlalpha;
	  }
	} /* end i loop */
      }
      if(icntbailout == 3){
	//
	//     only bailout if all exponents < threshold
	//
	lya_exit = TRUE;
      }
    } /* end Lyapunov branch */
    /* 
       
    evaluate success/failure of step

    */
    if (fabs(hdid - h) < EPS_PREC)
      ++(*nok); 
    else 
      ++(*nbad);
    if (((x-x2)*(x2-x1) >= 0.0) || (lya_exit)) {
      /* 
	 we are done, because of t final being reached or because of lyapnuov
	 exponents didn't change by much 
      */
      a_equals_b_vector(ystart,y,nvar);
      if (dp->kmax > 0) {
	/* save final step */
	dp->xp[dp->kount]=x;
	o1 = dp->kount * nvar;
	a_equals_b_vector((dp->yp+o1),y,nvar);
	dp->kount++;
      }
      return odeint_exit(&y,&dydx,&yscal,x,dp,0);
    }
    if (fabs(hnext) < hmin){
      /* 
	 ERROR: stepsize underflow 
      */
      fprintf(stderr,"nr_odeint: stepsize %.5e smaller than minimum, %.5e\n",
	      hnext,hmin);
      return odeint_exit(&y,&dydx,&yscal,x,dp,-1);
    }
    if (fabs(hnext) > hmax){
      fprintf(stderr,"nr_odeint: internal error: hnext: %g h_max: %g\n",
	      hnext,hmax);
      exit(-1);
    }
    h = hnext;
  }
  /* 
     ERROR: too many steps 
  */
  fprintf(stderr,"nr_odeint: too many steps, nok: %i nbad: %i\n",*nok,*nbad);
  return odeint_exit(&y,&dydx,&yscal,x,dp,-2);
}

#undef RK_MAXSTP

#define NR_RK_SAFETY 0.9
#define NR_RK_PGROW -0.2
#define NR_RK_PSHRNK -0.25
#define NR_RK_ERRCON 1.889568e-4	/* = (5/NR_RK_SAFETY)**(1/NR_RK_PGROW) */

/* 

take a fifth order quality controlled Runge Kutta step

uses the error scales yscal as input, htry as the stepsize to try, and
advances y[n] by htry

makes extensive use of the routines in deriv.c

on input, dydx has to hold the derivatives of y at time x

return 0 on success returns -1 on stepsize underflow

hmax should be >= 0 and indicate the largest absolute step taken,
regardless of error criterion. if dp->strain_rate_control is
activated, hmax will be limited to the strain rate determined step
locally


*/
int nr_rkqs(COMP_PRECISION *y,COMP_PRECISION *dydx,
	    int n,COMP_PRECISION *x,COMP_PRECISION htry,
	    COMP_PRECISION hmax,COMP_PRECISION eps,
	    COMP_PRECISION *yscal,
	    COMP_PRECISION *hdid,COMP_PRECISION *hnext,
	    struct der_par *dp)
{
  int i,nerrlim;
  COMP_PRECISION errmax,h,*yerr,*ytemp,tmp,htemp,hmax_loc;
#ifdef FSTRACK_DEBUG
  int ncount=0;
  COMP_PRECISION hold;
  int nerrmax;
#endif  
  my_vecalloc(&yerr,n,"nr_rkqs: 1");
  my_vecalloc(&ytemp,n,"nr_rkqs: 2");
  if(dp->strain_rate_control){
    /* adjust the maximum timestep to the strain rate determined time step */
    tmp = calc_max_dt_from_strain_rate(*x,y,dp);
    hmax_loc = MIN(hmax,tmp);
  }else{
    hmax_loc = hmax;
  }
  /* 
     check if timestep large enough and > 0
  */
  if(hmax_loc <= eps){
    fprintf(stderr,"nr_rkqs: error: h_max: %g at eps: %g (h_max should be >=0)\n",
	    hmax_loc,eps);
    if(dp->strain_rate_control)
      fprintf(stderr,"nr_rkqs: strain rate controlled hmax: %g\n",
	      calc_max_dt_from_strain_rate(*x,y,dp));
    exit(-1);
  }
  /* 
     initial trial stepsize 
  */
  h = MIN(fabs(htry),hmax_loc) * MY_SIGN(htry);
  /*
     determine the limit dimension for the error condition
     checks
  */
  nerrlim = determine_n_error(n,dp);
  for (;;) {
    /* 
       advance by one Cash Carp Runge Kutta step from y to ytemp and
       return the truncation error, yerr
    */
#ifdef FSTRACK_DEBUG
    ncount++;
    if(0)
      fprintf(stderr,"nr_rkqs: i: %3i t: %11g dt: %11g n: %i nerr: %i dt_max: %g\n",
	      ncount,*x,h,n,nerrlim,hmax_loc);
#endif
    /* take step here */
    nr_rkck(y,dydx,n,*x,h,ytemp,yerr,dp);
    /* 
       determine error condition 
    */
    for(errmax=0.0,i=0;i < nerrlim;i++) {
      tmp  =  fabs(yerr[i]/yscal[i]);
      if(tmp > errmax){
	errmax = tmp;
#ifdef FSTRACK_DEBUG
	nerrmax = i;
#endif
      }
#ifdef FSTRACK_DEBUG
      if(yscal[i] <= EPS_COMP_PREC)
	fprintf(stderr,"nr_rkqs: yscal[%i] too small: %.5e\n",i,yscal[i]);
#endif
    }
#ifdef FSTRACK_DEBUG
    if(eps <= EPS_COMP_PREC)
      PEE("nr_rkqs: eps too small");

#endif
    /* 
       scale max error 
    */
    errmax /= eps;
    /* 
       add physical, a priori knowledge based, error conditions,
       if any
    */
    check_physics_based_error(y,dydx,&n,&eps,&errmax,&h);
    if (errmax > 1.0) {
      /* 
	 error is too large, shrink stepsize
      */
#ifdef FSTRACK_DEBUG
      hold = h;
#endif
      htemp = NR_RK_SAFETY * h * pow(errmax,NR_RK_PSHRNK);
      /* 
	 but not more than a factor of ten 
      */
      tmp = 0.1 * h;
      h = ((h >= 0.0) ? (MAX(htemp,tmp)):(MIN(htemp,tmp)));
      /* test if h is too small */
      tmp = (*x) + h;
#ifdef FSTRACK_DEBUG
      if(0){
	/* 
	   report the max error on failure 
	*/
	fprintf(stderr,"nr_rkqs: max error: %11.4e at n: %i scale %.4e: rel: %.4e eps: %.4e hold: %.3e hnew: %.3e\n",
		fabs(yerr[nerrmax]),nerrmax+1,yscal[nerrmax],errmax,eps,hold,h);
      }
#endif
      if (fabs(tmp - *x)<EPS_PREC) {
	fprintf(stderr,"nr_rkqs: stepsize underflow: t: %.5e h: %5e y[%i]:\n",
		*x,h,n);
	print_vector(y,n,stderr);
	free(ytemp);free(yerr);return(-1);
      }
      continue;
    } else {
      /* 
	 error is OK, increase stepsize but no more than a factor of
	 five, and also to not more than hmax_loc
      */
      if (errmax > NR_RK_ERRCON) 
	*hnext = NR_RK_SAFETY * h * pow(errmax,NR_RK_PGROW);
      else 
	*hnext = 5.0*h;
      if(dp->strain_rate_control){
	/* adjust the maximum step */
	tmp = calc_max_dt_from_strain_rate(*x,y,dp);
	hmax_loc = MIN(hmax,tmp);
      }
      /* 
	 limit hnext to max step
      */
      *hnext = MIN(fabs(*hnext), hmax_loc) * MY_SIGN(*hnext);
      *hdid = h;		/* save actual used timestep */
      *x += *hdid;		/* increment time */
      /* 
	 save solution 
      */
      a_equals_b_vector(y,ytemp,n);
      break;
    }
  }
  free(ytemp);
  free(yerr);
  return 0;
}
#undef NR_RK_SAFETY
#undef NR_RK_PGROW
#undef NR_RK_PSHRNK
#undef NR_RK_ERRCON

/* 

compute fifth order Cash Karp Runge Kutta step. on input, dydx has to hold 
the value of the derivatives of y at time x 

*/
void nr_rkck(COMP_PRECISION *y,COMP_PRECISION *dydx,
	     int n,COMP_PRECISION x,COMP_PRECISION h,
	     COMP_PRECISION *yout,COMP_PRECISION *yerr,
	     struct der_par *dp)
{
  int i;
  /* static factors */
  static COMP_PRECISION 
    a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
    b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
    b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
    b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
    b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
    c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
    dc5 = -277.0/14336.0;
  static COMP_PRECISION dc1=(37.0/378.0)-(2825.0/27648.0),
    dc3=(250.0/621.0)-(18575.0/48384.0),
    dc4=(125.0/594.0)-(13525.0/55296.0),
    dc6=(512.0/1771.0)-0.25;
  /* 
     real variables 
  */
  COMP_PRECISION *ak2,*ak3,*ak4,*ak5,*ak6,*ytemp;

  /* allocate arrays */
  my_vecalloc(&ak2,n,"nr_rkck: ak2");
  my_vecalloc(&ak3,n,"nr_rkck: ak3");
  my_vecalloc(&ak4,n,"nr_rkck: ak4");
  my_vecalloc(&ak5,n,"nr_rkck: ak5");
  my_vecalloc(&ak6,n,"nr_rkck: ak6");
  my_vecalloc(&ytemp,n,"nr_rkck: ytemp");
  /* 
     compute Runge Kutta formulae 
  */
  for (i=0;i<n;i++)
    ytemp[i] = y[i] + b21*h*dydx[i];
  rk_check_phys_limit(ytemp,n,dp,TRUE);
  
  /* at t + a2*h */
  rk_derivs(x+a2*h,ytemp,ak2,n,dp);
  for (i=0;i<n;i++)
    ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
  rk_check_phys_limit(ytemp,n,dp,TRUE);

  /* at t+a3*h*/
  rk_derivs(x+a3*h,ytemp,ak3,n,dp);
  for (i=0;i<n;i++)
    ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
  rk_check_phys_limit(ytemp,n,dp,TRUE);

  /*  at t+a4*h */
  rk_derivs(x+a4*h,ytemp,ak4,n,dp);
  for (i=0;i<n;i++)
    ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
  rk_check_phys_limit(ytemp,n,dp,TRUE);

  /* at t+a5*h */
  rk_derivs(x+a5*h,ytemp,ak5,n,dp);
  for (i=0;i<n;i++)
    ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
  rk_check_phys_limit(ytemp,n,dp,TRUE);

  /* at t+a6*h */
  rk_derivs(x+a6*h,ytemp,ak6,n,dp);
  for (i=0;i<n;i++)
    yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
  rk_check_phys_limit(ytemp,n,dp,TRUE);

  /* error, defined as difference between fourth and fifth order method */
  for (i=0;i<n;i++)
    yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
  /* free arrays */
  free(ytemp);
  free(ak6);
  free(ak5);
  free(ak4);
  free(ak3);
  free(ak2);
}

/* 

free space, assign bailout time, and return the exit code

*/
int odeint_exit(COMP_PRECISION **y,
		COMP_PRECISION **dydx, 
		COMP_PRECISION **yscal,
		COMP_PRECISION tbailout, 
		struct der_par *dp, int exit_code)
{
#ifdef FSTRACK_DEBUG
  //fprintf(stderr,"odeint_exit: exiting at t: %11g with error code: %i\n",
  //tbailout,exit_code);
#endif
  dp->tbailout = tbailout;
  free(*dydx);free(*y);free(*yscal);
  return exit_code; 
}
