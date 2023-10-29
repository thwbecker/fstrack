/*


calculate the infinite strain axis at any given point following
Kaminski and Ribe (G^3, 2002GC000222, 2002)


$Id: calc_isa.c,v 1.12 2004/02/22 15:30:45 becker Exp $

*/
#include "fstrack.h"

#ifdef USE_DREX
/* 

calculate the Pi parameter and the grain orientation lag using 
the drex routines, and modelled after the drex example program

input is a position xp[3]

output:

       isac[3]: ISA axis at xp
       gol:     Grain Orientation lag, the Pi parameter
       

this routine is modeled after the pipar routine in Drex


*/

void calc_pipar(COMP_PRECISION *xp, COMP_PRECISION time,
		COMP_PRECISION *isac, COMP_PRECISION *golc,
		struct der_par *dp)
{
  COMP_PRECISION xstate[12],xstate_orig[12],dxstate[12],dt;
  COMP_PRECISION vpc[3],vpb[3],vpf[3];
  COMP_PRECISION vgmc[9],vgmb[9],vgmf[9],isab[3],
    isaf[3],char_strainc,tmp;
  COMP_PRECISION theta_isa_f,theta_isa_b,dtheta;
  static int nwork = 12;
  /* assign location */
  a_equals_b_vector3d(xstate,xp);	
  unity_matrix((xstate+3),3);	/* unnecessary, really */
  /* save center point state */
  a_equals_b_vector(xstate_orig,xstate,12);
  /* 
     compute velocity and velocity gradient matrix 
     at center point xp
  */
  fse_derivs_wrapper(time,xstate_orig,dxstate,nwork,dp,vpc,vgmc);
  /* compute dt as a small fraction of 1/|v| */
  dt = 0.01 * norm3d(vpc);
  /* compute the ISA at the center from the normalized vgm */
  char_strainc = char_strain_rate_from_vg(vgmc);
  scale_vector(vgmc,1.0/char_strainc,9);
  /* 
     compute the ISA axis and GOL parameter at the center from the
     normalized velocity gradient matrix
  */
  drex_isacalc(isac,golc,vgmc);
  /* 

     march back by dt 

  */
  rk_wrapper(xstate,time,time-dt,TRUE,FALSE,dp->rkeps,dp);
  /* compute velocity and gradient there */
  fse_derivs_wrapper(time-dt,xstate,dxstate,nwork,dp,vpb,vgmb);
  /* normalize velocity */
  normalize3d(vpb);
  /* normalize the velocity gradient matrix */
  scale_vector(vgmb,1.0/char_strain_rate_from_vg(vgmb),9);
  /* 
     compute the ISA axis and GOL parameter from the normalized 
     velocity gradient matrix
  */
  drex_isacalc(isab,golc,vgmb);
  if(fabs(vec_sum(isab,3)-3) < EPS_COMP_PREC){
    /* special case, assign velocity to ISA axis */
    a_equals_b_vector3d(isab,vpb);
  }
  /* angle between ISA and flow direction */
  theta_isa_b = acos(vec3ddotp(vpb,isab));
  /* 

     march forward by dt 

  */
  a_equals_b_vector(xstate,xstate_orig,12);
  rk_wrapper(xstate,time,time+dt,TRUE,FALSE,dp->rkeps,dp);
  /* compute velocity and gradient there */
  fse_derivs_wrapper(time+dt,xstate,dxstate,nwork,dp,vpf,vgmf);
  /* normalize velocity */
  normalize3d(vpf);
  /* normalize the velocity gradient matrix */
  scale_vector(vgmf,1.0/char_strain_rate_from_vg(vgmf),9);
  /* compute the ISA axis and GOL parameter */
  drex_isacalc(isaf,golc,vgmf);
  if(fabs(vec_sum(isaf,3)-3) < EPS_COMP_PREC){
    /* special case, assign velocity to ISA axis */
    a_equals_b_vector3d(isaf,vpf);
  }
  theta_isa_f = acos(vec3ddotp(vpf,isaf));
  /* 
     Pi parameter, compute from central difference 
     change between -dt and +dr isa axes
  */
  dtheta = fabs(theta_isa_b - theta_isa_f);
  if(dtheta > PI)
    dtheta -= PI;
  /* 
     now assign to the GOL parameter for the central location 
  */
  tmp = dtheta/(2.0 * dt * char_strainc);
  *golc = MIN(*golc, tmp);
}
#endif 

/*


old functions here

*/


// output of tracer ISA axis as a function of taumax
//#define TEST_TIME_EVOLUTION

void calc_isa(struct mod *model, int atracer)
{
  COMP_PRECISION evec[3],u[3],evec1[3],
    u1[3],char_strain,xloc1[12],xloc2[12],dummy,sloc[9],domega,
    xorig[12],dt,vg[9],vgm_er[3],vgm_ei[3],vgm_ev[9];
  BOOLEAN imaginary;
  /*
  COMP_PRECISION evecns[9],evalns[3],tmaxfrac,*vgnull,l2ns[9];
  int nnull;
  */
#ifdef TEST_TIME_EVOLUTION
  COMP_PRECISION azi;
  static int nt=0;
  nt++;
#endif
  //
  // preserve initial state
  copy_vector(model->tracer[atracer].state[0].x,xorig,12);
  // time at which to evaluate flow
  model->tracer[atracer].state[0].t = model->itime;
#ifdef TEST_TIME_EVOLUTION
  for(tmaxfrac=0.0;tmaxfrac<=1.0;tmaxfrac+=0.005){
    // calculate the infinite strain axis, evec
    calculate_isa_dir(model->tracer[atracer].state[0].t,
		      0.0,// dt
		      model->tracer[atracer].state[0].x,
		      TRUE,// calculate the rate of change
		      tmaxfrac,// fraction of tmax
		      u,evec,&char_strain,
		      model->tracer[atracer].state[0].s,
		      &model->tracer[atracer].state[0].clength,
		      model,vg,FALSE);
    // output of tmaxfrac azimuth(0<a<180) e_h bailout n_tracer lon lat
    azi = RAD2DEG(atan2(evec[PHI],-evec[THETA]));
    if(azi < 0)
      azi += 360.0;
    if(azi > 180)
      azi -= 180;
    printf("%g %g %g %g   %i %g %g\n",
	   tmaxfrac,azi,sqrt(1.0-evec[R]*evec[R]),
	   model->tracer[atracer].state[0].clength,nt,
	   PHI2LONGITUDE(xorig[PHI]),
	   THETA2LATITUDE(xorig[THETA]));
    copy_vector(xorig,model->tracer[atracer].state[0].x,12);
  }
  printf("\n\n");
#else
  calculate_isa_dir(model->tracer[atracer].state[0].t,
		    0.0,model->tracer[atracer].state[0].x,
		    TRUE,1.0,u,evec,&char_strain,
		    model->tracer[atracer].state[0].s,
		    &model->tracer[atracer].state[0].clength,
		    model,vg,FALSE);
#endif
  /*
    now figure out if the ISA exists, ie. if the rate of change
    in the largest eigenvector is small enough

    the cutoff criterion is rate of EV change < k * e_char
    since we advected by 1/e_char, the rate of change is
    clength * e_char, and we can simply check against k
  */
  if(model->tracer[atracer].state[0].clength > model->k )
    model->tracer[atracer].state[0].isa_exists = FALSE;
  else
    model->tracer[atracer].state[0].isa_exists = TRUE;
  /* 
     calculate the eigensystem of the VGM 
  */
  calc_eigensystem(vg,vgm_er,vgm_ei,vgm_ev,TRUE,&imaginary);
  if(0)
    fprintf(stderr,"r1: %12.5e i1: %12.5e r2: %12.5e i2: %12.5e r3: %12.5e i3: %12.5e  e:%1i i: %1i\n",
	    vgm_er[E1],vgm_ei[E1],vgm_er[E2],vgm_ei[E2],
	    vgm_er[E3],vgm_ei[E3],
	    model->tracer[atracer].state[0].isa_exists,
	    imaginary);
  if(imaginary)
    model->tracer[atracer].state[0].isa_exists = FALSE;
  if(model->tracer[atracer].state[0].isa_exists){
    // E1(exp(G tmax) )
    /*
      print_vector(vg,9,stderr);
      fprintf(stderr,"ev(F(tmax): %12g, %12g, %12g\n",
      evec[X],evec[Y],evec[Z]);
    */
    /*
      calculate the possible null space matrix for vg
    */
    /*
      vgnull=NULL;
      find_ABzero(vg,&vgnull,3,&nnull,1e-6);
      if(nnull){ 
      print_vector(vgnull,9,stderr);
      calc_l2_left_stretch(vgnull,l2ns);
      calc_eigensystem_sym(l2ns,evalns,evecns,TRUE);
      fprintf(stderr,"ev(NS(F)) : %12g, %12g, %12g\n",
	      evecns[E1*3+X],evecns[E1*3+Y],evecns[E1*3+Z]);
	      
	      }
	      free(vgnull);
	      fprintf(stderr,"\n");
    */
    /*
      
    now determine the grain orientation lag 
    (\Gamma) paramter. it is given by 
    
    \Gamma = 1/char_strain_rate | D_t \Omega | 
    
    where 
    
    D_t = d_t + \vec{u} . \grad  
    
    is the material derivative and 
    
    \Omega = cos^-1 (\vec{u} . \vec{e})
    
    with \vec{e} the local ISA direction and 
    \vec{u} the normalized velocity direction

    we shall determine D_t Omega numerically by central 
    differences
    */
    dt = model->gdt;
    //for(dt=10;dt>=1e-6;dt/=5.0){
    // ISA at t-dt
    copy_vector(xorig,xloc1,12);
    calculate_isa_dir(model->tracer[atracer].state[0].t,-dt,
		      xloc1,FALSE,1.0,u,evec,&dummy,sloc,
		      &dummy,model,vg,FALSE);
    // ISA at t+dt
    copy_vector(xorig,xloc2,12);
    calculate_isa_dir(model->tracer[atracer].state[0].t,dt,
		      xloc2,FALSE,1.0,u1,evec1,&dummy,sloc,
		      &dummy,model,vg,FALSE);
    // dOmega/dt
    domega = fabs((omega(u1,evec1) - omega(u,evec))/(2.0*dt));
    model->tracer[atracer].state[0].gamma = 
      domega / char_strain;
    /*
      fprintf(stdout,"gamma: %g dt: %g dx: %g domega: %g\n",
      model->tracer[atracer].state[0].gamma,
      dt,dist_on_sphere(xloc1,xloc2),
      domega);
      }
      printf("\n\n");		
    */
  }else{// ISA doesn't exist, gamma doesn't make sense
    model->tracer[atracer].state[0].gamma = 0.0;
  }
}
/*
  
routine to actually calculate the ISA directions at time 'time' and
location x[0..2], x[3..9] is strain.  if dt is non-zero, will advect
tracer by dt and determine ISA axis at this point

on output, u[3] is the velocity at the endpoint, evec[3] is the
largest eigendirection, char_strain the characteristic strain rate,
f[9] the F matrix at tmax, and l2[9] the l^2 left-stretch strain that
corresponds to F(tmax)

if calc_rate_of_change is set, will calculate
slength the change in the length of the eigenvector from 
tmax = taumax * char_strain to tmax = (taumax + 1) * char_strain

vg[9] will hold the velocity gradient matrix in C fashion on return

*/
void calculate_isa_dir(// input:
		       COMP_PRECISION time, // absolute starting time
		       COMP_PRECISION dt,// advection time
		       COMP_PRECISION *x,// location + F
		       BOOLEAN calc_rate_of_change, 
		       COMP_PRECISION tmaxfrac,/* fraction of 
						  model-wide tmax
						  to use, normally 1 
					       */
		       // output, besides x(3...9) which is also
		       // output
		       COMP_PRECISION *u,/* velocity at end, 
					    in velocity units */
		       COMP_PRECISION *evec,// evec at end
		       COMP_PRECISION *char_strain,
		       //
		       // normalized
		       COMP_PRECISION *l2,// L^2(F(tmax))
		       COMP_PRECISION *clength, /* change in the
						   largest eigenvector
						   at tmax */
		       struct mod *model,
		       COMP_PRECISION *vg, /* velocity gradient
					      matrix, will be returned
					      on output */
		       BOOLEAN verbose /* some output */
		       ) 
{
  static int nwork=12;
  COMP_PRECISION dx[12],val[3],vec[9],
    tmax,tmax1,ftmax1[9],l1[9],tloc,taumax,char_time,norm_strain;
  int i;
  if(tmaxfrac == 1.0)
    taumax = model->taumax;
  else
    taumax = model->taumax * tmaxfrac;
  if(dt != 0.0){// have to advect from time to tloc=time+dt
    tloc = time + dt;
    rk_wrapper(x,time,tloc,TRUE,FALSE,1e-12,MDP);
  }else{// no advection
    tloc = time;
  }
  //
  // for time = tloc and location xloc calculate the velocity (it's
  // dx[0..2] and the velocity gradient matrix G, it's dx[3...11], if
  // unity matrix is passed in. dx[3...1] is identical to MDP->vgm
  //
  unity_matrix((x+3),3);
  fse_derivs_wrapper(tloc,x,dx,nwork,model->dp,MDP->vp,MDP->vgm);
  /* velocities are output */
  copy_3dvector(MDP->vp,u);
  /* VGM matrix is output in C style */
  copy_vector(MDP->vgm,vg,9);
  *char_strain = char_strain_rate_from_vg(vg);
   if(*char_strain == 0.0){
    fprintf(stderr,"calculate_isa_dir: error: all strains are zero at x: (r:%g, t:%g, p:%g) t: %g\n",
	    x[R],x[THETA],x[PHI],time);
    exit(-1);
  }
  char_time = 1.0/(*char_strain);
  //
  // the "infinite" time is tau_max / char_strain rate
  // i.e., char_strain is 1 / characteristic time
  //
  tmax  = taumax * char_time;
  //fprintf(stderr,"calculate_isa_dir: strain rate: %g tmax: %g\n",
  //	  *char_strain,tmax);
  /*
   
  calculate the solution of d_t F = G . F at constant G,
  F(t) = exp(G t) = I + G.t + G^2/2! t^2 + G^3/3! t^3 + ...

  (x+3) = exp(vg tmax)

  */
  calc_exp_matrixt(vg, tmax, (x+3));
  if(verbose){
    print_vector_wl(vg,9,stderr,"VGM");
    print_vector_wl((x+3),9,stderr,"EVG");
  }
  /* 
     normalize x+3 strain 

  */
  norm_strain = max_vec((x+3),9);
  if(!norm_strain){
    fprintf(stderr,"calculate_isa_dir: error, max F(tmax) value is zero\n");
    exit(-1);
  }
  //
  // calculate the left-stretch tensor of F(tmax)
  //
  calc_l2_left_stretch((x+3),l2);
  // calculate eigenvectors of L^2 
  calc_eigensystem_sym(l2,val,vec,TRUE);
  for(i=0;i<3;i++)
    if(!finite(val[i])){
      fprintf(stderr,"calculate_isa_dir: error, eigenvalue %i not finite, will print L2 and vec\n",
	      i+1);
      print_vector(l2,9,stderr);
      print_vector(vec,3,stderr);
      exit(-1);
    }
  if(verbose)
    fprintf(stderr,"E1V: %12.5e %12.5e %12.5e  AZI: %12.5e E1: %12.5e E2: %12.5e E3: %12.5e\n\n",
	    vec[E1*3+R],vec[E1*3+PHI],vec[E1*3+THETA],
	    atan2(vec[E1*3+PHI],-vec[E1*3+THETA]),
	    val[E1],val[E2],val[E3]);
  //
  // assign to evec
  for(i=0;i<3;i++)
    evec[i] = vec[E1*3+i];
  if(calc_rate_of_change){
    // the next step for determining the rate of change in 
    // increment by one characteristic time
    tmax1 = tmax + char_time;
    // calculate F(tmax1) = exp(G tmax1 ), L^2(F(tmax1)), 
    // and eigensystem
    calc_exp_matrixt(vg,tmax1,ftmax1);
    calc_l2_left_stretch(ftmax1,l1);
    calc_eigensystem_sym(l1,val,vec,TRUE);
    if(vec3ddotp(evec,(vec+E1*3)) < 0){ /* flip direction 
					   if pointing 
					   in opposite ways */
      scale_vector(vec,-1.0,9);
    }
    // calculate the length of the difference vector between 
    // tmax and tmax1 (eigenvectors are normalized)
    *clength = vector_diff(evec,(vec+E1*3),3);
  }
}
/*
calculate

\Omega = cos^-1 (\vec{u} . \vec{e})

the input velocity u has units, has to be scaled to direction

*/
COMP_PRECISION omega(COMP_PRECISION *u, COMP_PRECISION *e)
{
  COMP_PRECISION unorm[3],ul;
#ifdef DEBUG
  if(fabs(norm3d(e)-1.0)>EPS_COMP_PREC){
    fprintf(stderr,"omega: error: eigenvector length %g not unity\n",
	    norm3d(e));
    exit(-1);
  }
#endif
  ul = norm3d(u);
  copy_3dvector(u,unorm);
  if(ul != 0.0)
    scale_vector3d(unorm,1.0/ul);
  return acos(vec3ddotp(unorm,e));
}
