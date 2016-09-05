/* 

structure to hold all variables that belong to derivative
calculations for Runge Kutta


*/
struct der_par{
  /* 
     general parameters for Runge Kutta integration 
     (those were global in Numerical recipes)

     kmax is the number of intermediate results stored
     dxsav is the time interval for storage
  */
  int kmax,kount,nok,nbad;
  COMP_PRECISION *xp,*yp,dxsav,tbailout,hmax;
  EXT_FTRN_PREC rkeps;		/* precision of stepper */
  int error_mode;		/* for determine_error_scale */
  my_boolean strain_rate_control; /* this uses -eps_strain/strain_rate as the max
				     time step for Runge Kutta */
  COMP_PRECISION eps_strain;


  /* 
     
  additional parameters for problem

  */
  /* 

     velocity field related 

  */
  VPREC *vr, *vt, *vp; /* r, theta, phi velocities */
  int n[5];/* n[0]: n_r n[1]: n_theta n[2]: n_phi 
	      n[3]: n_theta * n_phi   n[4]: n_theta * n_phi * n_r
	   */
  COMP_PRECISION dtheta,dphi;	/* spacing */
  COMP_PRECISION *rlevels;	/* radial levels */
  COMP_PRECISION velscale;// velocity scale
  /* timedependence of v */
  my_boolean history;/* if true, read in velocities at different
		     times as specified in times.dat */
  COMP_PRECISION *vtimes;/* [nvtimes*3]
			    time intervals (t_left t_mid t_right) 
			    at which velocities are given
			 */
  int nvtimes;/* nr of time intervals at which velocities are given, ie. 
		 vtimes is vtimes(nvtimes*3) */
  /* 
     operational quantities 
  */
  my_boolean remove_strain;	/* remove the strain, leave only the rotation? 
				   (this is typically FALSE)
				*/
  my_boolean strain_fraction_from_gamma; /* compute the strain-rate
					    fraction of the velocity
					    gradient tensor using
					    er.i.grd data with ratios
					    between dislocation and
					    diffusion creep */
#ifdef FSTRACK_USE_GGRD
  struct ggrd_gt ggrd_alpha[1];
#endif  
  /* 
     derived quantities 
  */
  COMP_PRECISION *vhdivmax;	/* max divergence at a given time */
  /* 
     spherical velocities
  */
  COMP_PRECISION v_loc[3];
  /* 
     velocity gradient matrix
  */
  COMP_PRECISION vgm[9];
  /* 

  calculation of Lyapunov exponents
  
  */
  my_boolean calc_lyapunov ;
  /* 
     renormalization limit and bailout crit and time renorm is the
     critical element of F size when renomalization takes place
  */
  COMP_PRECISION renorm,rlbailout,tmax;	
  COMP_PRECISION rlyap[3];	/* exponents */
  /* 

  DREX related
  
  */
  COMP_PRECISION drex_epsfac;	/* factor for RK errors, set to unity 
				   for equal error estimates, larger than 
				   unity for less precision, or smaller than zero
				   for -drex_epsfac/strain_rate for the max time step
				*/
  int drex_size3;		/* 
				   number of grains per direction 
				*/
  COMP_PRECISION drex_chi,drex_Mob,drex_Xol,drex_lambda; /* other parameters */
  int drex_type;	/*  
			    slip system type, see drex.h
			*/
  COMP_PRECISION drex_rss_tau[4]; /* resolved shear stresses on slip systems for free format */
  my_boolean drex_save_pole; 		/* save ODFs? */
  my_boolean drex_compute_avg_axes;	/* compute the average ODF axes? */
  int drex_module;	/* T,p dependent modules modes */
  COMP_PRECISION drex_module_p, /* pressure and temperature */
    drex_module_t;
  my_boolean drex_start_grains_oriented; /* start with a defined, rather than random distribution */
  int drex_np[2],drex_npole;	/* dimensions of ODF array */
  COMP_PRECISION drex_sstp;	/* transition pressure */
  struct drex_para *drex;

  /* 
     mode for texture development:
     0: off
     1: Kaminski and Ribe
  */
  int texture_mode;
  /* averaging: 0 hill 0.5 VHR 1: reuss*/  
  COMP_PRECISION vhr_favg;
  /* 

  grain orientation rotation (before gridding ODF functions)
  if TRUE, will rotate grains from globsal Cartesian system to 
  local E-N-U system. If false, will leave as is

  */
  my_boolean rotate_grain_coord_sys;
  /* 
     print the direction cosines for olivine [100] to file for 
     each timestep?
  */
  my_boolean print_o1xyz;
  int lpo_tracer_count;		/* count the number of times, LPOs are computed  */
  my_boolean constant_vgm_for_debugging;	/* this allows setting up 
						   a constant velocity gradient a
						   matrix and should be FALSE by 
						   default
						*/
  COMP_PRECISION cvfd_stw[3];	/* 
				   those are the S, T, and W values
				   from McKenzie & Jackson (1982) which will be 
				   used to assign a constant velocity gradient matrix 
				*/
};
