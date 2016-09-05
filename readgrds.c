#include "fstrack_flow.h"
/*

  read in velocities from a set of GMT grd files named 
  
  vr.i.grd, vt.i.grd, and vp.i.grd
  
  or 1/vr.i.grd 1/... through n/vr.i.grd ..
  
  where i runs from 1 to N, and N is the number of lines in the file dfilename
  which has the depth of each layer in positive numbers in units of km

  and the 1/ ... n/ directory mode is chosen if velocity fields at different
  times are specified

  the following variables refer to the model structure, eg. vr means MDP->vr

  on return, vr, vt, and vp will hold the velocities in r, theta, and phi direction
  on n[R] layers with n[PHI] and n[THETA] points in longitudinal and 
  latitudinal direction, resp

  n[5]: n[R], n[PHI], n[THETA], n[TPPROD] which is n[PHI]*n[THETA]
  and n[NRNTNP] which is n[R]*n[THETA]*n[R]

  r will hold the radial coordinate of each layer in ascending order

  all theta/phi coordinates of the grids will be from 0+dtheta/2 .. Pi-dtheta/2 
   and 0 .. 2Pi-dphi
  
  dtheta = Pi/n[THETA]
  dphi =  2Pi/n[PHI]
  
  velscale is the scaling velocity, output is in cm/yr NOT deg


  $Id: readgrds.c,v 1.20 2010/12/31 18:59:17 becker Exp $

  input_mode determines if we read GMT grd files or our own double prec binary 
  format

*/

void read_vel_grids(struct mod *model, my_boolean verbose,
		    my_boolean zero_boundary_vr)
{
  FILE *in,*out;
  //int dummy[4]={0,0,0,0};	/* GMT < 4.5.1 */
  GMT_LONG dummy[4]={0,0,0,0};	/* GMT >= 4.5.1 */
  int i,j,k,l,level,os,os1,ivt,*index,argc=1;
  char *argv[1];
  my_boolean init=FALSE,wraparound=FALSE,pixelreg=FALSE,frozen,
    weighted=FALSE;
  char sname[STRLEN],suffix[10],prefix[20],vsfile_loc[STRLEN];
  float *fgrd,fdummy,*weights;
  double *dgrd;
  COMP_PRECISION minphi,mintheta,omaxphi,maxtheta,std[4],rms[4],
    mean[4],theta,tmp=0.0,x[12],dx[12],vgm[9],vp[3];
  struct GRD_HEADER header[1];
  static my_boolean vel_init = FALSE;
  in = out = NULL;
  fgrd=NULL;dgrd=NULL;
  minphi = HC_FLT_MAX;omaxphi=HC_FLT_MIN;mintheta=HC_FLT_MAX;
  maxtheta = HC_FLT_MIN;weights=NULL;
  if(model->amode == ONLY_VEL_STATS)// only if we are interested in the velocity
    weighted = TRUE;                // statistics will we weigh the mean by area
  if(!vel_init){
    //
    // read time intervals for velocities from file
    read_time_intervals(model);
    //
    // depth layers from files
    read_depth_levels(model,&index);
    /*
      
      read the velocities in binary format, either GMT grd or double bin
      
    */
    if(MDP->n[FSTRACK_R] <= 1)
      PEE("read_vel_grids: error: should have more than one layer, check depth file");
    if(model->read_gmt){
      fprintf(stderr,"read_vel_grids: reading grd files\n");
      GMT_program="read_vel_grids";
      strcpy(suffix,"grd");
      argv[0] = GMT_program;
      GMT_begin (argc, argv);
      GMT_grd_init (header, argc, argv, FALSE);
    }else{
      fprintf(stderr,"read_vel_grids: reading bin files\n");
      strcpy(suffix,"bin");
    }
    if(model->amode == ONLY_VEL_STATS){
      sprintf(vsfile_loc,"%s.%s",VSFILE,model->ostring);
      fprintf(stderr,"read_vel_grids: writing z rms_vr rms_vt rms_vp rms_vh to %s\n",
	      vsfile_loc);
      out=myopen(vsfile_loc,"w","read_vel_grids");
    }
    for(ivt=0;ivt < MDP->nvtimes;ivt++){
      if(MDP->history)
	fprintf(stderr,"read_vel_grids: reading velocities for time [%12g, %12g] from %3i/\n",
		MDP->vtimes[ivt*3],MDP->vtimes[ivt*3+2],ivt+1);
      for(i=0;i<MDP->n[FSTRACK_R];i++){
	// determine number of grd file based on resorted arrays
	level = index[i]+1;// level numbers should go from 1 .. N 
	for(j=0;j<3;j++){
	  if(MDP->history)
	    sprintf(prefix,"%i/",ivt+1);
	  else
	    sprintf(prefix,"./");
	  // filenames
	  if(j==0)
	    sprintf(sname,"%svr.%i.%s",prefix,level,suffix);
	  else if(j==1)
	    sprintf(sname,"%svt.%i.%s",prefix,level,suffix);
	  else
	    sprintf(sname,"%svp.%i.%s",prefix,level,suffix);
	  if(model->read_gmt){
#ifdef USE_GMT3		/* old */
	    if(GMT_read_grd_info (sname,header) == -1){
	      fprintf(stderr,"read_vel_grids: error opening GMT grd file %s\n",sname);
	      exit(-1);
	    }
#else  /* new */
	    sprintf(header->name,"%s",sname);
	    if (GMT_read_grd_info (sname,header)) {
	      fprintf (stderr, "read_vel_grids: error opening file %s\n", sname);
	      exit(-1);
	    }
#endif
	  }else{
	    in=myopen(sname,"r","read_vel_grids");
	    // read header type of information
	    header->node_offset=FALSE;
	    fread(&header->x_min, sizeof(double), 1, in);
	    fread(&header->x_max, sizeof(double), 1, in);
	    fread(&header->y_min, sizeof(double), 1, in);
	    fread(&header->y_max, sizeof(double), 1, in);
	    fread(&header->x_inc, sizeof(double), 1, in);
	    fread(&header->y_inc, sizeof(double), 1, in);
	    fread(&header->nx, sizeof(int), 1, in);
	    fread(&header->ny, sizeof(int), 1, in);
	  }
	  if(!init){
	    /* 
	       obtain grid dimensions and check if they are the way we
	       like it, ie.  lon lat such that
	       
	       0 <= phi <= 2pi-dphi and 0+dtheta/2<=theta<=Pi-dtheta/2 */
	    pixelreg=(header->node_offset?TRUE:FALSE);
	    minphi=  LONGITUDE2PHI(header->x_min+(pixelreg?header->x_inc/2.0:0.0));	
	    omaxphi= LONGITUDE2PHI(header->x_max-(pixelreg?header->x_inc/2.0:0.0));
	    maxtheta=LATITUDE2THETA(header->y_min+(pixelreg?header->y_inc/2.0:0.0));
	    mintheta=LATITUDE2THETA(header->y_max-(pixelreg?header->y_inc/2.0:0.0));
	    MDP->dphi=  DEG2RAD( header->x_inc);
	    MDP->dtheta=DEG2RAD( header->y_inc);
	    if(DIFFERENT(minphi,0.0) || DIFFERENT(mintheta,MDP->dtheta*0.5) || 
	       DIFFERENT(maxtheta,PI-MDP->dtheta*0.5) || 
	       (DIFFERENT(omaxphi,TWOPI) && DIFFERENT(omaxphi,TWOPI-MDP->dphi))){
	      fprintf(stderr,"read_vel_grids: expecting 0/360/-90+dy*0.5/90-dy*0.5 range, problem with %s\n",
		      sname);
	      fprintf(stderr,"read_vel_grids: in radians: %g/%g/%g/%g\n",
		      mintheta,maxtheta,minphi,omaxphi);
	      fprintf(stderr,"read_vel_grids: xy extrame: %g %g %g %g\n",
		      header->x_min,header->x_max,header->y_min,header->y_max);
	      exit(-1);
	    }
	    //
	    // check if we should throw away double entries at 0 and 360
	    if(!DIFFERENT(omaxphi,TWOPI)){
	      MDP->n[FSTRACK_PHI] = header->nx - 1;
	      wraparound = TRUE;
	    }else{
	      MDP->n[FSTRACK_PHI] = header->nx;
	      wraparound = FALSE;
	    }
	    MDP->n[FSTRACK_THETA] = header->ny;
	    if(DIFFERENT(MDP->dtheta,PI/
			 ((COMP_PRECISION)(MDP->n[FSTRACK_THETA])))||
	       DIFFERENT(MDP->dphi,TWOPI/
			 ((COMP_PRECISION)(MDP->n[FSTRACK_PHI])))){
	      fprintf(stderr,"read_vel_grids: spacing error: ndx/dx phi: %g/%g theta: %g/%g\n",
		      TWOPI/MDP->n[FSTRACK_PHI],MDP->dphi,
		      PI/MDP->n[FSTRACK_THETA],MDP->dtheta);
	      exit(-1);
	    }
	    //
	    // set auxiliary grid dimensions
	    //
	    MDP->n[TPPROD] = MDP->n[FSTRACK_THETA]  * MDP->n[FSTRACK_PHI];// ny * nx
	    MDP->n[NRNTNP] = MDP->n[TPPROD] * MDP->n[FSTRACK_R];  // ny * nx * nr
	    os = MDP->n[NRNTNP] * MDP->nvtimes;//              ny * nx * nr *nt
	    //
	    // allocate space
	    my_svecalloc(&MDP->vr,os,"readgrds: vr");
	    my_svecalloc(&MDP->vt,os,"readgrds: vt");
	    my_svecalloc(&MDP->vp,os,"readgrds: vp");
	    if(model->read_gmt){
	      // this has to be of the original GRD file size
	      // NOT the new grid dimensions
	      fgrd = (float  *)malloc(sizeof(float)  * header->nx * header->ny);
	    }else{
	      dgrd = (double *)malloc(sizeof(double) * header->nx * header->ny);
	    }
	    if((model->read_gmt && !fgrd) ||(!model->read_gmt && !dgrd))
	      MEMERROR("read_vel_grids: velocity fields:");
	    if(weighted){
	      //
	      // need to construct 2-D array with area weights
	      //
	      my_svecalloc(&weights,MDP->n[TPPROD],"readgrds");
	      for(theta=mintheta,
		    k=0;k < MDP->n[FSTRACK_THETA];k++,theta += MDP->dtheta){
		tmp = sin(theta);
		for(l=0;l < MDP->n[FSTRACK_PHI];l++)
		  weights[k*MDP->n[FSTRACK_PHI]+l] = tmp;
	      }
	    }
	    fprintf(stderr,"read_vel_grids: x: %g/%g/%g nx: %i y: %g/%g/%g ny: %i wrap: %i v_c: %g\n",
		    PHI2LONGITUDE(minphi),PHI2LONGITUDE(omaxphi),
		    RAD2DEG(MDP->dphi),MDP->n[FSTRACK_PHI],
		    THETA2LATITUDE(maxtheta),THETA2LATITUDE(mintheta),
		    RAD2DEG(MDP->dtheta),MDP->n[FSTRACK_THETA],wraparound,
		    MDP->velscale);
	    init=TRUE;
	  }else{
	    if(DIFFERENT(minphi,LONGITUDE2PHI(header->x_min+(pixelreg?header->x_inc/2.0:0.0)))||
	       DIFFERENT(omaxphi,LONGITUDE2PHI(header->x_max-(pixelreg?header->x_inc/2.0:0.0)))||
	       DIFFERENT(maxtheta,LATITUDE2THETA(header->y_min+(pixelreg?header->y_inc/2.0:0.0)))||
	       DIFFERENT(mintheta,LATITUDE2THETA(header->y_max-(pixelreg?header->y_inc/2.0:0.0)))||
	       DIFFERENT(MDP->dphi,DEG2RAD(header->x_inc))||
	       DIFFERENT(MDP->dtheta,DEG2RAD( header->y_inc))){
	      fprintf(stderr,"read_vel_grids: grd files have different size, grd: %s\n",
		      sname);
	      exit(-1);
	    }
	  }
	  if(model->read_gmt)
	    // read the netcdf GRD file
#ifndef USE_GMT3		/* new */
	    GMT_read_grd (sname,header,fgrd, 0.0, 0.0, 0.0, 0.0, dummy, 0);
#else
	    GMT_cdf_read_grd (sname,header,fgrd, 0.0, 0.0, 0.0, 0.0, dummy, 0);
#endif
	  else{
	    fread(dgrd,sizeof(double),header->nx*header->ny,in);
	    fclose(in);
	  }
	  //
	  // leave velocities in cm/yr
	  //
	  // AND: leave those pointer calculations here, since we
	  // do not initially have the size of the arrays
	  //
	  os1  = MDP->n[NRNTNP] * ivt;
	  os1 += MDP->n[TPPROD] * i;
	  // these should theoretically be == zero
	  if(j == FSTRACK_R){
	    //
	    // vr
	    //
	    if((zero_boundary_vr)&&(1.0-MDP->rlevels[i] < EPS_PREC)){
	      fprintf(stderr,"read_vel_grids: WARNING: assuming level %3i is at surface and setting vr to zero\n",
		      level);
	      resort_and_check((MDP->vr+os1),fgrd,dgrd,MDP->n[FSTRACK_PHI],
			       MDP->n[FSTRACK_THETA],wraparound,1.0/MDP->velscale,
			       model->read_gmt,TRUE,0.0);
	    }else if((zero_boundary_vr)&&(MDP->rlevels[i]<CMB_RL)){
	      fprintf(stderr,"read_vel_grids: WARNING: assuming level %3i is at CMB     and setting vr to zero\n",
		      level);
	      resort_and_check((MDP->vr+os1),fgrd,dgrd,MDP->n[FSTRACK_PHI],
			       MDP->n[FSTRACK_THETA],wraparound,1.0/MDP->velscale,
			       model->read_gmt,TRUE,0.0);
	    }else
	      resort_and_check((MDP->vr+os1),fgrd,dgrd,MDP->n[FSTRACK_PHI],
			       MDP->n[FSTRACK_THETA],wraparound,1.0/MDP->velscale,
			       model->read_gmt,FALSE,fdummy);
	    calc_mean_and_stddev((MDP->vr+os1),&fdummy,MDP->n[TPPROD],
				 (mean+j),(std+j),(rms+j),FALSE,weighted,weights);
	  }else if(j==FSTRACK_THETA){
	    //
	    // vtheta
	    //
	    resort_and_check((MDP->vt+os1),fgrd,dgrd,MDP->n[FSTRACK_PHI],
			     MDP->n[FSTRACK_THETA],wraparound,1.0/MDP->velscale,
			     model->read_gmt,FALSE,fdummy);
	    calc_mean_and_stddev((MDP->vt+os1),&fdummy,MDP->n[TPPROD],
				 (mean+j),(std+j),(rms+j),FALSE,weighted,weights);
	  }else{
	    //
	    // vphi
	    //
	    if(j != FSTRACK_PHI)
	      PEE("read_vel_grds: index error");
	    resort_and_check((MDP->vp+os1),fgrd,dgrd,MDP->n[FSTRACK_PHI],MDP->n[FSTRACK_THETA],
			     wraparound,1.0/MDP->velscale,
			     model->read_gmt,FALSE,fdummy);
	    calc_mean_and_stddev((MDP->vp+os1),&fdummy,MDP->n[TPPROD],
				 (mean+j),(std+j),(rms+j),FALSE,weighted,weights);
	    //
	    // and horizontal stats, put those in the 4th element 
	    // of mean
	    //
	    calc_mean_and_stddev((MDP->vp+os1),(MDP->vt+os1),MDP->n[TPPROD],
				 (mean+3),(std+3),(rms+3),TRUE,weighted,weights);
	  }
	}
	if(verbose)
	  fprintf(stderr,"%13s: l: %2i i: %2i r: %9.2e z: %9.2e %s mean/RMS: vr: %9.2e/%9.2e vt: %9.2e/%9.2e vp: %9.2e/%9.2e\n",
		  sname,level,i,MDP->rlevels[i],ZDEPTH(MDP->rlevels[i]),
		  (weighted?"weighted":"unweighted"),mean[FSTRACK_R]*MDP->velscale,
		  rms[FSTRACK_R]*MDP->velscale,mean[FSTRACK_THETA]*MDP->velscale,
		  rms[FSTRACK_THETA]*MDP->velscale,mean[FSTRACK_PHI]*MDP->velscale,
		  rms[FSTRACK_PHI]*MDP->velscale);
	if(model->amode == ONLY_VEL_STATS)// velocity statistics output
	  fprintf(out,"%14.5e %14.5e %14.5e %14.5e %14.5e %5i %13.5f\n",
		  ZDEPTH(MDP->rlevels[i]),rms[FSTRACK_R]*MDP->velscale,
		  rms[FSTRACK_THETA]*MDP->velscale,rms[FSTRACK_PHI]*MDP->velscale,
		  rms[3]*MDP->velscale,ivt+1,
		  ((MDP->history)?(MDP->vtimes[ivt*3+1]):(0.0)));
      }
    }
    free(index);
    if(model->read_gmt)
      free(fgrd);
    else
      free(dgrd);
    if(MDP->remove_strain){
      if(model->verbose)
	fprintf(stderr,"read_vel_grids: computing max hor. divergence for srs mode\n");
      /* 
	 compute the max horizontal divergence at the second to top layer 
      */
      find_max_horizontal_divergence(model,MDP->rlevels[MDP->n[FSTRACK_R]-1]);
      if(0){
	/* 
	   test the frozen/non-frozen routine, this is for debugging only
	*/
	x[FSTRACK_R] = ND_RADIUS(100);
	for(x[FSTRACK_PHI]=LONGITUDE2PHI(0.0);x[FSTRACK_PHI]<=LONGITUDE2PHI(360.0)+1e-5;x[FSTRACK_PHI]+=DEG2RAD(1.0))
	  for(x[FSTRACK_THETA]=LATITUDE2THETA(89);x[FSTRACK_THETA]<=LATITUDE2THETA(-89)+1e-5;x[FSTRACK_THETA]+=DEG2RAD(1.0)){
	    fse_derivs_wrapper(0.0,x,dx,12,MDP,vp,vgm,TRUE,TRUE,&frozen,
			       MDP->strain_fraction_from_gamma);
	    fprintf(stdout,"%g %g %i\n",PHI2LONGITUDE(x[FSTRACK_PHI]),
		    THETA2LATITUDE(x[FSTRACK_THETA]),frozen);
	  }
      }
    }
    if(model->amode == ONLY_VEL_STATS){
      fclose(out);
      fprintf(stderr,"read_vel_grids: exiting after printing vel stats\n");
      free(weights);
      exit(0);
    }
    vel_init=TRUE;
  }else{
    fprintf(stderr,"read_vel_grids: using previously read velocity files\n");
  }
}
/*

  find the maximum horizontal divergences of the velocity field at all given
  times, each at depth xr

 */
void find_max_horizontal_divergence(struct mod *model, COMP_PRECISION xr)
{
  static my_boolean init=FALSE;
  my_boolean frozen;
  int i,j,k,nwork=12;
  COMP_PRECISION work[12],dwork[12],hdiv;
  if(!init)
    my_vecalloc(&MDP->vhdivmax,MDP->nvtimes,"fmhd");
  // initialize the deformation matrix, just in case
  unity_matrix((work+3),3);
  // location
  work[FSTRACK_R] = xr;
  for(i=0;i  < MDP->nvtimes;i++){// loop though all times
    MDP->vhdivmax[i] = 0.0;
    for(j=0;j < MDP->n[FSTRACK_PHI];j++){// loop through all longitudes
      work[FSTRACK_PHI] = (COMP_PRECISION)j * MDP->dphi;
      for(k=0;k < MDP->n[FSTRACK_THETA];k++){// loop through all latitudes
	work[FSTRACK_THETA] = (0.5 + (COMP_PRECISION)k) * MDP->dtheta;
	/* 
	   find the velocity gradient matrix in spherical coordinate
	   system
	*/
	fse_derivs_wrapper(MDP->vtimes[i*3+1],work,dwork,nwork,
			   model->dp,MDP->v_loc,MDP->vgm,
			   FALSE,FALSE,&frozen,
			   MDP->strain_fraction_from_gamma);
	/* 
	   horizontal part of the divergence 
	*/
	hdiv = fabs(MDP->vgm[TT] + MDP->vgm[PP]);
	if(hdiv > MDP->vhdivmax[i])
	  MDP->vhdivmax[i] = hdiv;
      }
    }
    if(model->verbose)
      fprintf(stderr,"find_max_horizontal_divergence: max hor. div. at time %11g is %11g (r: %g)\n",
	      MDP->vtimes[i*3+1],MDP->vhdivmax[i],
	      xr);
  }
}
/* 

 */
void resort_and_check(float *a,float *fb,double *db,
		      int m, int n,my_boolean wrap,COMP_PRECISION factor,
		      my_boolean read_gmt,my_boolean set_to_constant,
		      float constant)
{
  int i,j,nm,os1,os2,boff;
  static my_boolean warned=FALSE;
  nm=m*n;
  if(read_gmt){
    // check for NaNs
    for(i=0;i<nm;i++)
      if(!finite(fb[i])){
	fb[i]=0.0;
	if(!warned){
	  fprintf(stderr,"WARNING: at least one NaN entry in the data has been replaced with zero\n");
	  warned=TRUE;
	}
      }
  }else{
    // check for NaNs
    for(i=0;i<nm;i++)
      if(!finite(db[i])){
	db[i]=0.0;
	if(!warned){
	  fprintf(stderr,"WARNING: at least one NaN entry in the data has been replaced with zero\n");
	  warned=TRUE;
	}
      }
  }
  if(read_gmt){
    // see if we should average the 0 and 360 entries in b 
    // which might thus also be of dimension (m+1)*n, really
    if(wrap){
      boff = m+1;
      for(i=os1=os2=0;i<n;i++,os1+=m,os2+=boff){
	a[os1] = ((COMP_PRECISION)((fb[os2] + fb[os2+m])/2.0))*factor;
	for(j=1;j<m;j++)
	  a[os1+j] = ((COMP_PRECISION)fb[os2+j])*factor;
      }
    }else{
      for(i=os1=0;i<n;i++,os1+=m)
	for(j=0;j<m;j++)
	  a[os1+j] = ((COMP_PRECISION)fb[os1+j])*factor;
    }
  }else{// our own format, use doubles
    if(wrap){
      boff = m+1;
      for(i=os1=os2=0;i<n;i++,os1+=m,os2+=boff){
	a[os1] = ((COMP_PRECISION)((db[os2] + db[os2+m])/2.0))*factor;
	for(j=1;j<m;j++)
	  a[os1+j] = ((COMP_PRECISION)db[os2+j])*factor;
      }
    }else{
      for(i=os1=0;i<n;i++,os1+=m)
	for(j=0;j<m;j++)
	  a[os1+j] = ((COMP_PRECISION)db[os1+j])*factor;
    }
  }
  if(set_to_constant){
#ifdef DEBUG
    fprintf(stderr,"resort_and_check: WARNING: setting this field to constant: %g\n",
	    constant);
#endif
    for(i=os1=0;i<n;i++,os1+=m)
      for(j=0;j<m;j++)
	a[os1+j] = constant;
  }
}
/* 
     
   read in depth levels from file and assign to r vector
   
*/
void read_depth_levels(struct mod *model,int **index)
{
  FILE *in;
  int i;
  COMP_PRECISION *rnew;
  my_boolean warned = FALSE;
  in=myopen(DFILE,"r","read_depth_levels");
  MDP->n[FSTRACK_R]=0;
  MDP->rlevels=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION));
  if(!MDP->rlevels)
    MEMERROR("read_depth_levels");
  while(fscanf(in,FLT_FORMAT,(MDP->rlevels+MDP->n[FSTRACK_R]))==1){
    if(MDP->n[FSTRACK_R]>1)
      if(fabs(MDP->rlevels[MDP->n[FSTRACK_R]] -MDP->rlevels[MDP->n[FSTRACK_R]-1])<EPS_PREC)
	PEE("read_depth_levels: error: two radii are at same level");
    if(MDP->rlevels[MDP->n[FSTRACK_R]] < 0){
      /* flip sign */
      MDP->rlevels[MDP->n[FSTRACK_R]] = -MDP->rlevels[MDP->n[FSTRACK_R]];
      if(!warned){
	fprintf(stderr,"read_depth_levels: WARNING: flipping sign of depth levels in %s\n",
	       DFILE);
	warned=TRUE;
      }
    }
    MDP->rlevels[MDP->n[FSTRACK_R]] = ND_RADIUS(MDP->rlevels[MDP->n[FSTRACK_R]]);
    if((MDP->rlevels[MDP->n[FSTRACK_R]] > 1)||(MDP->rlevels[MDP->n[FSTRACK_R]] < CMB_R)){
      // check for above surface or below CMB
      fprintf(stderr,"read_depth_levels: radius %g out of range\n",MDP->rlevels[MDP->n[FSTRACK_R]]);
      exit(-1);
    }
    MDP->n[FSTRACK_R]++;
    MDP->rlevels=(COMP_PRECISION *)realloc(MDP->rlevels,sizeof(COMP_PRECISION)*
				       (MDP->n[FSTRACK_R]+1));
    if(!MDP->rlevels)
      MEMERROR("read_depth_levels");
  }
  fclose(in);
  // sort and create index
  *index=(int *)malloc(sizeof(int)*MDP->n[FSTRACK_R]);
  if(! *index)
    MEMERROR("read_depth_levels");
  indexx(MDP->n[FSTRACK_R],(MDP->rlevels-1),(*index-1));
  // reassign
  rnew=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*MDP->n[FSTRACK_R]);
  for(i=0;i<MDP->n[FSTRACK_R];i++)
    rnew[i] = MDP->rlevels[(*index)[i]];
  for(i=0;i<MDP->n[FSTRACK_R];i++)
    MDP->rlevels[i] = rnew[i];
  free(rnew);
  fprintf(stderr,"read_depth_levels: read %i levels from %s, r_min: %g r_max: %g expecting units of %g cm/yr\n",
	  MDP->n[FSTRACK_R],DFILE,MDP->rlevels[0],MDP->rlevels[MDP->n[FSTRACK_R]-1],
	  (R_E*1e5/TIMESCALE)/MDP->velscale);
}
/*
  
  read in times for time history of velocities, if needed
  those are in the format t_left^1 t_mid^1 t_right^1 t_left^2 ...
  
*/
void read_time_intervals(struct mod *model)
{
  FILE *in;
  int n,n3;
  my_vecalloc(&MDP->vtimes,3,"rti: 1");
  if(MDP->history){
    in=myopen(THFILE,"r","read_time_intervals");
    n=n3=0;
    while(fscanf(in,TWO_FLT_FORMAT,(MDP->vtimes+n3),
		 (MDP->vtimes+n3+2))==2){
      if(n > 0){
	if((MDP->vtimes[n3+2] < MDP->vtimes[n3])||
	   (MDP->vtimes[n3] < MDP->vtimes[(n-1)*3])||
	   (MDP->vtimes[n3+2] < MDP->vtimes[(n-1)*3+2])||
	   (fabs(MDP->vtimes[(n-1)*3+2] - MDP->vtimes[n3])>EPS_PREC)){
	  PE("read_time_intervals: error, expecting ascending time intervals");
	  PE("read_time_intervals: which should have smaller time first and have no gaps");
	  fprintf(stderr,"read_time_intervals: %i: %g %g %i: %g %g\n",
		  n-1, MDP->vtimes[(n-1)*3],MDP->vtimes[(n-1)*3+2],
		  n,MDP->vtimes[n3],MDP->vtimes[n3+2]);
	  exit(-1);
	}
      }
      // mid point
      MDP->vtimes[n3+1] =  
	(MDP->vtimes[n3] + MDP->vtimes[n3+2])/2.0;
      n++;
      n3+=3;
      my_vecrealloc(&MDP->vtimes,n3+3,"rti: 2");
    }
    MDP->nvtimes = n;
    if(!MDP->nvtimes){
      fprintf(stderr,"read_time_intervals: error, no times read from %s\n",
	      THFILE);
      exit(-1);
    }else{
      fprintf(stderr,"read_time_intervals: read %i time intervals from %s\n",
	      MDP->nvtimes,THFILE);
      fprintf(stderr,"read_time_intervals: t_min: %g t_max: %g\n",
	      MDP->vtimes[0],MDP->vtimes[(MDP->nvtimes-1)*3+2]);
    }
    fclose(in);
    if(-MDP->vtimes[0] < model->bstime){
      model->bstime = -MDP->vtimes[0];
      /* adjust backward advection time */
      fprintf(stderr,"read_time_intervals: WARNING: adjusting backward advection time accordingly to %g\n",
	      model->bstime);
    }
  }else{
    MDP->nvtimes = 1;
    MDP->vtimes[0] = MDP->vtimes[1] = MDP->vtimes[2] = 0.0;
  }
}
