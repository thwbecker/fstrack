#include "fstrack.h"
/* 

   program to test a cartesian version of the particle tracker/strain
   integration routine as built into fstrack

   $Id: cart_tester.c,v 1.2 2004/04/19 18:40:51 becker Exp $

*/
#define ASPECT_RATIO 2.0

void main(int argc,char **argv)
{
  int i,j,k,n[3],nvtimes=1,izero=0,iunity=1,nxny,nxnynz,ngrid,index,
    mode,nsteps=100;
  COMP_PRECISION **v,dx[3],x[12],*zlevels,ddummy,ti,dt,eps=1e-10,
    lyap[3],time,tp,xstart,ystart,tfinal,xmax[3],xmin[3],
    T,W,S,xloc[3],a,b,c,d,vel,alpha,veccyl[3],
    veccart[3],xcyl[3];
  FILE *out,*out2,*out3,*out4;
  
  // general defaults
  xstart=ystart=0.5;
  mode=0;
  nsteps=10;
  tfinal=10.0;
  //
  // mode 0: velocities are specified for pure/simple shear combinations
  //
  // see mckenzie and jacqkson (1982)
  T=S=0.0;W=1.0;// simple shear
  //T=W=0.0;S=1.0;// pure shear
  //T=1.0;W=S=0.0;// isotropic contraction
  //T=1.0;S=1.0;W=0.0;//uniaxial contraction
  //
  //
  // modes 1 through 5 are cornerflow solutions corresponding to 
  // Mc Kenzie (GJRAS, 1979) cases,  a - d in table 1
  //
  vel=1.0;alpha=45;
  for(i=1;i<argc;i++){
    if(strcmp(argv[i],"-h")==0 || strcmp(argv[i],"-?")==0){// help
      help(argv);
    }
    else if(strcmp(argv[i],"-alpha")==0){
      advance_argument(&i,argc,argv);// assign alpha angle for cornerflow
      sscanf(argv[i],"%lf",&alpha);
    }
    else if(strcmp(argv[i],"-vel")==0){// assign velocity for cornerflow
      advance_argument(&i,argc,argv);
      sscanf(argv[i],"%lf",&vel);
    }
    else if(strcmp(argv[i],"-mode")==0){// operation mode
      advance_argument(&i,argc,argv);
      sscanf(argv[i],"%i",&mode);
    }
    else if(strcmp(argv[i],"-vel")==0){// assign velocity for cornerflow
      advance_argument(&i,argc,argv);
      sscanf(argv[i],"%lf",&vel);
    }
    else if(strcmp(argv[i],"-tfinal")==0){// assign final time
      advance_argument(&i,argc,argv);
      sscanf(argv[i],"%lf",&tfinal);
    }
    else if(strcmp(argv[i],"-nsteps")==0){// number output steps
      advance_argument(&i,argc,argv);
      sscanf(argv[i],"%i",&nsteps);
    }
    else if(strcmp(argv[i],"-xstart")==0){// tracer start location x
      advance_argument(&i,argc,argv);
      sscanf(argv[i],"%lf",&xstart);
    }
    else if(strcmp(argv[i],"-ystart")==0){// tracer start location y
      advance_argument(&i,argc,argv);
      sscanf(argv[i],"%lf",&ystart);
    }
    else if(strcmp(argv[i],"-S")==0){//  pure shear type velocity
      advance_argument(&i,argc,argv);
      sscanf(argv[i],"%lf",&S);
    }
    else if(strcmp(argv[i],"-W")==0){// simple shear type vorticity
      advance_argument(&i,argc,argv);
      sscanf(argv[i],"%lf",&W);
    }
     
    else{
      fprintf(stderr,"%s: option ""%s"" not understood\n\n",argv[0],argv[i]);
      help(argv);
    }
  }
  fprintf(stderr,"%s: tracer location x/y(t=0): %g/%g t_final: %g nsteps: %i mode: %i\n",
	  argv[0],xstart,ystart,tfinal,nsteps,mode);
  if(mode == 0)
    fprintf(stderr,"%s: pure/simple shear type velocities: S: %g W: %g T: %g\n",
	    argv[0],S,W,T);
  else{
    fprintf(stderr,"%s: corner flow, type %i: alpha: %g vel: %g\n",
	    argv[0],mode,alpha,vel);
    alpha = DEG2RAD(alpha);  
    calc_cornerflow_constants(mode-1, vel, alpha, &a, &b, &c, &d,argv);
  }
  // set up velocity grid
  ngrid = 50;
  // horizontal grid goes from 0 ... 1-dx, periodic
  for(i=0;i<2;i++){
    if(i==0){
      n[i] = ngrid*(int)ASPECT_RATIO;
      dx[i] = ASPECT_RATIO/(COMP_PRECISION)(n[i]);
    }else{
      n[i] = ngrid;
      dx[i] = 1.0/(COMP_PRECISION)(n[i]);
    }
  }
  // vertical grid goes from 0 ... 1 
  n[FSTRACK_Z] = 4;
  zlevels=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*n[FSTRACK_Z]);
  dx[FSTRACK_Z] = 1.0/((COMP_PRECISION)n[FSTRACK_Z]-1.0);
  for(i=0;i<n[FSTRACK_Z];i++)
    zlevels[i] = (COMP_PRECISION)i * dx[FSTRACK_Z];
  nxny = n[FSTRACK_X] * n[FSTRACK_Y]; 
  nxnynz = nxny * n[FSTRACK_Z];
  v=(COMP_PRECISION **)malloc(sizeof(COMP_PRECISION *)*3);
  for(i=0;i<3;i++){
    v[i] = (COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*nxnynz);
    if(! *(v+i)){fprintf(stderr,"%s: memory error\n",argv[0]);exit(-1);
    }
  }
  // bounds of box
  for(i=0;i<3;i++){
    xmin[i]=0.0;
    xmax[i]=dx[i]*(n[i]-1);
  }
  //
  // initialize velocities 
  //
  for(i=0;i<n[FSTRACK_Z];i++){
    x[FSTRACK_Z] = xmin[FSTRACK_Z] + (COMP_PRECISION)i * dx[FSTRACK_Z];
    for(j=0;j<n[FSTRACK_Y];j++){
      x[FSTRACK_Y] = xmin[FSTRACK_Y] + (COMP_PRECISION)j * dx[FSTRACK_Y];
      for(k=0;k<n[FSTRACK_X];k++){
	x[FSTRACK_X] = xmin[FSTRACK_X] + (COMP_PRECISION)k * dx[FSTRACK_X];
	index = i * nxny + j*n[FSTRACK_X] + k;
	if(mode==0){// pure/simple shear solution
	  xloc[FSTRACK_X] = x[FSTRACK_X] - (ASPECT_RATIO/2.0);
	  xloc[FSTRACK_Y] = x[FSTRACK_Y] - 0.5;
	  xloc[FSTRACK_Z] = x[FSTRACK_Z] - 0.5;
	  *(v[FSTRACK_X]+index) =  (S-T) * xloc[FSTRACK_X] + W * xloc[FSTRACK_Y];
	  *(v[FSTRACK_Y]+index) = -(S+T) * xloc[FSTRACK_Y];
	  *(v[FSTRACK_Z]+index) =                                2.0*T*xloc[FSTRACK_Z];
	}else{// cornerflow solution
	  cart_to_shift_cyl(x,xcyl,mode);
	  cylvel(xcyl,veccyl,a,b,c,d);/* obtain cylindrical velocity 
					 components */
	  cylvec2cartvec(xcyl, veccyl, veccart);// convert to cartesian vel
	  if(mode<=2){
	    *(v[FSTRACK_X]+index) =  veccart[FSTRACK_X];
	    *(v[FSTRACK_Y]+index) =  veccart[FSTRACK_Y];
	    *(v[FSTRACK_Z]+index) =  0.0;
	  }else{
	    *(v[FSTRACK_X]+index) =   veccart[FSTRACK_Y];
	    *(v[FSTRACK_Y]+index) =  -veccart[FSTRACK_X];
	    *(v[FSTRACK_Z]+index) =  0.0;
	  }
	}
      }
    }
  }
  fprintf(stderr,"%s: velocities at z=0 in vx.bin, vy.bin: nx: %i ny: %i dx:%g/%g max:%g/%g\n",
	  argv[0],n[FSTRACK_X],n[FSTRACK_Y],dx[FSTRACK_X],dx[FSTRACK_Y],xmax[FSTRACK_X],xmax[FSTRACK_Y]);
  out=fopen("vx.bin","w");fwrite(v[FSTRACK_X],sizeof(double),nxny,out);fclose(out);
  out=fopen("vy.bin","w");fwrite(v[FSTRACK_Y],sizeof(double),nxny,out);fclose(out);
  fprintf(stderr,"%s: vel.hdr holds: nx ny dx dy xmin ymin xmax ymax\n",
	  argv[0]);
  out=fopen("vel.hdr","w");
  fprintf(out,"%i %i %g %g  %g %g %g %g\n",
	  n[FSTRACK_X],n[FSTRACK_Y],dx[FSTRACK_X],dx[FSTRACK_Y],xmin[FSTRACK_X],xmin[FSTRACK_Y],xmax[FSTRACK_X],xmax[FSTRACK_Y]);fclose(out);

  // iniitialize tracer
  x[FSTRACK_X]=xstart;
  x[FSTRACK_Y]=ystart;
  x[FSTRACK_Z]=0.0;
  for(i=3;i<12;i++)// F maxtrix (x[3...11]) starts as unity
    x[i]=0.0;
  x[3]=x[7]=x[11]=1.0;

  // set up time stepping
  time = 0.0;
  tp = (tfinal-time)/(COMP_PRECISION)(nsteps-1);

  // output files
  out=fopen("tracer.xyzt","w");// location and time
  out2=fopen("tracer.tF","w");// time and F matrix
  out3=fopen("tracer.tLe","w");// time and eigensystem of L=(F.F^T)^0.5
  out4=fopen("tracer.tftheta","w");// time theta[deg](cyl coord) log_10(a/b)
  //
  fprintf(stderr,"%s: writing to tracer.xyzt, tracer.tF, tracer.tftheta and tracer.tLe\n",
	  argv[0]);
  //
  strain_eigen_out(out3, x, time);
  fprintf(out,"%20.10e %20.10e %20.10e %20.10e\n",x[FSTRACK_X],x[FSTRACK_Y],x[FSTRACK_Z],time);
  fprintf(out2,"%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n",
	  time,x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11]);// time F_ij
  cyl_strain_out(out4,x,time,mode);
  while(time < tfinal){
    ti = time;
    time += tp;
    dt = (time - ti)/10.0;
    ellsphere_cart(v[FSTRACK_X],v[FSTRACK_Y],v[FSTRACK_Z],n,(n+1),(n+2),dx,(dx+1),zlevels,&ddummy,
		   &nvtimes,x,&ti,&time,&dt,&iunity,&eps,&izero,&ddummy,lyap,
		   &ddummy,&ddummy,&izero,&ddummy);
    fprintf(out,"%20.10e %20.10e %20.10e %20.10e\n",x[FSTRACK_X],x[FSTRACK_Y],x[FSTRACK_Z],time);
    strain_eigen_out(out3, x, time);
    fprintf(out2,"%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n",
	    time,x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11]);// time F_ij
    cyl_strain_out(out4,x,time,mode);
  }
  fclose(out);fclose(out2);fclose(out3);fclose(out4);
  for(i=0;i<3;i++)
    free((v+i));
}

// calculate theta(cyl system)[deg] log_10(a/b) output 
void cyl_strain_out(FILE *out,COMP_PRECISION *x,COMP_PRECISION time,int mode)
{
  COMP_PRECISION xcyl[3],l2[9],evec[9],eval[3];
  cart_to_shift_cyl(x,xcyl,mode);  
  // get L^2 from F
  calc_l2_left_stretch((x+3),l2);
  calc_eigensystem_sym(l2,eval,evec,FALSE);
  fprintf(out,"%g %g %g\n",time,RAD2DEG(xcyl[THETA]),
	  log10(eval[2]/eval[0])*0.5);
}
//
// go from cartesian to shifted and/or flipped cartesian, 
// then cylindrical coordinate system
//
void cart_to_shift_cyl(COMP_PRECISION *x,COMP_PRECISION *xcyl,int mode)
{
  COMP_PRECISION xloc[3];
  if(mode<=2){
    // shift origin
    xloc[FSTRACK_X] = x[FSTRACK_X] - ASPECT_RATIO/2.0;
    xloc[FSTRACK_Y] = x[FSTRACK_Y] - 1.0;
    xloc[FSTRACK_Z] = x[FSTRACK_Z] ;
  }else{
    // shift and flip
    xloc[FSTRACK_X] = -(x[FSTRACK_Y] - 1.0);
    xloc[FSTRACK_Y] =  (x[FSTRACK_X] - ASPECT_RATIO/2.0);
    xloc[FSTRACK_Z] = x[FSTRACK_Z];

  }
  cart2cyl(xloc,xcyl);// go from cartesian to cylindrical system
}
	  
void strain_eigen_out(FILE *out, COMP_PRECISION *x, COMP_PRECISION time)
{
  COMP_PRECISION eval[3],evec[9],l2[9];
  // calculate the L^2 left stretch matrix
  calc_l2_left_stretch((x+3),l2);
  calc_eigensystem_sym(l2,eval,evec,TRUE);
  fprintf(out,"%5g %11g %11g %11g %15.8e %15.8e %15.8e %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n",
	  time,x[FSTRACK_X],x[FSTRACK_Y],x[FSTRACK_Z],
	  sqrt(eval[2]),sqrt(eval[1]),
	  sqrt(eval[0]),
	  evec[6+FSTRACK_X],evec[6+FSTRACK_Y],evec[6+FSTRACK_Z],
	  evec[3+FSTRACK_X],evec[3+FSTRACK_Y],evec[3+FSTRACK_Z],
	  evec[  FSTRACK_X],evec[  FSTRACK_Y],evec[  FSTRACK_Z]);
}
void advance_argument(int *i,int argc, char **argv)
{
  if(argc <= *i + 1){// no arguments left
    fprintf(stderr,"%s: input parameters: error: option ""%s"" needs a value\n",
	    argv[0],argv[*i]);
    exit(-1);
  }
  *i += 1;
}
void help(char **argv)
{
  fprintf(stderr,"%s\ncalculates tracer path and strains for analytical solutions\n",
	  argv[0]);
  fprintf(stderr,"\noptions:\nseveral\n");
  exit(-1);
  
}
