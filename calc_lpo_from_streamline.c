#include "fstrack.h"
/* 


   read in velocity gradients along a streamline and compute texture 

*/
int main(int argc, char **argv)
{
  
  struct drex_para *drex=NULL;
  COMP_PRECISION tau[4]={1,2,3,1e60}; /* slip systems for olivine */
  COMP_PRECISION tau_ens = 1.0; /* slip systems for enstatite, simple scale factor */
  COMP_PRECISION Mob = 125;	/* grain boundary mobility, ~125 */
  COMP_PRECISION Xol = 70.;	/* olivine fraction in %, 70 */
  COMP_PRECISION chi = 0.3;	/* threshold for grain boundary sliding, ~0.3*/
  COMP_PRECISION lambda = 5.0;	/* nucleation parameter, ~5 */
  int size3 = 17;		/* number of grains**(1/3) */
  int odf_np[2] = {180,90};	/* dimension of pole figure array */
  int n,i,ii,tc,poleout,gc;
  my_boolean save_gamma;	/*  */
  my_boolean o1xyz_out = 2;	/* 1: only [100] 2: all three */
  COMP_PRECISION vgm[9],*drex_dx,x[3],v[3],strain[9],xn[3],l2[9],l[9],div,dt,dstrain,t,
    dtu;
  char ofilename[300];
  int itmp;
  FILE *out;
  
  int tsubd = 50;  /* subdivide time step? */
  //int npole_save = 40;			/* save pole every */
  int npole_save = 20;			/* save pole every */


  save_gamma = TRUE;
  
  if(argc>1)			/* mobility */
    sscanf(argv[1],"%lf",&Mob);

  if(argc>2)
    sscanf(argv[2],"%lf",&chi);	/* GBS */
  
  if(argc>3)
    sscanf(argv[3],"%lf",&Xol);	/* Olivine */
  
  if(argc>4)
    sscanf(argv[4],"%lf",&lambda);	/* lambda */
  
  if(argc>5)
    sscanf(argv[5],"%i",&itmp);	/* type of output, 1: [100], 2: all three */
  o1xyz_out = (my_boolean)itmp;
  if(argc>6)
    sscanf(argv[6],"%i",&size3);	/* number of grains */
  

  fprintf(stderr,"%s: using Mob: %g chi: %g Xol: %g lambda: %g o1xyz_out: %i size3: %i\n",
	  argv[0],Mob,chi,Xol,lambda,o1xyz_out,size3);
  
  
  /* init once */
  drex_initialize(&drex,0,size3,Xol,tau,tau_ens,Mob,chi,lambda,odf_np,FALSE,save_gamma);
  my_vecalloc(&drex_dx,drex->bsize,argv[0]);
  
  /*  */
 

  /* for every streamline */
  drex_init_pathline(drex);	/* drex */
  for(i=0;i<9;i++)		/* FSE */
    strain[i]=0.0;
  strain[XX] = strain[YY] = strain[ZZ] = 1.0;

  // velocity gradient matrix is transpose of grad v
  //
  //     | d_x(vx) d_y(vx) d_z(vx) |
  //     | d_x(vy) d_y(vy) d_z(vy) |
  //     | d_x(vz) d_y(vz) d_z(vz) |


  fprintf(stderr,"%s: reading vgm info from stdin\n",argv[0]);
  n=poleout=0;
  t=0.;gc=0;

  //  
  // x (km), z(km), ux(mm/yr), uz(mm/yr), T (deg C), d ux/dx, d ux/dz, d uz/dx, d uz/dz dstrain dtime
  //
  while(fscanf(stdin,"%lf %lf %lf %lf %*f %lf %lf %lf %lf %lf %lf",
	       (x+FSTRACK_X),(x+FSTRACK_Z),(v+FSTRACK_X),(v+FSTRACK_Z),
	       (vgm+XX),(vgm+XZ),(vgm+ZX),(vgm+ZZ),
	       &dstrain,&dt) == 10){
    x[FSTRACK_Y] =0.0;v[FSTRACK_Y]=0.0;
    vgm[XY] = vgm[YX] = vgm[YY] = vgm[YZ] = vgm[ZY] = 0.0;

    /* make sure we have no divergence */
    div = remove_trace(vgm);
    if(fabs(div) >1e-5)
      fprintf(stderr,"%s: WARNING: removed divergence of %.4e\n",argv[0],div);
    /* update the FSE strain, dF/dt = VGM */
    for(i=0;i<9;i++)
      strain[i] += vgm[i] * dt;
    /* location */
    if(n){
      /* 
	 Euler scheme for testing 
	 
      */    
      for(i=0;i<3;i++)
	xn[i] += v[i] * dt;
    }else{			/* first step */
      for(i=0;i<3;i++)
	xn[i] = x[i];
    }
    /* 
       move the LPO parameters in smaller steps 
    */
    dtu = dt/(COMP_PRECISION)tsubd;
    for(tc = 0;tc < tsubd; tc++){
      /* 
	 compute the DREX derivatives from the velocity gradient matrix
	 by calling the driver for drex_deriv_ftrn
      */
      drex_deriv(vgm,drex->x,drex_dx,drex);
      /* derivatives for all other components */
      for(i=0;i < drex->bsize;i++)
	drex->x[i] += drex_dx[i] * dtu;
    }
    t += dt;
    /* 
       output 
    */
    /* left stretch matrix */
    calc_l2_left_stretch(strain,l2);calc_sqrt_sym_matrix(l2,l); 
    fprintf(stdout,"%i\t%11g\t%11g %11g\t%11g %11g\t%11g %11g %11g %11g %11g %11g\n",
	    n,t,x[FSTRACK_X],x[FSTRACK_Z],xn[FSTRACK_X],xn[FSTRACK_Z],l[XX],l[XY],l[XZ],l[YY],l[YZ],l[ZZ]);
    if(n%npole_save == 0){
      /* save pole */
      drex_compute_pdens(drex,x,t,FALSE,o1xyz_out,0,TRUE);
      poleout++;
      /* print gamma? */
      if(save_gamma){
	gc++;
	sprintf(ofilename,"gamma.%i.dat",gc);
	out = fopen(ofilename,"w");
	for(ii=0;ii<drex->size;ii++)
	  fprintf(out,"%11g\t%g\t %11g %11g\n",
		  drex->odf[ii],
		  drex->gamma_save[ii*3],
		  drex->gamma_save[ii*3+1],
		  drex->gamma_save[ii*3+2]);
	fclose(out);
	fprintf(stderr,"%s: written gamma factors to %s\n",
		argv[0],ofilename);

      }
    }
    n++;
  }
  
  fprintf(stderr,"%s: processed %i steps\n",argv[0],n);
  
}




