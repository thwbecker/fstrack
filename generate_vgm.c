#include "fstrack.h"
/* 


generate a synthetic VGM matrix streamline

*/
int main(int argc, char **argv)
{
  COMP_PRECISION vgm[9],e[9],x[3],v[3],eval[3],cstrainr,l[9],
    dt,dstrain,t,tmax,rot,ss,ps,co,div,strain[9],l2[9],lval[3],
    lvec[9],s,w;
  int i,simmode;
  FILE *out1,*out2;

  out1 = fopen("fse.dat","w");
  out2 = fopen("vgm.dat","w");

  x[FSTRACK_X]=x[FSTRACK_Y]=x[FSTRACK_Z] = 0.;
  v[FSTRACK_X]=v[FSTRACK_Y]=v[FSTRACK_Z] =0.;
  for(i=0;i<9;i++)		/* FSE */
    strain[i]=0.0;
  strain[XX] = strain[YY] = strain[ZZ] = 1.0;

  simmode=0;			/* plate modes

				   0: simple shear all the way 
				   1: simple shear, reversed
				   2: simple shear, change by 90 deg
				   3: pure shear
				   4: pure shear, then simple shear
				   5: simple shear, continuously rotating
				   6: pure shear, continuously rotating
				   7: compression, extension

				   2-D modes:
				   10:
				   

				*/
  if(argc > 1)
    sscanf(argv[1],"%i",&simmode);
  

  tmax = 1.;
  dt = .005;
  for(t=0;t<tmax;){
    /* plate modes */
    ss = 0.;			/* simple shear */
    ps = 0;			/* pure shear */
    co = 0;			/* contraction */
    /* 2-D s/w modes */
    s = w = 0;
    rot = 0;			/* rotation of coordinate system,
				   deg */
    switch(simmode){
    case 0:
      ss = 10.;			/* simple shear */
      generate_vgm_plate(vgm,ss,ps,co,rot);
      break;
    case 1:			/* simple shear, reverse in direction */
      if( t < 0.5)
	ss = 10;
      else
	ss = -10;
      generate_vgm_plate(vgm,ss,ps,co,rot);
      break;
    case 2:			/* simple shear, then rotate */
      ss = 10;
      if(t < 0.5)
	rot = 0;
      else
	rot = 90;
      generate_vgm_plate(vgm,ss,ps,co,rot);
      break;
    case 3:
      ps = 1;
      generate_vgm_plate(vgm,ss,ps,co,rot);
      break;
    case 4:			
      if(t < 0.5)
	ps = 1;
      else
	ss = 10;
      generate_vgm_plate(vgm,ss,ps,co,rot);
      break;
    case 5:
      ss = 10.;
      rot = t*90;
      generate_vgm_plate(vgm,ss,ps,co,rot);
      break;
    case 6:
      ps = 1.;
      rot = t*90;
      generate_vgm_plate(vgm,ss,ps,co,rot);
      break;
    case 7:
      if(t < .5)
	co = 1.;
      else
	co = -1.;
      generate_vgm_plate(vgm,ss,ps,co,rot);
      break;
      /*  */
    case 10:
      s = w = 1;
      generate_vgm_sw(vgm,s,w,rot);
      break;
    case 11:
      s = 2;
      w = 1;
      generate_vgm_sw(vgm,s,w,rot);
      break;
    case 12:
      s = 1;w= 0;
      generate_vgm_sw(vgm,s,w,rot);
      break;
    case 13:
      w = -1;
      s = -2;
      generate_vgm_sw(vgm,s,w,rot);
      break;
    case 14:
      s = 1;
      w = -s;
      generate_vgm_sw(vgm,s,w,rot);
      break;
      /*  */
    default:
      fprintf(stderr,"%s: simulation mode %i undefined\n",argv[0],
	      simmode);
      exit(-1);
      break;
    }
    /* generate a velocity gradient matrix */

    /*  */
    div = remove_trace(vgm);	/* remove divergence, if any */
    if(fabs(div)>1e-7)
      fprintf(stderr,"%s: t: %g removing divergence %g\n",
	      argv[0],t,div);
    
    calc_cd_symm_part(vgm,e);	/* compute symmetric part (strain-rate) */
    cstrainr = char_strain_rate_from_e(e,eval); /* characteristic (max) strain-rate */
    dstrain = cstrainr * dt;	/* strain step = strain_rate * dt */

    /* update the FSE strain, dF/dt = VGM */
    for(i=0;i<9;i++)
      strain[i] += vgm[i] * dt;
    /* compute sqrt(left-stretch) */
    calc_l2_left_stretch(strain,l2);calc_sqrt_sym_matrix(l2,l); 
    calc_eigensystem_sym(l,lval,lvec,TRUE);
    /* check for plausibility */
    if(lval[0] <= 0){
      fprintf(stderr,"%s: error, smallest eigenvalue <= 0\n",
	      argv[0]);
      exit(-1);
    }
    if(pow(lval[0]*lval[1]*lval[2],1./3) <= 0){
      fprintf(stderr,"%s: error, volume eigenvalue <= 0\n",
	      argv[0]);
      exit(-1);
    }
    
    t += dt;

    /* output of eigenvectors, smallest to highest */
    fprintf(out1,"%11g\t",t);
    for(i=0;i<3;i++)
      fprintf(out1,"%12.8e %7.4f %7.4f %7.4f\t",
	      lval[i],lvec[i*3+FSTRACK_X],lvec[i*3+FSTRACK_Y],lvec[i*3+FSTRACK_Z]);
    fprintf(out1,"\n");
    
    fprintf(out2,"%.2f %.2f %.2f %.2f\t%12.8e %12.8e %12.8e %12.8e\t%12.8e %.5f\n",
	    x[FSTRACK_X],x[FSTRACK_Y],v[FSTRACK_X],v[FSTRACK_Y],
	    vgm[XX],vgm[XY],vgm[YX],vgm[YY],dstrain,dt);

  }
  fclose(out1);
  fclose(out2);
}

/* 

compute a velocity gradient tensor with ss, ps, co strike slip, pure
shear, contraction components and a coordinate system that is rotated
by angle(deg) rot around the z axis

*/
void generate_vgm_plate(COMP_PRECISION *vgm,
			COMP_PRECISION ss,
			COMP_PRECISION ps,
			COMP_PRECISION co,
			COMP_PRECISION rot)
{
  COMP_PRECISION s,t,w;		/* follow mckenzie EPSL 1983 
				   y is up, x is right
				
				*/
  COMP_PRECISION l[9];
  s = t = w = 0;
  w += ss;			/* simple shear, right lateral is
				   positive ---> 
				            <---
				*/
  s += ps;			/*
				  pure shear
				     |
				     v
				  <-- --> is positive
				     ^
				     |
				*/
  t += co;			/* 
				   contraction
				     |
				     v
				   -> <-  positive
				     ^
				     |

				*/
  l[XX] = s-t;l[XY] = w;     l[ZZ]=0;
  l[YX] =   0;l[YY] = -(s+t);l[YZ]=0;
  l[ZX] =   0;l[ZY] =      0;l[ZZ]=2*t;
  /* rotate matrix */
  rotate_cart_mat(l,vgm,DEG2RAD(rot),0.0,0.0);
  zero_small_entries(vgm,9);
}

void generate_vgm_sw(COMP_PRECISION *vgm,
		     COMP_PRECISION s,
		     COMP_PRECISION w,
		     COMP_PRECISION rot)
{
  int i;
  COMP_PRECISION l[9];
  for(i=0;i<9;i++)
    l[i] = 0;
  l[XX]=0;  l[XY]=s-w;
  l[YX]=s+w;l[YY]=0;
  /* rotate matrix */
  rotate_cart_mat(l,vgm,DEG2RAD(rot),0.0,0.0);
  zero_small_entries(vgm,9);


}
