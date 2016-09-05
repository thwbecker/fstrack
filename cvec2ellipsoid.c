/*
  
  convert input strain ellipsoid in the catesian eigenvector form as output 
  by output.c for the STRAIN_EVAL_CART mode into a cload of points

  reads data from stdin

  terrible coding!!

  $Id: cvec2ellipsoid.c,v 1.2 2004/04/19 18:41:02 becker Exp $

*/
#include "fstrack.h"


#define N 1000

int main(int argc, char **argv)
{
  COMP_PRECISION a[3][3],e[3],x,y,r[3],sp[3],t,lim,dx,dz,xs[N],zs[N],
    xscale=200;
  int i,j,ninc=14,l1,l2,n;
  if(argc==2)
    sscanf(argv[1],"%i",&ninc);
  // ninc is the number of subdivisions
  if(ninc>=N){fprintf(stderr,"increase array size\n");exit(-1);}
  // get the z array
  zs[0]=0.0;zs[1]=1.0;zs[2]=-1.0;
  dz=2.0/(COMP_PRECISION)ninc;
  ninc=3;
  for(y=dz;y<1.0;y+=dz){
    zs[ninc] = y;
    zs[ninc+1] = -zs[ninc];
    ninc+=2;
  }
  // start reading in ellipses
  while(fscanf(stdin,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
	       &t,(r+FSTRACK_X),(r+FSTRACK_Y),(r+FSTRACK_Z),(e),(e+1),(e+2),
	       &a[0][FSTRACK_X],&a[0][FSTRACK_Y],&a[0][FSTRACK_Z],
	       &a[1][FSTRACK_X],&a[1][FSTRACK_Y],&a[1][FSTRACK_Z],
	       &a[2][FSTRACK_X],&a[2][FSTRACK_Y],&a[2][FSTRACK_Z])==16){
    for(i=0;i<3;i++){
      if(e[i]<0){
	fprintf(stderr,"%s: error, e_%i (%g) < 0\n",argv[0],i+1,e[i]);
	exit(-1);
      }
      for(j=0;j<3;j++)// scale eigenvectors with eigenvealues
	a[i][j] *= e[i]/xscale;
    }
    for(l2=0;l2<ninc;l2++){// z loop
      lim = 1.0 - zs[l2]*zs[l2];
      xs[0]=0.0;
      n=1;
      if(lim > 0.0){
	xs[1]=lim;xs[2]=-lim;n+=2;
	dx=lim/((COMP_PRECISION)ninc*lim+2);
	for(y=dx,i=n;y<lim;i+=2,y+=dx){
	  xs[i] = y;
	  xs[i+1] = -xs[i];
	  n+=2;
	}
      }else{
	if(zs[l2]>0){
	  for(i=0;i<3;i++){
	    sp[i] = r[i];
	    sp[i] += a[2][i];
	  }
	  printf("%11g %11g %11g\n",sp[FSTRACK_X],sp[FSTRACK_Y],sp[FSTRACK_Z]);
	  for(i=0;i<3;i++){
	    sp[i] = r[i];
	    sp[i] -= a[2][i];
	  }
	  printf("%11g %11g %11g\n",sp[FSTRACK_X],sp[FSTRACK_Y],sp[FSTRACK_Z]);
	}
      }
      for(l1=0;l1<n;l1++){// loop over x
	y = lim - xs[l1]*xs[l1];
	if(y > 0.0){
	  y = sqrt(y);
	  for(i=0;i<3;i++){
	    sp[i] = r[i];
	    sp[i] += xs[l1] * a[0][i];
	    sp[i] += y * a[1][i];
	    sp[i] += zs[l2] * a[2][i];
	  }
	  printf("%11g %11g %11g\n",sp[FSTRACK_X],sp[FSTRACK_Y],sp[FSTRACK_Z]);
	  for(i=0;i<3;i++){
	    sp[i] = r[i];
	    sp[i] += xs[l1] * a[0][i];
	    sp[i] -= y * a[1][i];
	    sp[i] += zs[l2] * a[2][i];
	  }
	  printf("%11g %11g %11g\n",sp[FSTRACK_X],sp[FSTRACK_Y],sp[FSTRACK_Z]);
	}
	x = lim - xs[l1]*xs[l1];
	if(x > 0.0){
	  x = sqrt(x);
	  for(i=0;i<3;i++){
	    sp[i] = r[i];
	    sp[i] += x * a[0][i];
	    sp[i] += xs[l1] * a[1][i];
	    sp[i] += zs[l2] * a[2][i];
	  }
	  printf("%11g %11g %11g\n",sp[FSTRACK_X],sp[FSTRACK_Y],sp[FSTRACK_Z]);
	  for(i=0;i<3;i++){
	    sp[i] = r[i];
	    sp[i] -= x * a[0][i];
	    sp[i] += xs[l1] * a[1][i];
	    sp[i] += zs[l2] * a[2][i];
	  }
	  printf("%11g %11g %11g\n",sp[FSTRACK_X],sp[FSTRACK_Y],sp[FSTRACK_Z]);
	}
      }
    }
    printf("\n");
  }

  return 0;
}
