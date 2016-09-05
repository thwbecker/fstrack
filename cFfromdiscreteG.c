#include "fstrack.h"
/* 
   for testing purposes, read cartesian velocity gradient matrix from
   stdin and calculate strain by forward Euler method.

   d_t F = G . F  

   F(t+dt) = F(t) + (G . F(t)) dt

*/

int main(int argc, char **argv)
{
  /* initialize the deformation matrix as the identity matrix */
  COMP_PRECISION g[3][3],x[3],time, oldtime,dt,
    f[3][3]={{1,0,0},{0,1,0},{0,0,1}},df[3][3],l2[3][3];
  int n,i,j;
  oldtime = 0.0;
  n=0;
  /* 
     read in VGM in cartesian coordinates 
     
     G[i][j] = d_j v_i , i.e. we read in

     x y z time d_x(v_x) d_y(v_x) d_z(v_x) d_x(v_y) ....

  */
  while(fscanf(stdin,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
	       (x+FSTRACK_X),(x+FSTRACK_Y),(x+FSTRACK_Z),&time,
	       &g[FSTRACK_X][FSTRACK_X],&g[FSTRACK_X][FSTRACK_Y],&g[FSTRACK_X][FSTRACK_Z],
	       &g[FSTRACK_Y][FSTRACK_X],&g[FSTRACK_Y][FSTRACK_Y],&g[FSTRACK_Y][FSTRACK_Z],
	       &g[FSTRACK_Z][FSTRACK_X],&g[FSTRACK_Z][FSTRACK_Y],&g[FSTRACK_Z][FSTRACK_Z]) == 13){
    n++;
    if(n > 1){
      dt = time - oldtime;
      /* 
	 dF = G . F, do matrix multiplication
      */
      calc_a_times_b_3x3mat(g,f,df); 
      /* 
	 increment the deformation matrix

	 F = F + dF*dt 
	 
      */
      for(i=0;i<3;i++)
	for(j=0;j<3;j++)
	  f[i][j] += df[i][j] * dt;
      /* calculate L^2 left-stretch of F 
	 L2 = F . F^T 
      */
      calc_a_dot_at_3x3(f,l2);
      /* output of location and upper right half of L^2 */
      printf("%g %g %g %g %g %g %g %g %g\n",
	     x[FSTRACK_X],x[FSTRACK_Y],x[FSTRACK_Z],
	     l2[FSTRACK_X][FSTRACK_X],l2[FSTRACK_X][FSTRACK_Y],l2[FSTRACK_X][FSTRACK_Z],
	     l2[FSTRACK_Y][FSTRACK_Y],l2[FSTRACK_Y][FSTRACK_Z],l2[FSTRACK_Z][FSTRACK_Z]);
    }
    oldtime = time;
  }
  return 0;
}
