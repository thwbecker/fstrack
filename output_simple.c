/*
  
simple output routines that do not require other subroutines

$Id: output_simple.c,v 1.1 2010/12/31 18:59:17 becker Exp $
  
*/


#include "fstrack.h"


//
// print vector in row format
//
void print_vector(COMP_PRECISION *t, int n, FILE *out)
{
  int i;
  for(i=0;i<n;i++)
    fprintf(out,"%11g ",t[i]);
  fprintf(out,"\n");
}
/* to be called from FORTRAN for debugging purposes */
void print_vector_stderr__(COMP_PRECISION *t, int *n)
{
  print_vector(t,*n,stderr);
}
//
// print vector with label in row format
//
void print_vector_wl(COMP_PRECISION *t, int n, FILE *out,
		     char *label)
{
  fprintf(out,"%s: ",label);
  print_vector(t,n,out);
}

/* print a symmetric 6 by 6 stiffness matrix, C sorted, in upper right
   hand side format */
void print_sym_6by6(COMP_PRECISION *sav,FILE *out)
{
  int j,k;
  for(j=0;j < 6;j++)
    for(k=j;k < 6;k++)
      fprintf(out,"%9.5e ",sav[j*6+k]);
  fprintf(out,"\n");
}
void print_6by6_nice(COMP_PRECISION *sav,FILE *out)
{
  int j,k;
  for(j=0;j < 6;j++){
    for(k=0;k < 6;k++)
      fprintf(out,"%10.4f ",sav[j*6+k]);
    fprintf(out,"\n");
  }
}
void print_6by6_hprec(COMP_PRECISION *sav,FILE *out)
{
  int j,k;
  for(j=0;j < 6;j++){
    for(k=0;k < 6;k++)
      fprintf(out,"%23.14e ",sav[j*6+k]);
    fprintf(out,"\n");
  }
}
/* 

print a sav[6,6] stifness matrix sorted for Vera's programs

*/
void print_cij_vera_sorted(COMP_PRECISION *sav, FILE *out)
{
  fprintf(out,"%2i %i %7.2f\n",1,1,sav[(1-1)*6]);
  fprintf(out,"%2i %i %7.2f\n",2,2,sav[(2-1)*6+1]);
  fprintf(out,"%2i %i %7.2f\n",3,3,sav[(3-1)*6+2]);
  fprintf(out,"%2i %i %7.2f\n",4,4,sav[(4-1)*6+3]);
  fprintf(out,"%2i %i %7.2f\n",5,5,sav[(5-1)*6+4]);
  fprintf(out,"%2i %i %7.2f\n",6,6,sav[(6-1)*6+5]);
  fprintf(out,"%2i %i %7.2f\n",1,2,sav[(1-1)*6+1]);
  fprintf(out,"%2i %i %7.2f\n",1,3,sav[(1-1)*6+2]);
  fprintf(out,"%2i %i %7.2f\n",1,4,sav[(1-1)*6+3]);
  fprintf(out,"%2i %i %7.2f\n",1,5,sav[(1-1)*6+4]);
  fprintf(out,"%2i %i %7.2f\n",1,6,sav[(1-1)*6+5]);
  fprintf(out,"%2i %i %7.2f\n",2,3,sav[(2-1)*6+2]);
  fprintf(out,"%2i %i %7.2f\n",2,4,sav[(2-1)*6+3]);
  fprintf(out,"%2i %i %7.2f\n",2,5,sav[(2-1)*6+4]);
  fprintf(out,"%2i %i %7.2f\n",2,6,sav[(2-1)*6+5]);
  fprintf(out,"%2i %i %7.2f\n",3,4,sav[(3-1)*6+3]);
  fprintf(out,"%2i %i %7.2f\n",3,5,sav[(3-1)*6+4]);
  fprintf(out,"%2i %i %7.2f\n",3,6,sav[(3-1)*6+5]);
  fprintf(out,"%2i %i %7.2f\n",4,5,sav[(4-1)*6+4]);
  fprintf(out,"%2i %i %7.2f\n",4,6,sav[(4-1)*6+5]);
  fprintf(out,"%2i %i %7.2f\n",5,6,sav[(5-1)*6+5]);
}
