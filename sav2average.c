/* 



$Id: sav2average.c,v 1.3 2007/07/07 00:19:37 becker Exp $


NOT SURE WHAT THIS WAS SUPPOSED TO BE, SEE AVERAGE_RPHI FILES FOR KERNEL AVERAGING


*/
#include "fstrack.h"

int main(int argc, char **argv)
{
  VERA_PREC sav[36],cij[81],rho,x[3],alpha,beta,gamma,*savr[36];
  struct tensor{
    VPREC sav[36];
  };
  int i,j;
  FILE *in;
  if(argc==1){


  }
  
  n=0;
  while(fscanf(in,THREE_FLT_FORMAT,x,(x+1),(x+2))==3){
    if(!read_sym_6by6(sav,in)){
      fprintf(stderr,"%s: error\n",argv[0]);
      exit(-1);
    }
    n++;
  }
  
  fprintf(stderr,"%s: read %i tensors\n",argv[0],n);
}

