#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <rpc/rpc.h>
#include "SHahhead.h"
#include <sacio.h>
#include "proto.h"
char *progname;

int read_sac(ahhed *,float **,char *);



int main(int argc, char **argv)
{
  ahhed head1, head2;
  float *data1, *data2;
  int i;
  if( argc != 3 ) {
    fprintf(stderr,"usage: %s radial_sac_file tangential_sac_file\n",argv[0]);
    exit(-1);
  }
  progname=argv[0];
  /* input traces */
  if(read_sac(&head1,&data1,argv[1])){
    fprintf(stderr,"error: %s: cant open radial (or north) file <%s> for read\n", argv[0], argv[1]);
    exit(-1);
  }
  if(read_sac(&head2,&data2,argv[2])){
    fprintf(stderr,"error: %s: cant open tangential (or east) file <%s> for read\n", argv[0], argv[2]);
    exit(-1);
  }
  for(i=0;i< head1.record.ndata;i++)
    printf("%g %g %g\n",head1.record.delta*i,data1[i],data2[i]);
  
  
  free(data1); free(data2); 
}


