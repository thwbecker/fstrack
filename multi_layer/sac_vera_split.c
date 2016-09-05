/* 

compute splitting using the levin & park method as implemented by
vera

*/
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>
#include <sacio.h>
#include <string.h>

/* from vera  */
extern void vera_split_(float *,float *,float *,int *,float *,float *,float *,int *);
/* form menke routines */
int menke_distaz(float, float, float, float, float *, float *, float * );
/*  */
int read_sac(float **,int *,float *,float *,char *);

int main(int argc, char **argv)
{
  float *rdata,*tdata,*time,bazi1,bazi2,delta1,delta2;
  int nlen1,nlen2,i;
  float fastaz,deltat,misfit;
  int null;
  
  if(argc != 3){
    fprintf(stderr,"%s radial.sac transverse.sac\n",argv[0]);
    exit(-1);
  }

  if(read_sac(&rdata,&nlen1,&delta1,&bazi1,argv[1])){
    fprintf(stderr,"%s: error: cant open radial file %s for read\n", 
	    argv[0], argv[1]);
    exit(-1);
  }
  if(read_sac(&tdata,&nlen2,&delta2,&bazi2,argv[2])){
    fprintf(stderr,"%s: error: cant open transverse file %s for read\n", 
	    argv[0], argv[2]);
    exit(-1);
  }
  if(nlen1 != nlen2){fprintf(stderr,"length mismatch\n");exit(-1);}
  if(delta1 != delta2){fprintf(stderr,"delta mismatch\n");exit(-1);}
  if(bazi1 != bazi2){fprintf(stderr,"bazi mismatch\n");exit(-1);}
  time = (float *)malloc(sizeof(float)*nlen1);
  for(i=0;i < nlen1;i++)
    time[i] = i * delta1;

  vera_split_(time,rdata,tdata,&nlen2,&fastaz,&deltat,&misfit,&null);
  fastaz = bazi1 + fastaz;

  fprintf(stderr,"bazi fast_azi deltat misfit\n");

  if(fastaz < 0)
    fastaz += 360;
  if(fastaz > 360)
    fastaz -= 360;
  if(fastaz > 180)
    fastaz -= 180;

  if(null)
    printf("%g 0 0 %g\n",bazi1,misfit);
  else
    printf("%g %g %g %g\n",bazi1,fastaz,deltat,misfit);

  free(rdata); free(tdata); free(time);
}



/* 
read the main info from a SAC file

*/
#define RS_MAX 1000000
#define CHECK_ERROR {if(nerr > 0){fprintf(stderr,"read_sac: error %i\n",nerr);exit(nerr);}}

int read_sac(float **data,int *nlen, float *delta, float *baz,char *name)
{
  int nerr , max = RS_MAX;
  float beg , del ,slon,slat,elon,elat,distance,az;
  *data = (float *)malloc(sizeof(float)*RS_MAX);
  rsac1(name, *data, nlen, &beg, &del, &max, &nerr, strlen(name) );
  if(nerr>0)return nerr;
  *data = (float *)realloc(*data,sizeof(float)*(*nlen));
  getfhv ( "DELTA" , delta , &nerr , 5 ) ; CHECK_ERROR; /* get spacing */
  /* compute back-azimuth */
  getfhv ( "STLO" , &slon , &nerr , 5 ) ; CHECK_ERROR; 
  getfhv ( "STLA" , &slat , &nerr , 5 ) ; CHECK_ERROR; 
  getfhv ( "EVLO" , &elon , &nerr , 5 ) ; CHECK_ERROR; 
  getfhv ( "EVLA" , &elat , &nerr , 5 ) ; CHECK_ERROR; 
  /* compute back-azimuth */
  menke_distaz( elat, elon,slat, slon,&distance, &az, baz );
  fprintf(stderr,"read_sac: read %i entries (%g s, dt: %g) from %s, event: %g, %g station %g, %g bazi %g\n",
	  *nlen,*nlen * (*delta), *delta,name,elon,elat,slon,slat,*baz);
  return 0;
}
