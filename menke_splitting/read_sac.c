#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <sacio.h>
#include <string.h>

#include "SHahhead.h"
int read_sac(ahhed *,float **,char *);

#define RS_MAX 100000
#define CHECK_ERROR {if(nerr > 0){fprintf(stderr,"read_sac: error %i\n",nerr);exit(nerr);}}



/* 


read the main info from a SAC file

*/
int read_sac(ahhed *head,float **data,char *name)
{
  float beg , del ,delta,slon,slat,elon,elat;
  int nlen , nerr , max = RS_MAX ;
  
  *data = (float *)malloc(sizeof(float)*RS_MAX);
  rsac1(name, *data, &nlen, &beg, &del, &max, &nerr, strlen(name) ) ;
  if(nerr>0)
    return nerr;
  *data = (float *)realloc(*data,sizeof(float)*nlen);
  getfhv ( "DELTA" , &delta , &nerr , 5 ) ; CHECK_ERROR; /* get spacing */
  getfhv ( "STLO" , &slon , &nerr , 5 ) ; CHECK_ERROR; 
  getfhv ( "STLA" , &slat , &nerr , 5 ) ; CHECK_ERROR; 
  getfhv ( "EVLO" , &elon , &nerr , 5 ) ; CHECK_ERROR; 
  getfhv ( "EVLA" , &elat , &nerr , 5 ) ; CHECK_ERROR; 
  head->event.lat = elat;
  head->event.lon = elon;
  head->station.slat = slat;
  head->station.slon = slon;
  head->record.ndata = nlen;
  head->record.delta = delta;

  return 0;
}
