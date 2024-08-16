#include "fstrack.h"
/* 

   convert a file with 

   back_azi[deg] fast_azi[deg] dt[s]

   into splitting statistics, can select if delay times are used for
   averaging or not (default is not)
   
   fazi2splitstat file nharm use_dt_for_avg

*/
#define MAX_NR_HARM 10
int main(int argc, char **argv)
{
  FILE *in;
  int n,nan,i;
  COMP_PRECISION *x=NULL,*fazi=NULL,*dt=NULL,mean_fazi,std_fazi,
    mean_dt,std_dt, /* magnitude of the fitting terms */
    *ma_fazi,*msa_fazi,
    *a_fazi,*a_dt,
    *ma_dt,*msa_dt;
  int nharm = 6,npara,use_dt_for_avg = 0;
  static my_boolean verbose = FALSE;
  
  if((argc>1 && strcmp(argv[1],"-h")==0)||(argc>4)){
    fprintf(stderr,"%s: usage\n%s [stdin/FILENAME] [nharm, %i] [use_dt_for_avg, %i]\n\n",argv[0],argv[0],nharm,use_dt_for_avg);
    exit(-1);
  }
  if((argc > 1)&&(strcmp(argv[1],"stdin")!=0))
    in = myopen(argv[1],"r",argv[0]);
  else
    in=stdin;
  if(argc > 2)
    sscanf(argv[2],"%i",&nharm);
  if(argc > 3)
    sscanf(argv[3],"%i",&use_dt_for_avg); /*  */
  if(nharm > MAX_NR_HARM){
    fprintf(stderr,"%s: too many harmonics (%i)\n",
	    argv[0],nharm);
    exit(-1);
  }
  n=nan=0;
  my_vecrealloc(&x,1,argv[0]);
  my_vecrealloc(&fazi,1,argv[0]);
  my_vecrealloc(&dt,1,argv[0]);
  while(fscanf(in,THREE_FLT_FORMAT,(x+n),(fazi+n),(dt+n))==3){
    if(!finite(fazi[n]))nan++;
    n++;
    my_vecrealloc(&x,(n+1),argv[0]);
    my_vecrealloc(&fazi,(n+1),argv[0]);
    my_vecrealloc(&dt,(n+1),argv[0]);
  }
  if(verbose){
    fprintf(stderr,"%s: read in %i split measurements (%i nan)\n",
	    argv[0],n,nan);
    if(use_dt_for_avg)
      fprintf(stderr,"%s: using delay time for averaging\n",argv[0]);
    else
      fprintf(stderr,"%s: not using delay time for averaging\n",argv[0]);
  }
  /* for fitting */
  npara = 2 + nharm * 2;
  my_vecalloc(&a_fazi, npara,"");
  my_vecalloc(&a_dt,   npara,"");
  my_vecalloc(&ma_fazi,nharm,"");
  my_vecalloc(&msa_fazi,nharm,"");
  my_vecalloc(&ma_dt,  nharm,"");
  my_vecalloc(&msa_dt,  nharm,"");

  analyze_splitting_azi_dep(n,nharm,x,fazi,dt,
			    &mean_fazi,&std_fazi,
			    &mean_dt,&std_dt,
			    a_fazi,a_dt,
			    ma_fazi,msa_fazi,ma_dt,msa_dt,
			    TRUE,use_dt_for_avg);
  if(verbose){
    /* 
       print mean fast azimuth and splitting time as well as std. 
    */
    fprintf(stderr,"%s: fast_azi: %6.1f (+/-%5.1f%% of 180) dt: %6.2fs (+/-%5.1f%%)\n",
	    argv[0],mean_fazi,std_fazi/180.*100.,
	    mean_dt,std_dt/mean_dt*100.);
    /* harmonic terms */
    fprintf(stderr,"%s: fit to dazi: ",argv[0]);
    for(i=0;i < nharm;i++)	/* normalized strength of terms and
				   relative to uncertainties */
      fprintf(stderr,"%1it: %6.2f (%6.2f) ",
	      i+1,ma_fazi[i],log(ma_fazi[i]/msa_fazi[i]));
    fprintf(stderr,"\n");
    fprintf(stderr,"%s: fit to   dt: ",argv[0]);
    for(i=0;i < nharm;i++)
      fprintf(stderr,"%1it: %6.2f (%6.2f) ",
	      i+1,ma_dt[i],log(ma_dt[i]/msa_dt[i]));
    fprintf(stderr,"\n");
  }
  printf("%6.2f %6.2f\t %6.2f %6.4f\t",
	 mean_fazi,std_fazi,mean_dt,std_dt);
  for(i=0;i < nharm;i++)
    printf("%8.4f ",ma_fazi[i]);
  printf("\n");

}
