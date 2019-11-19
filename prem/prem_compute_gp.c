/* 

compute the gravity and pressure in PREM
by performing a simple (midpoint) integration

*/
#include "prem.h"
#define CPREC double

int main(int argc, char **argv)
{
  struct prem_model prem;
  CPREC drhodr,dr,*g,*rho,sum,G,gfac,r2,rl,rfac,p;
  int nr,i,nprint;
  FILE *out;
  /*  */
  nr = 100000;			/* number of steps */
  nprint = 100;			/* print every nprint */
  /* constant */
  G = 6.6726e-11;		/* G constant */
  gfac = 3.14159265358979323846264 * 4.0 * G;
  rfac = PREM_RE_KM * 1e3;
  /* 
     read prem 
  */
  prem_read_model(PREM_MODEL_FILE,&prem,TRUE);
  /* 
     need to save g and rho 
  */
  g = (CPREC *)calloc(nr,sizeof(CPREC));
  rho = (CPREC *)calloc(nr,sizeof(CPREC));
  if(!g || !rho){
    fprintf(stderr,"memerror\n");exit(-1);
  };
  /* step size */
  dr = 1./((CPREC)nr);
  /* 
     compute integral for g at mid points 
  */
  sum = 0.0;
  out = fopen("prem.g.dat","w");
  for(rl=dr/2,i=0;i<nr;rl+=dr,i++){
    prem_get_rhodrho((rho+i),&drhodr,rl,&prem);
    r2 = rl*rl;
    sum += rho[i] * r2 * dr;
    /* m(r) = int_0^r 4 \pi r^2 \rho(r) dr */
    /* g(r) = G m(r)/r^2
            = 4\piG/r^2 \int_0^r r^2 \rho(r) dr
    */
    g[i] = gfac / r2 * sum * rfac;
    if(i%nprint == 0)
      fprintf(out,"%g %g\n",rl*rfac/1e3,g[i]);
  }
  fclose(out);
  /* 
     integrate pressure 
  */
  p = 1e5;			/* at surface */
  rl = 1;
  out = fopen("prem.p.dat","w");
  fprintf(out,"%g %g\n",rl*rfac/1e3,p);
  for(i=1;i<nr-1;i++){
    /* dp/dr = -g(r) rho(r) */
    rl -= dr;
    p += dr * rfac * g[i] * rho[i];
    if(i % nprint == 0)
      fprintf(out,"%g %g\n",rl*rfac/1e3,p);
  }
  fclose(out);

  return 0;
}
