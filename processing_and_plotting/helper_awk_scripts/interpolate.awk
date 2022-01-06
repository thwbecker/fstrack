#
# 
# read in x_i y_i values and interpolates from xmin to xmax
#
BEGIN{
  n=0;
  if(xmax==0)
    xmax=2850.0;
  if(xmin==0)
    xmin=50;
  if(dx==0)
    dx=50.0;
}
{
  if((substr($1,1,1)!="#")&&(NF>=2)){
    n++;
    x[n]=$1;
    y[n]=$2;
    if(n>=2){
      if(x[n]<x[n-1]){
	printf("error, need sorted x values") > "/dev/stderr";
      }
    }
  }
}
END{
  for(x1=xmin;x1<=xmax+1e-5;x1+=dx)
    print(x1,interpolate(x,y,n,x1));
}
function interpolate ( x , y, n , x1 ) {
  j=1;
  while( (j<n) && (x[j] < x1))
    j++;
  if(j==1)
    j=2;
  i=j-1;

  fac=(x1 - x[i])/(x[j]-x[i]);
  fac2=1.0-fac;
  tmp=  fac  * y[j] + fac2 * y[i];

  return (tmp);
}
