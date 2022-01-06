#
# 
# read in x_i y_i values and interpolates to xint
#
BEGIN{
  n=0;
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
    print(interpolate(x,y,n,xint));
}
function interpolate ( x , y, n , x1 ) {
  j=1;
  while( (j<n) && (x[j] < x1))
    j++;
  if(j==1)
    j=2;
  i=j-1;
  #print(i,j,x[i],x[j])
  fac=(x1 - x[i])/(x[j]-x[i]);
  fac2=1.0-fac;
  tmp=  fac  * y[j] + fac2 * y[i];

  return (tmp);
}
