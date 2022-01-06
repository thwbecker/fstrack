#
# calculate mean and standard deviation from col
# if col is not set, use col=1
#
# $Id: meanandstddev.awk,v 1.1 2004/10/28 20:41:27 becker Exp becker $
#
BEGIN {
  sum = sum2 = 0.0;
  n = 0;
  if(col==0)
    col=1;
}
{
  if((substr($1,1,1)!="#") && ($col != "")  && (tolower($col) != "nan")){
    n++;
    x[n] = $col;
    sum += x[n];
  }
}
END {
  if(n){
    mean = sum/n;
    for(i=1;i<=n;i++){
      fac = x[i]-mean;
      fac *= fac;
      sum2 += fac;
    }
    printf("%20.15e %20.15e\n",mean,sqrt(sum2/(n-1)));
  }
}
