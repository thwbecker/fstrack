#
# calculate the minimum, maximum, mean, and standard, and N
# deviation of col
#
# $Id: minmaxmeanstd.awk,v 1.1 2002/12/27 22:01:31 becker Exp becker $
#
BEGIN {
  sum = sum2 = 0.0;
  n = 0;
  min = 1e30;
  max =-1e30;

  if(col==0)
    col=1;
}
{
  if((NF>=col)&&(substr($1,1,1)!="#")&&(tolower($col)!="nan")){
    sum += $col;
    sum2 += $col * $col;
    if(min>$col)min=$col;
    if(max<$col)max=$col;
    n++;
  }
}
END {
  if(n){
    if(n>1)
      std= sqrt ((n * sum2 - sum * sum) / ((n*(n-1))));
    else
      std="nan";
    mean = sum/n;
    print(min,max,mean,std,n);
  }
}
