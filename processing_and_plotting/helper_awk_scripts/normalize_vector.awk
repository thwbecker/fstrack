#
# normalize a vector in row form
#
BEGIN{
  n=0; sum=0.0;
}
{
  if((NF>=1) && (substr($1,1,1)!="#")){
    n++;
    x[n]=$1;
    sum +=  x[n]* x[n];
  }
}
END{
  sum=(sum <= 0.0)?(1.0):(sqrt(sum));
  for(i=1;i<=n;i++)
    print(x[i]/sum);
}
