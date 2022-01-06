# count number of repose intervals larger than
# t
BEGIN{
  n=0;
  xmax=-1e20;
  xmin=1e20;
}
{
  if($1 != "" && (substr($1,1,1)!="#")){
    n++;
    x[n]=$1;
    if(x[n]<xmin)xmin=x[n];
    if(x[n]>xmax)xmax=x[n];
  }
}
END{
  m=0;
  for(x1=0;x1 <= xmax+1;x1+=1.0){
    m++;
    xl[m]=x1;
    b[m]=0;
    for(j=1;j<=n;j++)
      if(x[j] >= x1)
	b[m]++;
  }
  for(i=1;i<=m;i++)
    print(xl[i],b[i],b[i]/n);
}
