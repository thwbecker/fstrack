BEGIN{
  n=0;
  ymax=-1e99;
}
{
  if((substr($1,1,1)!="#")&&(NF>=2)){
    n++;
    x[n]=$1; 
    y[n]=$2;
    if(y[n]>ymax)
      ymax=y[n];
   
  }
}
END{
  for(i=1;i<=n;i++)
    printf("%g %g\n",x[i],y[i]/ymax);
}
