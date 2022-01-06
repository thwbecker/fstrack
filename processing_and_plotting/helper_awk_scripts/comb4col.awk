BEGIN{
  m=2*n;
  o=3*n;
}
{
if(NR>o)
  v[NR-o]=$1;
else if(NR>m)
  z[NR-m]=$1;
else if(NR>n)
  y[NR-n]=$1;
else
  x[NR]=$1;
}
END{
  for(i=1;i<=n;i++)
    print(x[i],y[i],z[i],v[i]);

}
