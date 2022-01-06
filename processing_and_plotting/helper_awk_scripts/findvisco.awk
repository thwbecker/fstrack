# reads a adfvect type viscosity file and determines the viscosity 
# at depth z[km]
{
  if(NR==1)n=$1;
  else{
    r[NR-1]=$2;
    v[NR-1]=$1*1e21;
    if(NR-1>n)print("error",n,NR-1) > "/dev/stderr";
  }
}
END{
  rc=(6371-z)/6371.;
  i=n;
  while((i>=1)&&(rc<r[i]))
    i--;
  if(i==n)print(v[i]);else print(v[i+1]);
}
