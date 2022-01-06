#
# convert a list of lon lat points to triangle node file format
#
{
  if((substr($1,1)!="#")&&(NF>=2)){
    n++;
    lon[n]=$1;lat[n]=$2;
  }
}
END{
  print(n,2,0,0);
  for(i=1;i<=n;i++){
    print(i,lon[i],lat[i]);


  }
}
