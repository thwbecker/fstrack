#
# reads in lon lat data and finds duplicates
#
BEGIN{
  n=0;
}
{
  if((substr($1,1,1)!="#")&&(NF>=2)){
    n++;
    lon[n]=$1;lat[n]=$2;
    line[n]=$0;
  }
}
END{
  for(i=1;i<=n;i++)
    for(j=1;j<i;j++){
      dlon = lon[i] - lon[j];
      if(dlon < 0)
	dlon = -dlon;
      if(dlon > 180.0)	
	dlon = 360.0 - dlon;
      if(dlon < 1e-4){
	dlat = lat[i] - lat[j];
	if(dlat < 0)
	  dlat = -dlat;
	if(dlat < 1e-4){
	  print("found double",lon[i],lat[i],lon[j],lat[j],i,j);
	}
      }
    }
}
