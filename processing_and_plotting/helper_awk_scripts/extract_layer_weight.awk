#
# given a vedpth.dat file with depths in [km] every row, extract the respective layer thickness
# if by_volume is set to non-zero, will use r^2 weighting
BEGIN{
    n=0;
}
{
  n++;
  z[n]=$1;
  if(z[n]>0)
    z[n] = -z[n];
  r[n] = (6371+z[n])/6371;		# normalized radius

}
END{

  zrange=z[n]-z[1];

  dz[1] = (z[2]-z[1])/2;
  for(i=2;i<n;i++){
    zl = (z[i]-z[i-1])/2;
    zu = (z[i+1]-z[i])/2;
    dz[i]=zu+zl;
  }
  dz[n] = (z[n]-z[n-1])/2;
  if(by_volume){		# by volume
      wsum=0;
      for(i=1;i<=n;i++){
	  dz[i] *= r[i]**2;
	  wsum += dz[i];
      }
      
  }else{
      wsum = zrange;		# linear avg
  }
  for(i=1;i<=n;i++){
    print(dz[i]/wsum);
  }
}
