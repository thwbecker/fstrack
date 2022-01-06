#
#
# read lon1 lat1 lon2 lat2 in degrees and give azimuth CW from north in degrees
#
BEGIN{
  pif=57.295779513082320876798154814105;
  twopi=6.2831853071795864769252867665590;
}
{

  if((NF >= 4) && (substr($1,1,1)!="#")){
    lon1=$1/pif;
    lat1=$2/pif;
    lon2=$3/pif;
    lat2=$4/pif;
    if(NF==4)
	printf("%g\n",azimuth(lon1,lat1,lon2,lat2)*pif)
    else
	printf("%g ",azimuth(lon1,lat1,lon2,lat2)*pif)
	
  }
  for(i=5;i<=NF;i++)
      printf("%s ",$i);
  if(NF>4)
      printf("\n");
}



#
# compute azimuth of lon lat in radians
#
function azimuth (lon1, lat1,lon2,lat2) {
  

  slat[1]=sin(lat1);clat[1]=cos(lat1);
  slat[2]=sin(lat2);clat[2]=cos(lat2);
  
  tmp_azi = atan2(sin(lon2-lon1)*clat[2],
		  clat[1]*slat[2]-slat[1]*clat[2]*cos(lon1-lon2));
  if(tmp_azi < 0)
    tmp_azi += twopi;

  return tmp_azi;
}
