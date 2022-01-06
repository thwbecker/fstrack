#
# for set lon0, lat0 reads in file with 
#
# lon1 lat1 in degrees
#
# returns distance in degrees, and backazimuth from lon0 lat1 to point in degrees
# output is 
#
# lon1 lat1 distance azimuth
#
# lon: longitude(0 ..  360)
# lat: latitude (-90 .. 90)
#
# $Id: dist_on_sphere.awk,v 1.1 2004/08/27 18:57:30 becker Exp $
#
BEGIN{

    pif=57.295779513082320876798154814105;
    twopi=6.2831853071795864769252867665590;

    rlon0=lon0/pif;
    rlat0=lat0/pif;


}
{

  if((NF >= 2)&&((substr($1,1,1)!="#"))){
      rlon1=$1/pif;
      rlat1=$2/pif;
      printf("%s %s %11g %11g \t",$1,$2,
	     distance(rlon0,rlat0,rlon1,rlat1)*pif,
	     azimuth(rlon0,rlat0,rlon1,rlat1)*pif);
    for(i=3;i<=NF;i++)
      printf("%s ",$i);
    printf("\n");
  }
}
END{
}
#
# input in degrees, output in radians
#
function distance(lon1,lat1,lon2,lat2)
{
# d=2*asin(sqrt((sin((lat1-lat2)/2))^2 + 
#                 cos(lat1)*cos(lat2)*(sin((lon1-lon2)/2))^2))


  tmp1=sin((lat1-lat2)/2.0);
  tmp1=tmp1*tmp1;
  
  tmp2=sin((lon1-lon2)/2.0);
  tmp2=tmp2*tmp2;
  tmp2*=cos(lat1);
  tmp2*=cos(lat2);

  tmp3=sqrt(tmp1+tmp2);
  return 2.0*asin(tmp3);
}
  
 
function asin( x ) 
{
  tmp=atan2(x,sqrt(1.-x*x+1.0e-15));
  return(tmp);
}

#
# compute azimuth of lon lat in radians
#
function azimuth (lon1, lat1,lon2,lat2) {
  
  slat[1]=sin(lat1);
  clat[1]=cos(lat1);
  slat[2]=sin(lat2);
  clat[2]=cos(lat2);
  
  tmp_azi = atan2(sin(lon2-lon1)*clat[2],
		  clat[1]*slat[2]-slat[1]*clat[2]*cos(lon1-lon2));
  if(tmp_azi < 0)
    tmp_azi += twopi;

  return tmp_azi;
}
