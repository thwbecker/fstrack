#
# lon0 lat0 lon1 lat1 in degrees
#
# returns distance in degrees, 
# azimuth from lon0 lat0 to lon1 lat1
# and mid point location
#
# output is 
#
# distance[deg] azimuth[CW from N in degree] lonm latm
#
#
#
BEGIN{

    pif=57.295779513082320876798154814105;
    twopi=6.2831853071795864769252867665590;

    frac  = 0.5;
}
{

  if((NF >= 4)&&((substr($1,1,1)!="#"))){

      rlon0=$1/pif;
      rlat0=$2/pif;

      rlon1=$3/pif;
      rlat1=$4/pif;

      d = distance(rlon0,rlat0,rlon1,rlat1);
      if(d<1e-14){
	  print("zero distance for line", NR,d) > "/dev/stderr"
	  printf("NaN NaN NaN NaN ")
      }else{
	  A=sin((1-frac)*d)/sin(d);
	  B=sin(frac*d)/sin(d);
	  
	  x = A*cos(rlat0)*cos(rlon0) +  B*cos(rlat1)*cos(rlon1);
	  y = A*cos(rlat0)*sin(rlon0) +  B*cos(rlat1)*sin(rlon1);
	  z = A*sin(rlat0)            +  B*sin(rlat1);
	  
	  latm=atan2(z,sqrt(x^2+y^2))*pif;
	  lonm=atan2(y,x)*pif;
	  if(lonm<0)
	      lonm+=360;

	  printf("%11g %11g %g %g\t",
		 d*pif,azimuth(rlon0,rlat0,rlon1,rlat1)*pif,
		 lonm, latm);
      }
      for(i=5;i<=NF;i++)
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
