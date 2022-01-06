#
# given lon lat azi and distance[km], compute new location
#
BEGIN{
  pi = 3.141592653589793238;
  twopi = 2.0*pi;
  f = pi/180;
}
{
  if((substr($1,1,1)!="#")&&(NF >= 4)){
      # reads in 
      lon[1] = $1*f;		# longitude in deg
      lat[1] = $2*f;		# latitude in deg
      tc = $3*f;		# azimuth in degree CW from North
      d = $4/6371.0087714;		# distance in km
    
      pfromazi(lon,lat,d,tc);

      
      printf("%20.10e %20.10e ",lon[2]/f,lat[2]/f);
      for(i=5;i<=NF;i++)
	  printf("%s ",$i);
      printf("\n");
  }
  
}
# 
# compute lon[2], lat[2] in distance dist1 from lon[1],lat[1] in direction azi1
# all in radian, azimuth is clockwise from north
#
# need to pass array 


function pfromazi(lon,lat,d,tc)
{
    lat[2] =asin(sin(lat[1])*cos(d)+cos(lat[1])*sin(d)*cos(tc));
    dlon=atan2(sin(-tc)*sin(d)*cos(lat[1]),
	       cos(d)-sin(lat[1])*sin(lat[2]));
    lon[2]=mod(lon[1]-dlon+pi,twopi )-pi;

}




function asin( x ) {
    tmp=atan2(x,sqrt(1.0-x*x));
    return(tmp);
}

function mod (y, x){
  tmp_mod = y - x * int(y/x);
  if ( tmp_mod < 0) 
      tmp_mod += x;
  return tmp_mod;
}
