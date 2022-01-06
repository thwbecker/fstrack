#
# reads x y z ... 
#
# Euler poles and writes 
#
# lon lat R ...
#
# Euler poles. The difference to a normal xyz2lonlat conversion
# is that the lon range is limited to 0 ... 180, with flip of the pole 
# magnitude sign
#
BEGIN{

  twopi = 6.28318530717958647;
  f = 360.0/twopi;
}
{
  # check for GMT -M files
  if($1==">")
    print($1);
  else
    if((substr($1,1,1)!="#") && (NF>=3)){
      x=$1;
      y=$2;
      z=$3;
      
      tmp1 = x*x + y*y;
      tmp2=tmp1 + z*z;
      if(tmp2 > 0.0)
	r=sqrt(tmp2);
      else
	r=0.0;
      theta=atan2(sqrt(tmp1),z);
      phi=atan2(y,x);
      if(phi < 0)
	phi += twopi;
      if(phi >= twopi)
	phi -= twopi;
      lon = phi *f;
      lat = 90 - theta*f;
      mag = r;
      if(lon > 180){
	lon -= 180;
	lat = -lat;
	mag = -mag;
      }
      if(lon < 0){
	lon += 180;
	lat = -lat;
	mag = -mag;
      }
      printf("%11g %11g %11g ",lon,lat,mag)
      for(i=4;i<=NF;i++)
	printf("%s ",$i);
      printf("\n");
    }
}
END{
}
