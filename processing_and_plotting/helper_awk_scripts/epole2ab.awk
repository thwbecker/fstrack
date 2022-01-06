#
# convert lon[deg] lat[deg] deg/My net rotation into torooidal coefficients a10 a11 b11
#
BEGIN{
  f=32.184314659;
  pif=57.295779513082321;	# 180/pi
}
{
  if((substr($1,1,1)!="#")&&(NF<=3)){
    lon=$1;
    lat=$2;

    rot=$3 * f;

    a10 =  sin(lat/pif)*rot;
    a11 = -cos(lon/pif)*cos(lat/pif)*rot;
    b11 = -sin(lon/pif)*cos(lat/pif)*rot;

    printf("%9.5e %9.5e %9.5e\n",a10,a11,b11);
  }

}