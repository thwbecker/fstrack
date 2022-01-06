# reads in euler poles in 
#  plate lon lat magnitude format
# output is 
#   plate lon lat magnitude wx wy wz
BEGIN{
  pi=3.14159265358979323;
  f=pi/180.0;
}
{
  if((substr($1,1,1)!="#")){
    plate=$1;
    lon=$2;
    lon *= f;
    lat=$3;
    lat *= f;
    w=$4;

    wz=w*sin(lat);
    p =w*cos(lat);
    wx=cos(lon)*p;
    wy=sin(lon)*p;
    for(i=1;i<=NF;i++)
      printf("%s ",$i);
    printf("%12.6e %12.6e %12.6e\n",wx,wy,wz);
  }
}

