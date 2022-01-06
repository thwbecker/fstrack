#
# inserts a zero when greenwich is crossed
#
BEGIN{
  
  am360=360.;
}
{
  lon=$1;
  lat=$2;
  if(NR>1){
    if(lon>330 && oldlon < 30){
      dx=oldlon+(360.0-lon);
      dy=lat-oldlat;
      midlat=oldlat+dy*oldlon/dx;
      printf("%20.10lf %20.10lf\n",0,midlat);
#      printf("%20.10lf %20.10lf\n",firstlon,firstlat);
      print(">");
      firstlon=am360;firstlat=midlat;
      printf("%20.10lf %20.10lf\n",am360,midlat);
      printf("%20.10lf %20.10lf\n",lon,lat);
    }else if(lon<30 && oldlon > 330){
      dx=lon+(360.0-oldlon);
      dy=lat-oldlat;
      midlat=oldlat+dy*(360.0-oldlon)/dx;
      printf("%20.10lf %20.10lf\n",am360,midlat);
#      printf("%20.10lf %20.10lf\n",firstlon,firstlat);
      print(">");
      firstlon=0.0;firstlat=midlat;
      printf("%20.10lf %20.10lf\n",0,midlat);
      printf("%20.10lf %20.10lf\n",lon,lat);
    }else{
      printf("%20.10lf %20.10lf\n",lon,lat);
    }
    oldlon=lon;
    oldlat=lat;
  }else{
    oldlon=lon;
    oldlat=lat;
    firstlon=lon;
    firstlat=lat;
    printf("%20.10lf %20.10lf\n",lon,lat);
  }
}END{
  
#  printf("%20.10lf %20.10lf\n",firstlon,firstlat);

}
