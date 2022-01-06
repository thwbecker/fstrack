# convert geographical coordinates from degree:minutes.secondsNWES 
# into decimal format
{
if((substr($1,1,1)!="#")){

  l1=length($1);
  lat_field=substr($1,1,l1-1);
  lat_sign=substr($1,l1,1);
  if(lat_sign == "S")
    lat_sign = -1.0;
  else if(lat_sign == "N")
    lat_sign = 1.0;
  else
    printf("last character of 1. field in line %i is neither N nor S\n",NR);
  split(lat_field,lat,":");
  latitude=lat_sign * (lat[1] + lat[2]/60.0);
  


  l2=length($2);
  lon_field=substr($2,1,l2-1);
  lon_sign=substr($2,l2,1);
   if(lon_sign == "E")
    lon_sign = -1.0;
   else if(lon_sign == "W")
    lon_sign = 1.0;
   else
     printf("last character of 2. field in line %i is neither W nor E\n",NR);
   split(lon_field,lon,":");
   longitude=lon_sign * (lon[1] + lon[2]/60.0);
  
# uncomment for 0...360 instead of -180 ... 180

#   if(longitude<0.0)
#     longitude += 360.0;
   
   print(longitude,latitude);
 
  
  
  
}
}
