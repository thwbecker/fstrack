#
#
# convert  116w24.79 34n11.01  format of lon lat  to decimals
#
# $Id: lonlatverbose2lonlat.awk,v 1.1 2003/03/11 00:59:02 becker Exp $
#
BEGIN{

}
{
  if((substr($1,1,1)!="#")&&(NF>=2)){
    lon=tolower($1);
    lat=tolower($2);
    if(match(lon,"w")){
      lon_split="w";lon_sign=-1;
    }else if(match(lon,"e")){
      lon_split="e";lon_sign=1;
    }else{
      print("error, no w or e found in lon string") > "/dev/stderr";
    }
    if(match(lat,"s")){
      lat_split="s";lat_sign=-1;
    }else if(match(lat,"n")){
      lat_split="n";lat_sign=1;
    }else{
      print("error, no n or s found in lat string") > "/dev/stderr";
    }
    nlon=split(lon,lonarr,lon_split);
    nlat=split(lat,latarr,lat_split);
    if((nlon != 2)||(nlat != 2)){
      print("format error",lon,lat) > "/dev/stderr";
    }
    longitude = lon_sign*(lonarr[1] + lonarr[2]/60.);
    latitude  = lat_sign*(latarr[1] + latarr[2]/60.);
    printf("%20.10e %20.10e ",longitude, latitude);
    for(i=3;i<=NF;i++)
      printf("%s ",$i);
    printf("\n");

  }
}
