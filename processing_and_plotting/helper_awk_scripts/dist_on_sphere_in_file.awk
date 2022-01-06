#
# for set lon0, lat0 reads in file with 
#
# lon1 lat1 in degrees
#
# returns distance in radians, if km=1, in km using a spherical Earth
#
# lon: longitude(0 ..  360)
# lat: latitude (-90 .. 90)
#
# $Id: dist_on_sphere.awk,v 1.1 2004/08/27 18:57:30 becker Exp $
#
BEGIN{
    if((lon0=="")||(lat0=="")){
	print("error, set lon0 and lat0") > "/dev/stderr"
    }

}
{
  if((NF >= 2)&&((substr($1,1,1)!="#"))){
      printf("%g\t\t",distance(lon0,lat0,$1,$2)*((km==1)?(6371.0):(1.0)));
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

  tmplat1=lat1*0.0174532925199433;
  tmplat2=lat2*0.0174532925199433;
  tmplon1=lon1*0.0174532925199433;
  tmplon2=lon2*0.0174532925199433;

  tmp1=sin((tmplat1-tmplat2)/2.0);
  tmp1=tmp1*tmp1;
  
  tmp2=sin((tmplon1-tmplon2)/2.0);
  tmp2=tmp2*tmp2;
  tmp2*=cos(tmplat1);
  tmp2*=cos(tmplat2);

  tmp3=sqrt(tmp1+tmp2);
  return 2.0*asin(tmp3);
}
  
 
function asin( x ) 
{
  tmp=atan2(x,sqrt(1.-x*x+1.0e-15));
  return(tmp);
}

