#
# read in lon0 lat0 lon1 lat1 in degrees
#
# returns distance in radians, if km=1, in km using a spherical Earth
#
# lon: longitude(0 ..  360)
# lat: latitude (-90 .. 90)
#
# $Id: dist_on_sphere.awk,v 1.1 2004/08/27 18:57:30 becker Exp becker $
#
BEGIN{
    R_EARTH = 6371.0087714;
    if(km == 1)
	fac = R_EARTH;
    else if(deg == 1)
	fac = 57.295779513082320876798154814105
    else
	fac = 1;
    if(robust=="")
	robust = 1;		# use faster (?) method or more robust?
	
}
{
  if((NF >= 4)&&((substr($1,1,1)!="#"))){
      d = distance($1,$2,$3,$4,robust)*fac; # using the arcsin function
      printf("%20.17e ",d);
      for(i=5;i<=NF;i++)
	  printf("%s ",$i);
      if(repeat_coordinates)
	  printf("%s %s %s %s ",$1,$2,$3,$4);
      printf("\n");
  }
}
END{
}
#
# input in degrees, output in radians
#
function distance(lon1,lat1,lon2,lat2,robust)
{
    tmplat1=lat1*0.0174532925199433;
    tmplat2=lat2*0.0174532925199433;
    tmplon1=lon1*0.0174532925199433;
    tmplon2=lon2*0.0174532925199433;
    dist_lat = tmplat1-tmplat2;
    dist_lon = tmplon1-tmplon2;
    if(robust){
	# this one is better for antipodes
	cos_phi1 = cos(tmplat1);
	cos_phi2 = cos(tmplat2);
	sin_phi1 = sin(tmplat1);
	sin_phi2 = sin(tmplat2);
	cos_dlambda = cos(dist_lon);
	
	tmp1 = cos_phi2 * sin(dist_lon);
	tmp1 *= tmp1;
	tmp2 = cos_phi1*sin_phi2 - sin_phi1*cos_phi2*cos_dlambda;
	tmp2 *= tmp2;
	
	tmp3 = sin_phi1*sin_phi2 + cos_phi1*cos_phi2*cos_dlambda;

	return atan2(sqrt(tmp1+tmp2),tmp3);
    }else{
	# but tihs one is OK up to 1e-14
	# d=2*asin(sqrt((sin((lat1-lat2)/2))^2 + cos(lat1)*cos(lat2)*(sin((lon1-lon2)/2))^2))
	tmp1=sin(dist_lat/2.0);
	tmp1 *=tmp1;
	
	tmp2=sin(dist_lon/2.0);
	tmp2 *= tmp2;
	tmp2 *= cos(tmplat1);
	tmp2 *= cos(tmplat2);
	
	tmp3=sqrt(tmp1+tmp2);
	return 2.0*asin(tmp3);
    }
}
function asin( x ) 
{
    #tmp=atan2(x,sqrt(1.-x*x+1.0e-16)); # not sure why I used that, seems stable
    tmp=atan2(x,sqrt(1.-x*x));
    return(tmp);
}

