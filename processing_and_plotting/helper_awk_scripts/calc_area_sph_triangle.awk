#
# calculate the area of a spherical triangle on the Earth's 
# surface, assuming the Earth is a sphere. 
#
# input:
#
# lon1 lat1
# lon2 lat2
# lon3 lat3
#
# output:
#
# area
#
# $Id: calc_area_sph_triangle.awk,v 1.1 2006/08/26 18:33:18 becker Exp becker $
#
BEGIN{
  twopi=6.2831853071795864769252867665590;
  pi = twopi/2.0
  pif=57.2957795130823;
  n=0;
#  R = 6371.;		# turn output into km^2
  R = 1.0;
  t = 0;
}
{
  if((substr($1,1,1)!="#")&&($1!=">")&&(NF>=2)){
    if(n == 3){
      t++;
      print(tri_area(lon,lat)*R**2);
      n=0;
    }
    n++;
    lon[n]=$1/pif;
    lat[n]=$2/pif;
  }
}
END{
 if(n == 3){
   t++;
   print(tri_area(lon,lat)*R**2);
 }
 print("processed ",t,"triangles") > "/dev/stderr";
}

function tri_area (lon, lat)
{
# which one is faster/better?
  mode = 2;

  if(mode == 1){
# use angles between sides of triangle
    a1 = azimuth(lon[1],lat[1],lon[3],lat[3]);
    a2 = azimuth(lon[1],lat[1],lon[2],lat[2]);
    an[1] = a1 - a2;if(an[1]<0)an[1] = -an[1];
    
    a1 = azimuth(lon[2],lat[2],lon[1],lat[1]);
    a2 = azimuth(lon[2],lat[2],lon[3],lat[3]);
    an[2] = a1 - a2;if(an[2]<0)an[2] = -an[2];
    
    a1 = azimuth(lon[3],lat[3],lon[2],lat[2]);
    a2 = azimuth(lon[3],lat[3],lon[1],lat[1]);
    an[3] = a1 - a2;if(an[3]<0)an[3] = -an[3];
    for(i=1;i<=3;i++){
      if(an[i] > pi)
	an[i] = twopi - an[i];
      print(an[i]*pif);
    }
    area = (an[1] + an[2] + an[3] - pi);
  }else{
# use triangle side lengths
    an[1] = dist_on_sphere(lon[1],lat[1],lon[2],lat[2]);
    an[2] = dist_on_sphere(lon[2],lat[2],lon[3],lat[3]);
    an[3] = dist_on_sphere(lon[3],lat[3],lon[1],lat[1]);
    s = (an[1]+an[2]+an[3])/2.0;
    area = 4.0*atan(sqrt(tan(s/2.0)*tan((s-an[1])/2.0)*tan((s-an[2])/2.0)*tan((s-an[3])/2.0)));
  }
  
  return area;
}

#
# compute azimuth of lon lat in radians
#
function azimuth (lon1, lat1,lon2,lat2) {
  

  slat[1]=sin(lat1);clat[1]=cos(lat1);
  slat[2]=sin(lat2);clat[2]=cos(lat2);
  
  azi = atan2(sin(lon2-lon1)*clat[2],
		  clat[1]*slat[2]-slat[1]*clat[2]*cos(lon1-lon2));
  if(azi < 0)
    azi += twopi;

  return azi;
}

#
#  distance on sphere, input in radians
#
function dist_on_sphere(lon1,lat1,lon2,lat2)
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

function atan( x )
{
  tmp = atan2(x, 1.0);
  return tmp;

}

function tan ( x )
{
  return sin(x)/cos(x);

}
