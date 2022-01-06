#
# compute the latitude, lat3, and longitude, lon3 of an
# intersection formed by the crs13 true bearing from point 1 and the
# crs23 true bearing from point 2:
#
# input
#
# lon1 lat1 course_13 lon2 lat2 course_23 [all in degree]
#

BEGIN{
  
  pi = 3.141592653589793238;
  twopi = 2.0*pi;
  f = pi/180;
}
{
  if((substr($1,1,1)!="#")&&(NF >= 6)){
      # reads in 
      # first point
      lon[1] = $1*f;		# longitude in deg
      lat[1] = $2*f;		# latitude in deg
      crs13 =  $3*f;		# course from 1 to 3, in degrees
      
      # second point
      lon[2] = $4*f;		# longitude in deg
      lat[2] = $5*f;		# latitude in deg
      crs23 =  $6*f;		# course from 2 to 3, in degrees
 #
 # Returns the point of intersection of two paths defined by point and bearing
 #
 #   see http://williams.best.vwh.net/avform.htm#Intersection
 #

      dLat = lat[2]-lat[1]; dLon = lon[2]-lon[1];
  
      xtmp = sqrt(sin(dLat/2)*sin(dLat/2) + cos(lat[1])*cos(lat[2])*sin(dLon/2)*sin(dLon/2) );
      dist12 = 2.*asin(xtmp);
      if (dist12 == 0) 		# zero distance
	  print("NaN NaN");
      else{
	  # initial/final bearings between points
	  xtmp = ( sin(lat[2]) - sin(lat[1])*cos(dist12) ) / ( sin(dist12)*cos(lat[1]) ) ;
	  brngA = acos(xtmp );
	  #print(brngA) > "/dev/stderr"
	  #if (isNaN(brngA)) brngA = 0;  // protect against rounding
	  xtmp =  ( sin(lat[1]) - sin(lat[2])*cos(dist12) ) / ( sin(dist12)*cos(lat[2]) ) ;
	  brngB = acos(xtmp);
  
	  if (sin(lon[2]-lon[1]) > 0) {
	      brng12 = brngA;
	      brng21 = twopi - brngB;
	  } else {
	      brng12 = twopi - brngA;
	      brng21 = brngB;
	  }
  
	  alpha1 = (crs13 - brng12 + pi) % (twopi) - pi;  # angle 2-1-3
	  alpha2 = (brng21 - crs23 + pi) % (twopi) - pi;  # angle 1-2-3

	  # infinite intersections or ambiguous interaction  
	  if ((sin(alpha1)==0 && sin(alpha2)==0) || 
	      (sin(alpha1)*sin(alpha2) < 0)){
	      print(alpha1,alpha2) > "/dev/stderr"
	      print("NaN NaN") ;
	  }else{
	      
	      #alpha1 = abs(alpha1);
	      #alpha2 = abs(alpha2);
	      # ... Ed Williams takes abs of alpha1/alpha2, but seems to break calculation?
	      
	      xtmp =  -cos(alpha1)*cos(alpha2) + sin(alpha1)*sin(alpha2)*cos(dist12) ;
	      alpha3 = acos(xtmp);
	      dist13 = atan2( sin(dist12)*sin(alpha1)*sin(alpha2), 
			      cos(alpha2)+cos(alpha1)*cos(alpha3) )
	      xtmp = sin(lat[1])*cos(dist13) + cos(lat[1])*sin(dist13)*cos(crs13) ;
	      lat3 = asin(xtmp );
	      dLon13 = atan2( sin(crs13)*sin(dist13)*cos(lat[1]), 
			      cos(dist13)-sin(lat[1])*sin(lat3) );
	      lon3 = lon[1]+dLon13;
	      lon3 = (lon3+3*pi) % (twopi) - pi;  # normalise to -180..+180º
	      
	      print(lon3/f,lat3/f)
	  }
      }
  }
}

function abs( x ) {
    if(x < 0)return -x; else return x;
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

function acos( x ) {
  tmp = 1.0-x*x;
  if(tmp<0)
    tmp = atan2(x,0.0);
  else
    tmp=atan2(x,sqrt(tmp));
  tmp=1.570796326794897 - tmp;
  return(tmp);
}

#
# compute azimuth of lon lat in radians
#
function azimuth (lon1, lat1,lon2,lat2) {
  

  slat[1]=sin(lat1);clat[1]=cos(lat1);
  slat[2]=sin(lat2);clat[2]=cos(lat2);
  
  tmp_azi = atan2(sin(lon2-lon1)*clat[2],
		  clat[1]*slat[2]-slat[1]*clat[2]*cos(lon1-lon2));
  if(tmp_azi < 0)
    tmp_azi += twopi;

  return tmp_azi;
}
function adjust_azi(course, azi) {
    if(course < 0)
	course += twopi;
    if(azi < 0)
	azi += twopi;
    if(azi > course) {
	if(azi - course > pi)
	    course += pi;
    }else{
	if(course - azi > pi)
	    course -= pi;
    }
    return course;
}
	

