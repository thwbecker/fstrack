#
# compute velocities given rotations vector wx,wy,wz; reads
# 
# lon lat ...
#
# outputs of velocities in 
#
# lon lat v_p -v_t |v| azimuth ...
#
# velocities are in [cm/yr], azimuth is direction CW from North in degrees
#
BEGIN{
  f=0.0174532925199433;#pi/180.0;
  pihalf=1.5707963267949;#pi/2.0;
  R=6371000.0;
  fvel=11.1194926644559;#R/1.e6*1.e2*f;
  pcnt=0;
  
}
{
  if((NF>=2) && (substr($1,1,1) !="#")){
    lon=$1;
    lat=$2;

    phi=lon*f;
    lambda=lat*f;
    theta = pihalf - lambda;

# coordinates of point in x y z frame with unit length
# readius (velocity factor fvel has the R included)
    sin_theta=sin(theta);cos_theta=cos(theta);
    sin_phi  =sin(phi);  cos_phi=cos(phi);
    x = cos_phi * sin_theta;
    y = sin_phi * sin_theta;
    z = cos_theta;

# velocities in xyz frame in degrees/my
    vx=wy*z-wz*y;
    vy=wz*x-wx*z;
    vz=wx*y-wy*x;
# velocities in spherical coordinates in rad/my
    vr= sin_theta*cos_phi*vx + sin_theta*sin_phi*vy + cos_theta*vz;
    vt= cos_theta*cos_phi*vx + cos_theta*sin_phi*vy - sin_theta*vz;
    vp=-          sin_phi*vx +           cos_phi*vy;
# output of velocities in 
#
# lon lat v_p -v_t |v| azimuth 
#
    v=sqrt(vp*vp+vt*vt);
    azimuth=atan2(vp,-vt)/f;
    if(azimuth<0)azimuth += 360.0;
    printf("%g %g %g %g %g %g ",lon,lat,vp*fvel,-vt*fvel,v*fvel,azimuth);
    for(i=3;i<=NF;i++)
      printf("%s ",$i);
    printf("\n");
  }else{
    print($0);
  }
}
