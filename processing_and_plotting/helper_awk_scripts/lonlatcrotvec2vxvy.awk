# reads in 
#
# lon lat w_x w_y w_z 
#
# where w_i is a cartesian rotation vector and write velocities to stdout as
# 
# vx vy
#
# if w_i is in deg/Myr, then v_i will be in cm/yr
#
BEGIN{
  f=0.0174532925199433;#pi/180.0;
  pihalf=1.5707963267949;#pi/2.0;
  fvel=11.1194926644559;#R/1.e6*1.e2*f;
  pcnt=0;
}
{
  if((substr($1,1,1)!="#")&&(NF>=5)){
    phi=$1*f;
    theta=pihalf-$2*f;
    wx=$3;
    wy=$4;
    wz=$5;
# coordinates of point in x y z frame with unit length
# radius (velocity factor fvel has the R included)
    sin_theta=sin(theta);cos_theta=cos(theta);
    sin_phi=sin(phi);cos_phi=cos(phi);
    x=cos_phi*sin_theta;
    y=sin_phi*sin_theta;
    z=cos_theta;
# velocities in xyz frame in degrees/my
    vx = wy*z - wz*y;
    vy = wz*x - wx*z;
    vz = wx*y - wy*x;
# velocities in spherical coordinates in rad/my
    vr= sin_theta*cos_phi*vx + sin_theta*sin_phi*vy + cos_theta*vz;
    vt= cos_theta*cos_phi*vx + cos_theta*sin_phi*vy - sin_theta*vz;
    vp=-          sin_phi*vx +           cos_phi*vy;
# 
    printf("%g %g ",vp*fvel,-vt*fvel);
    for(i=6;i<=NF;i++)
      printf("%s ",$i);
    printf("\n");
  }
}
