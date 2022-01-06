# reads from two files:
#
# a) rotation poles, determined by .rotvec.dat suffix:
#    euler vectors, format:
#    plate wx wy wz [deg/My]
#
# b) points with two plate codes, format:
#    lon lat plate_1 (fixed) plate_2 (moving with respect to plate 1)
#
# output of velocities in 
#
#     lon lat v_p -v_t plate_1 plate_2 |v| azimuth 
#
# velocities are in [cm/yr], azimuth is direction CW 
# from North in degrees
#
BEGIN{
  f=0.0174532925199433;#pi/180.0;
  pihalf=1.5707963267949;#pi/2.0;
  R=6371000.0;
  fvel=11.1194926644559;#R/1.e6*1.e2*f;
  pcnt=0;
}
{
  namepart=substr(FILENAME,length(FILENAME)-9,length(FILENAME));
  
  if(namepart=="rotvec.dat"){
# read in the euler poles
    if(($1 != "") && (substr($1,1,1) !="#")){
      pcnt++;
# reading the rotation vecs 
      name[pcnt]=$1;
      wx[pcnt]=$2;
      wy[pcnt]=$3;
      wz[pcnt]=$4;
    }
  }else{
# reading the points and plate codes    
    if(pcnt == 0){
      print("WARNING: no plates yet read") > "/dev/stderr";
    }
    if((NF>=4) && (substr($1,1,1) !="#")){
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
# plates 
      plate1=$3;plate2=$4;
# search for plates
      fo1=0;fo2=0;
      for(i=1;i<=pcnt;i++){
	if(name[i]==plate1){
	  wx1=wx[i];wy1=wy[i];wz1=wz[i];
	  fo1=1;
	}
	if(name[i]==plate2){
	  wx2=wx[i];wy2=wy[i];wz2=wz[i];
	  fo2=1;
	}
      }
      if((fo1==0) || (fo2==0)){
	print("could not find plates") > "/dev/stderr";
      }else{
# relative Euler vector
	wxr=wx2-wx1;
	wyr=wy2-wy1;
	wzr=wz2-wz1;
#	print("rel vector:",wxr,wyr,wzr) > "/dev/stderr";
# velocities in xyz frame in degrees/my
	vx=wyr*z-wzr*y;
	vy=wzr*x-wxr*z;
	vz=wxr*y-wyr*x;
# velocities in spherical coordinates in rad/my
	vr= sin_theta*cos_phi*vx + sin_theta*sin_phi*vy + cos_theta*vz;
	vt= cos_theta*cos_phi*vx + cos_theta*sin_phi*vy - sin_theta*vz;
	vp=-          sin_phi*vx +           cos_phi*vy;
# output of velocities in 
#
# lon lat v_p -v_t plate_1 plate_2 |v| azimuth ...
#
# format
	v=sqrt(vp*vp+vt*vt);
	azimuth=atan2(vp,-vt)/f;
	if(azimuth < 0)
	    azimuth += 360.0;
	printf("%g %g %g %g %s %s %g %g ",
	       lon,lat,vp*fvel,-vt*fvel,plate1,plate2,
	       v*fvel,azimuth);
	for(i=5;i<=NF;i++)
	  printf("%s ",$i);
	printf("\n");
      }
    }}
}
