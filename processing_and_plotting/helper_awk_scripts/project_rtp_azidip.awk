#
# convert r,theta,phi UNITY VECTOR into local azimuth 
# and dip coordinates. if lower_half == 1, make sure it's the 
# lower hemisphere, likewise for upper_half==1
#
# in degrees, dip is the same as latitude
#
# note that this assumes that azimuth goes from 0 (N) CW 
# to 90 E 180 S ... 
#
BEGIN{
  f=57.2957795130823;
}
{
  if((substr($1,1,1)!="#")&&(NF>=3)){
    er=$1;
    et=$2;
    ep=$3;
    azi = atan2(ep,-et)*f;
    if(azi<0)
      azi += 360.0;
    dip = asin(er)*f;		# colatitude (theta) would be acos 
    if(lower_half == 1){
      if(dip > 0){
	dip = -dip;
	azi += 180;
      }
    }
    if(upper_half == 1){
      if(dip < 0){
	dip = -dip;
	azi += 180;
      }
    }
    while(azi < 0.0)
      azi += 360.0;
    while(azi > 360.0)
      azi -= 360.0;
    print(azi,dip);
  }
}
END{

}
function asin( x ) {
  xx = x * x;
  if(xx < 1)
    tmp=atan2(x,sqrt(1.0-x*x));
  else
    tmp=atan2(x,0);
  return(tmp);
}

