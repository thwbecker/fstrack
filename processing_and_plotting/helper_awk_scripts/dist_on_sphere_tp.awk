#
# read in lon1 lat1 lon2 lat2
# returns distance in radians
#
# lon: longitude (0..360), lat: latitude (-90..90)
#
# JUST FOR TESTING PURPOSES
#
# $Id: dist_on_sphere_tp.awk,v 1.1 2001/08/06 15:41:33 becker Exp $
#
BEGIN{
}
{
  if((NF >= 4)&&((substr($1,1,1)!="#"))){
    phi1=$1*0.0174532925199433;
    phi2=$3*0.0174532925199433;
    theta1=(90.-$2)*0.0174532925199433;
    theta2=(90.-$4)*0.0174532925199433;
    printf("%g ",distancetp(theta1,phi1,theta2,phi2));
    for(i=5;i<=NF;i++)
      printf("%s ",$i);
    printf("\n");
  }
}
END{
}
#
# input in radians, output in radians
#
# theta: colatidude 0..Pi
# phi:   longitude  0...2Pi
#
#
function distancetp(theta1,phi1,theta2,phi2)
{
  tmp1 = (theta2-theta1)/2.0;
  tmp1 = sin(tmp1);
  tmp1 = tmp1 * tmp1;
  
  tmp2 = (phi1-phi2)/2.0;
  tmp2 = sin(tmp2);
  tmp2 = tmp2 * tmp2;

  tmp3  = tmp1;
  tmp3 += sin(theta1) * sin(theta2) * tmp2;
  tmp3 = sqrt(tmp3);
  return 2.0*asin(tmp3);
}

function asin( x ) {
  tmp=atan2(x,sqrt(1.0-x*x+1.0e-17));
  return(tmp);
}
  
 


