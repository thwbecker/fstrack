#
# calculates major stresses
# expects that extension is positive, input format:
#
# input:
#
# s11 s12 s22 ....
#
# output:
#
# fms sms azi .... 
#
# first major stress fms will be the most extensive
# second major stress the least extensive, most compressive
#
# azi is the angle to the first major stress axis (most extensive)
# in degrees clockwise from north (direction 2, y-axis)
#
# $Id: calcms.awk,v 1.1 2004/11/01 20:07:21 becker Exp becker $
#
BEGIN{
#
# 90/pi from 2\theta = atan((2\sigmax_{xy})(\sigma_{xx}-\sigma_{yy})
#
  fac = 28.64788975654116;
}
{
  if((substr($1,1,1)!="#")&&(NF>=3)){
    s11=$1;
    s12=$2;
    s22=$3;
    x1 = (s11 + s22)/2.0;
    x2 = (s11 - s22)/2.0;
    r = x2 * x2 + s12 * s12;
    if(r > 0.0){
      r = sqrt(r);
      fms = x1 + r;
      sms = x1 - r;
    }else{
      fms = sms = x1;
    }
# first method
    if(x2 != 0.0)
      deg = fac * atan2(s12,x2);
    else if(s12 <= 0.0)
      deg= -45.0;
    else
      deg=  45.0;
    azi = 90.0 - deg;

    printf("%g %g %g ",fms,sms,azi); 
    for(i=4;i<=NF;i++)
      printf("%s ",$i);
    printf("\n");
  }
}
