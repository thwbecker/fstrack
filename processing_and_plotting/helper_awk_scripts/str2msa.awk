#
# convert s_xx s_xy s_yy input to 
#
# s_1, s_2, azi
#
# output. azi is in degrees CW from North
#
BEGIN{

}
{
  if(NF>=3 && (substr($1,1,1)!="#")){
    s11=$1;
    s12=$2;
    s22=$3;
    
    x1 = (s11 + s22)/2.0;
    x2 = (s11 - s22)/2.0;
    r = sqrt(x2 * x2 + s12 * s12 );
    fms = x1 + r ;
    sms = x1 - r;
    deg = 45.0;
    if(x2 != 0.0)
      deg= 22.5*(atan2(s12,x2)/0.78539816339744831);
    else{
      if(s12<0.0)
	deg = deg * -1.0;
    }
# go from clockwise from east (x) to clockwise from north    
    deg = deg + 90.0;
    print(fms,sms,deg);
  }
}
