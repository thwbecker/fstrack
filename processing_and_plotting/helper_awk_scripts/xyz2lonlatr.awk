#
# reads x y z ...
#
# writes lon lat R ...
#
BEGIN{
#
# convert cartesian system coordinates to lon (phi) lat (latitude, 90- theta)
# and z
# if lower_half==1, will project into 0 ... -90 hemisphere and assign a nagative r value
# 
    pi = 3.1415926535897932384626433832795;
    twopi = 2.0*pi;
    pi = pi/2.0;
    f = 360.0/twopi;
}
{
  # check for GMT -M files
  if($1==">")
    print($1);
  else{
    if((substr($1,1,1)!="#") && (NF>=3)){
      x=$1;
      y=$2;
      z=$3;
      
      tmp1 = x*x + y*y;
      tmp2 = tmp1 + z*z;

      if(tmp2 > 0.0)
	r=sqrt(tmp2);
      else
	r=0.0;
      
      theta=atan2(sqrt(tmp1),z);
      phi=atan2(y,x);

      if(lower_half == 1){
	if(theta <= pihalf){
	  theta = pi - theta;
	  phi += pi;
	}
      }
      if(upper_half == 1){
	if(theta > pihalf){
	  theta = pi - theta;
	  phi += pi;
	}
      }
      while(phi < 0)
	phi += twopi;
      while(phi > twopi)
	phi -= twopi;
      printf("%20.15e  %20.15e  %20.15e ",
	     phi*f,90.-theta*f,r);
      for(i=4;i<=NF;i++)
	printf("%s ",$i);
      printf("\n");
    }
  }
}
END{
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
