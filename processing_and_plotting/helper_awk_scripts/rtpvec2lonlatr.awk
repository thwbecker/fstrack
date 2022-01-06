#
# convert r t p vector to lonlat r
#
BEGIN{
  twopi = 6.28318530717958647;
  pi = twopi / 2.0;
  pihalf = twopi/4.0;
  f = 360.0/twopi;

}
{
  if($1==">")
    print($0);
  else{
    if((substr($1,1,1)!="#") && (NF>=3)){
      # vector components
      x = $2;			# theta
      y = $3;			# phi
      z = $1;			# r
      
      tmp1 = x*x + y*y;
      tmp2=tmp1 + z*z;

      if(tmp2 > 0.0)
	r=sqrt(tmp2);
      else
	r=0.0;
      # 
      theta=atan2(sqrt(tmp1),z);
      phi=atan2(y,x);
      
      if(lower_half == 1){	# limit to lower half
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
