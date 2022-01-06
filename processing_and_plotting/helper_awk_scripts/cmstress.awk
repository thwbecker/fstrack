BEGIN{
}
{
  x=$1;
  y=$2;
  s11=$3;
  s12=$4;
  s22=$5;
  
  x1 = (s11 + s22)/2.0;
  x2 = (s11 - s22)/2.0;
  r = sqrt(x2 * x2 + s12 * s12 );
  fms = x1 + r ;
  sms = x1 - r;
  deg = 45.0;
  if(x2 !=  0.0)
    deg= 22.5*(atan2(s12,x2)/atan2(1.0,1.0));

  printf("%g %g %g %g %g\n",x,y,fms,sms,deg);
}
