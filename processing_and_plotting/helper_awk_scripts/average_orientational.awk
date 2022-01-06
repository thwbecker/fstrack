#
# compute the arithmetic average of an ASCII file with orientational (0...180 deg) data
#
BEGIN{
  f=57.2957795130823;
  x=y=0.0;
  n=0;
}
{
  if((substr($1,1,1)!="#")&&(NF>=1)){
    ttheta = $1/f*2.0;
    x+=sin(ttheta);y+=cos(ttheta);
    n++;
  }
}
END{
  x/=n;y/=n;
  ttmean = atan2(x,y);
  tmean = ttmean/2*f;
# limit to 0..180 range
  if(tmean < 0)
    tmean += 360.0;
  if(tmean > 180)
    tmean -= 180;
  printf("%g\n",tmean);

}
