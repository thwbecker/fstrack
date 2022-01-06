#
# determine the horizontal center of a GMT region
#
{
  a=substr($1,3);
  split(a,s,"/");
  print((s[2]+s[1])/2.0);


}
