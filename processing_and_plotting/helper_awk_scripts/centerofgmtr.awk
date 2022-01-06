#
# determine the center of a GMT region
#
{
  a=substr($1,3);
  split(a,s,"/");
  printf("%g/%g",(s[2]+s[1])/2.0,(s[4]+s[3])/2);


}
