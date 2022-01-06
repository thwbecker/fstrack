#
# convert lon lat r(normalized) input to x y z output
# input:
#
# lon lat r [ v1 v2 v3 ...]
#
# output:
#
# x y z     [ v1 v2 v3 ...]
#
# $Id: lonlat2xyz.awk,v 1.3 2003/12/01 02:42:50 becker Exp $
#
BEGIN{
  f=0.017453292519943296;
}
{
  if((substr($1,1,1)!="#") && (NF>=3)){
    lambda=$2*f;
    phi=$1*f;
    r = $3;
    
    tmp=cos(lambda)*r;
    
    x=tmp * cos(phi);
    y=tmp * sin(phi);
    z=sin(lambda)*r;
    
    
    printf("%20.16e %20.16e %20.16e ",x,y,z);
    for(i=4;i<=NF;i++)
      printf("%s ",$i);
    printf("\n");
  }
}
