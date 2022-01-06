# convert 
#
# lon lat mag(radial component) 
# 
# input to 
#
# x y z 
#
# output
#
# $Id: lonlatmag2xyz.awk,v 1.2 2002/02/23 22:41:36 becker Exp $
#
BEGIN{
  f=0.017453292519943296;
}
{
  if((substr($1,1,1)!="#") && (NF>=3) && (tolower($3)!="nan")){
    theta=(90.0-$2)*f;
    phi=$1*f;
    r=$3;
    
    tmp=sin(theta)*r;
    
    x=tmp * cos(phi);
    y=tmp * sin(phi);
    z=cos(theta)*r;
    
    printf("%20.16e %20.16e %20.16e ",x,y,z);
    for(i=4;i<=NF;i++)
      printf("%s ",$i);
    printf("\n");
  }
}
