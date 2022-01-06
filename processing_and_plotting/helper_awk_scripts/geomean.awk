#
# $Id: geomean.awk,v 1.4 2002/05/28 14:04:49 becker Exp $
#
# calculates the area weighted mean of data given in 
#
# lon lat z
#
# format
#
BEGIN{
  pi=3.141592653589793238462643383;
  fac=pi/180.0;
}
{
  if((NF>=3) && (substr($1,1,1)!="#")){
    if(tolower($3)!="nan"){
      w   = cos($2*fac);
      s  += w*$3;
      sw += w;
    }
  }
}
END{
  if(sw != 0.0){
#    printf("%20.10e %20.10e\n",s,sw)
    printf("%20.16e\n",s/sw);
  }
}

