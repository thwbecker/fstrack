BEGIN{
#
# split polygon if distance between points is larger than 
# dlim (in km, needs to be specified) is 
#
  twopi = 6.2831853071795864769252867;
  if(dlim == "")
    dlim = twopi;
  else
    dlim /= 6371;
  n=0;
}
{
  if((substr($1,1,1)!="#")&&($1!="")){
    if($1 == ">")
      print($0);
    else{
      if(NF>=2)
	n++;
      if(n==1){
	print($0);
      }else{
	d = distance($1,$2,xold,yold);
	if(d > dlim)
	  print(">");
	print($0);
      }
      xold=$1;yold=$2;
    }
  }
}
#
# input in degrees, output in radians
#
function distance(lon1,lat1,lon2,lat2)
{
# d=2*asin(sqrt((sin((lat1-lat2)/2))^2 + 
#                 cos(lat1)*cos(lat2)*(sin((lon1-lon2)/2))^2))

  tmplat1=lat1*0.0174532925199433;
  tmplat2=lat2*0.0174532925199433;
  tmplon1=lon1*0.0174532925199433;
  tmplon2=lon2*0.0174532925199433;

  tmp1=sin((tmplat1-tmplat2)/2.0);
  tmp1=tmp1*tmp1;
  
  tmp2=sin((tmplon1-tmplon2)/2.0);
  tmp2=tmp2*tmp2;
  tmp2*=cos(tmplat1);
  tmp2*=cos(tmplat2);

  tmp3=sqrt(tmp1+tmp2);
  return 2.0*asin(tmp3);
}
  
 
function asin( x ) 
{
  tmp=atan2(x,sqrt(1.-x*x+1.0e-15));
  return(tmp);
}

