#
# given two lon lat locations, compute point at relative distance along 
# great circle
#
# input:
# lon1 lat1 lon2 lat2 frac
#

BEGIN{
    pi = 3.141592653589793238;
    twopi = 2.0*pi;
    f = pi/180;
}
{
    if((substr($1,1,1)!="#")&&(NF>=5)){
	lon1 = $1;
	lat1 = $2;
	lon2 = $3;
	lat2 = $4;
	frac = $5;
	d = distance(lon1,lat1,lon2,lat2,f);
	if(d > 0){
	    lon1r=lon1*f;
	    lat1r=lat1*f;
	    lon2r=lon2*f;
	    lat2r=lat2*f;
	    
	    A=sin((1-frac)*d)/sin(d);
	    B=sin(frac*d)/sin(d);
	    
	    x = A*cos(lat1r)*cos(lon1r) +  B*cos(lat2r)*cos(lon2r);
	    y = A*cos(lat1r)*sin(lon1r) +  B*cos(lat2r)*sin(lon2r);
	    z = A*sin(lat1r)            +  B*sin(lat2r);

	    latr=atan2(z,sqrt(x^2+y^2));
	    lonr=atan2(y,x);
	    #if(lonr < 0)lonr+=twopi;
	}else{
	    lonr=lon1;
	    latr=lat1;
	}
	printf("%20.15e %20.15e ",lonr/f,latr/f);
	for(i=6;i<=NF;i++)
	    printf("%s ",$i)
	printf("\n");
    }
}
END{

}

#
# input in degrees, output in radians
#
function distance(lon1,lat1,lon2,lat2,fac)
{
# d=2*asin(sqrt((sin((lat1-lat2)/2))^2 + 
#                 cos(lat1)*cos(lat2)*(sin((lon1-lon2)/2))^2))
    
  tmplat1=lat1*fac;
  tmplat2=lat2*fac;
  tmplon1=lon1*fac;
  tmplon2=lon2*fac;

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

