# major stresses to vector plot
#
# uses x y first_major_stress second_major_stress angle
# input for GMT vector output which is
# x y azimuth (from North) amplitude
#
# assumes extension is positive unless 'geo' flag is set
# 
# angle goes anticlockwise in degrees from the west-east axis
#
# th
#
BEGIN{
# if no scaling factor is given, use unity as scale
  if(factor == 0)
    factor=1.0;
# if geo is set, ie. not equal zero
# stresses are in compression positive format
  if(geo)
    factor *= -1.0;
}
{
  if((substr($1,1,1)!="#") && $1!=""){
# read in data
    x=$1;
    y=$2;
    fms=$3*factor;
    sms=$4*factor;
    deg=90.-$5;
# right angle     
    deg2=deg+90.0;
  
    if(1){
# first major stress is extension
      if(fms > 0){
	print(x,y,deg,fms);
	print(x,y,180.0+deg,fms);
# or compression
      }else{
	fms *= -1.0;
	hfms = fms/2.0;
	degrad=deg*0.017453293;
	sdf=sin(degrad)*hfms;
	cdf=cos(degrad)*hfms;
	x1 = x+sdf;
	y1 = y+cdf;
	x2 = x-sdf;
	y2 = y-cdf;
	deg += 180.0;
	print(x1,y1,deg,fms);
	print(x2,y2,180.0+deg,fms);
      }
# second major stress extension
      if(sms > 0.0){
	print(x,y,deg2,sms);
	print(x,y,180+deg2,sms);
# or compression
      }else{
	sms *= -1.0;
	hsms = sms/2.0;
	deg2rad=deg2*0.017453293;
	sds=sin(deg2rad)*hsms;
	cds=cos(deg2rad)*hsms;

	x1 = x+sds;
	y1 = y+sds;
	x2 = x-sds;
	y2 = y-cds;
	deg2 += 180.0;
	print(x1,y1,deg2,sms);
	print(x2,y2,180.0+deg2,sms);
      }
    }
  }
}

