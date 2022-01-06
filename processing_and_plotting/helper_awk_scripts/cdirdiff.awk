# calculate the average deviation in stress directions assuming they are 180\deg
# periodic, input is 
#
# lon lat a1 a2 [weight]
#
# WHERE A1 AND A2 ARE IN DEGREES
#
# if input has five columns, will use the last for weights,
# otherwise, 
#
# IF AREA_WEIGHTS is set, will (additionally) determine weight from 
# latitude
#
# if allow_negative is set, WILL OUTPUT DEVIATIONS AS a2-a1
#
# the random mean would be 45. with sigma of 0.68
#
# if add_lon_lat is set, will print lon lat after weight diff
#
# OUTPUT is
#
# weight angle_diff[deg] [lon lat] ($6,...$NF)
#
#
# $Id: cdirdiff.awk,v 1.1 2005/08/02 17:51:11 becker Exp becker $
#
BEGIN{
  s=ws=0.0;
  if(area_weights=="")
      area_weights = 0;
  if(allow_negative=="")
      allow_negative = 0;
  if(add_lon_lat == "")
      add_lon_lat = 0;
}
{
  if((substr($1,1,1)!="#")&&(NF>=4)){
      x=$1;y=$2;
      if((tolower($3)!="nan") && (tolower($4)!="nan")){
	  for(i=1;i <= 2;i++){
# make sure angles are within [0;360];
	      a[i] = $(2+i);
	      while(a[i] < 0)
		  a[i] += 360.0;
	      while(a[i] > 360.0)
		  a[i] -= 360.0;
	  }
# weighting
	  if(NF >= 5){
	      if(area_weights){
		  weight = $5 * cos(y*0.017453293);
	      }else{
		  weight = $5;
	      }
	  }else{
	      if(area_weights){
		  weight = cos(y*0.017453293);
	      }else{
		  weight = 1.0;
	      }
	  }
# calculate difference in angles assuming that they are Pi periodic
	  da = a[2] - a[1];
	  if(allow_negative){
# allow negative (counter vs. clockwise deviations)
	      if(da > 0){
		  if(da > 180)
		      da -= 180;
		  if(da > 90)
		      da -= 180;
	      }else{
		  if(da<-180)
		      da += 180;
		  if(da<-90)
		      da += 180;
	      }
	  }else{
#
# deviations in both ways count the same
#
	      if(da < 0)
		  da = -da;
	      while(da > 180.0)
		  da -= 180.0;
	      if(da > 90.)
		  da = 180.-da;
	  }
#
# output of weight and angle difference in degrees
#
	  printf("%11g %11g ",weight,da);
#
# output of other columns
#
      }else{
	  printf("0 NaN ");
      }
      if(add_lon_lat)
	  printf("%g %g ",x,y);
      
      if(NF > 5){
	  for(i=6;i<=NF;i++)
	      printf("%s ",$i);
      }
      printf("\n");
  }
}
