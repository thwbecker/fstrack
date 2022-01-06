#
# convert lon lat coordinates to 80N 34S type description of location
#
{
  if(substr($1,1,1)!="#"){
    if(NF>=2){
      lon = $1;
      if(lon > 180)
	lon = -(360-lon);
      if(lon < -180)
	lon = 360+lon;
      if(lon >= 0)
	label1="E";
      else{
	label1="W";
	lon=-lon;
      }
      lat = $2;
      if(lat > 90)
	print("error");
      if(lat < -90)
	print("error");
      if(lat >= 0)
	label2="N";
      else{
	label2="S";
	lat = -lat;
      }
      printf("%.2f%s %.2f%s ",lon,label1,lat,label2);
      for(i=3;i<=NF;i++)
	printf("%s ",$i);
      printf("\n");
    }else{
      print("error");
    }
  }

}
