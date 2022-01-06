{
  if(NR == linenumber){
    
    hour=substr($1,7,2);min=substr($1,9,2);sec=substr($2,1,2);
    lat=$3;if($4=="S")lat *= -1.0;
    lon=$5;if($6=="W")lon *= -1.0;

    year=substr($1,1,2);mon=substr($1,3,2);day=substr($1,5,2);
    date=substr($1,5,2)+"."+substr($1,3,2)+"."+substr($1,1,2);
    mag=$8;
    if(file)
      printf("imgf%2s%2s%2s%2s%2s%2s.gif\n",day,mon,year,hour,min,sec);
    else
      if(region)
	{
	  west  = lon - 15.0;
	  east  = lon + 15.0;
	  south = lat - 15.0;
	  north = lat + 15.0;
	  if(south < -90)south= -90;
	  if(north > 90)north= 90;
	  printf("%g %g %g %g\n",west,east,south,north);
	}
      else
	if(file2)
	  printf("imgf%2s%2s%2s%2s%2s%2s.html\n",day,mon,year,hour,min,sec);
	else
	  printf("%g %g %g\n",lon,lat,mag);
  };
}
