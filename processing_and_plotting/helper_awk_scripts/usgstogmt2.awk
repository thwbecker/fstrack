{
  lat=$3;
  lon=$5;
  if($4 == "S")lat *= -1;
  if($6 == "W")lon *= -1;
  if(lon < -160)
    printf("%g %g 16 0 4 6 %s\n",lon+5,lat-3,NR);
  else
    printf("%g %g 16 0 4 6 %s\n",lon-5,lat-3,NR);
}
