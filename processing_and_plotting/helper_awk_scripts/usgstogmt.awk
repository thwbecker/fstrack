{
  lat=$3;
  lon=$5;
  if($4 == "S")lat *= -1;
  if($6 == "W")lon *= -1;
  print(lon,lat,$8/25);
}
