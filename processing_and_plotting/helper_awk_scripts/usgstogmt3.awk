{
  hour=substr($1,7,2);min=substr($1,9,2);sec=substr($2,1,2)
  lat=$3;if($4=="S")lat *= -1.0;
  lon=$5;if($6=="W")lon *= -1.0;
  year=substr($1,1,2);mon=substr($1,3,2);day=substr($1,5,2);
  date=substr($1,5,2)+"."+substr($1,3,2)+"."+substr($1,1,2);
  printf("<tr align=center><td>%s</td><td><a href=http://www.geophysik.uni-frankfurt.de/seismograms/%02d%02d%02dtnsseismogram.gif>%s.%s.%s</a></td><td>%s:%s:%s</td>\n",NR,day,mon,year,day,mon,year,hour,min,sec);
  printf("<td><a href=http://pubweb.parc.xerox.com/map/color=1/ht=90.00/lat=%s/lon=%s/mark=%s,%s>%6s%1s, %6s%1s</a></td><td>%4s</td><td> %3s</td><td><a href=http://www.geophysik.uni-frankfurt.de/topomaps/imgf%2s%2s%2s%2s%2s%2s.html>%1s%s %1s%s %1s%s</a></td> </tr>\n",
	 lat,lon,lat,lon,
	 $3,$4,$5,$6,$7,$8,
	 day,mon,year,hour,min,sec,
	 substr($9,1,1),tolower(substr($9,2,length($9))),
	 substr($10,1,1),tolower(substr($10,2,length($10))),substr($11,1,1),tolower(substr($11,2,length($11))));
}
