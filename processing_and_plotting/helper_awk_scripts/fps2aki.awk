# convert FPFIT format haukkson catalog to Aki & Richards convention
#   take out desired info, adjust format of date, lat & lon.
#   take out flagged multiple solutions, since by inspection these have the higher misfit 
#
# throws out all with misfit > 0.2 and station dist < 0.5
#
BEGIN{
    misfit_crit = 0.2;
    stdr_crit = 0.5;
}
{
#get time info.
  year=substr($0,1,4);
  month=substr($0,5,2);
  day=substr($0,7,2);
  hour=substr($0,10,2);
  min=substr($0,12,2);
  sec=substr($0,15,2);
# printf("%4i %2i %2i %2i %2i %2i\n",year,month,day,hour,min,sec) > "/dev/stderr";   
	
#convert time to years after 1970.
  time=sprintf("%4i %2i %2i %2i %2i %2i",year,month,day,hour,min,sec);
  tsec=sprintf("%.12e",(mktime(time)-28800)/31557600);    #NOTE: correct from UTC, convert to proportion of a 365.25 day year.

#get location info., reformat
  lon=substr($0,29,4);
  east=tolower(substr($0,33,1));
  lon2=substr($0,34,5);
  long=(lon+lon2/60.)*((east=="e") ? 1 : -1);
  if(long < 0)long += 360.0;
  lat=substr($0,20,3);
  lat2=substr($0,24,5);
  south=tolower(substr($0,23,1));
  lati=(lat+lat2/60.)*((south=="s") ? -1 : 1);
           
  depth=substr($0,39,7);

#get parameters required for potency tensor calculation.
  mag=substr($0,48,6);
  strike=substr($0,84,3) - 90;
  if(strike<0)strike+=360.;
  dip=substr($0,88,2);
  rake=substr($0,90,4);

  # misfit 0..1
  misfit = substr($0,96,4);
  
  # number of observations
  nobs=substr($0,101,3);

  # azimuthal gap
  agap=substr($0,56,4)
  
  # station distribution ratio
  stdr=substr($0,111,4);

#don't include non-converged solutions
  conv=substr($0,130,1);
#don't include multiple solutions.  Take the one with the lowest misfit.           
  repeat=substr($0,131,1);
  if((repeat != "*")&&(conv != "C")&&(misfit <= misfit_crit)&&(stdr>=stdr_crit))
    {
	#printf("%7.4f %9.4f %8.4f %6.2f\t%4.0f %4.0f %4.0f\t%4.2f \n",tsec,long,lati,depth,strike,dip,rake,mag);

#	print(stdr,stdr_crit) > "/dev/stderr"
	printf("%9.4f %8.4f %6.2f\t%4.0f %4.0f %4.0f\t%4.2f\t%11g %11g\t%20s\n",\
	       long,lati,depth,strike,dip,rake,mag,long,lati,tsec);
#	print(mag,misfit,nobs)
    }

#  print(nobs);
}
END{}
