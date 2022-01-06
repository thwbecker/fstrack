#
# compute mean of orientational data with 180 deg 
# periodicity
#
# input expects azimuth in degree in column col
# if col is not set, will use the first column
# if wcol is not set, will use all data evenly, else weigh by that column
#
BEGIN{
  if(col == "")
    col=1;
  if(wcol == "")
      wcol = -1;		     # no weighting

  pif = 1.0/57.2957795130823208767; # pi/180
  ws = sx = sy = 0.0;
}

{
  if((substr($1,1,1) != "#") && (NF >= col) && (tolower($col) != "nan")){
    azi2 = $col * pif * 2;
    if(wcol == -1){
	w = 1.0;			# in case we want to weight, say by amplitude, leave unity for now
    }else{
	w = $(wcol);
    }
    sx += w * sin(azi2);
    sy += w * cos(azi2);
    ws += w;			
  }
}
END{
  if(w != 0.0){
    #
    # go back to azimuth, ie. divide by two and go to degree
    azi = atan2(sx,sy)/2/pif;
    #
    # go to 0 ... 360 convention
    #if(azi< 0)azi += 360;
    printf("%lg\n",azi);
  }else{
    print("could not find any data") > "/dev/stderr";
  }
}

