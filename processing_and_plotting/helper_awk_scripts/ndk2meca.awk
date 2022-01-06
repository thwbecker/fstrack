#
# gawk script to read new Global (Harvard) CMT moment tensor solutions
# in the new "ndk" format
#
# set parmaters like 
#
#     gawk --assign select=0 --assign format=0 ndk2meca.awk jan76_dec05.ndk
#
#
# if select flag is set to
#
# 0: (default) will print all events
#
# 1: will select thrust mechanisms only
# 2: will select strike-slip mechanisms only
# 3: will select normal mechanism only
# 4: select eventid `event'
#
#
# format selects the output:
#
# 0: CMT best double couple in psmeca -Sc format (default) 
# 1: moment tensor          in psmeca -Sm
# 2: lon lat z moment[Nm]  date time 
# 3: e_plunge[1] e_plunge[2] 
# 4: lon lat z mw year time_sec_since_epoch
# 5:  moment tensor          in psmeca -Sm with time since epoch instead of event ID

# 6: date         time       elat    elon   edep  Mag
# 7:  moment tensor          in psmeca -Sm with Mw-date-time instead of event ID
# $Id: ndk2meca.awk,v 1.4 2010/07/06 22:06:57 becker Exp becker $
#
BEGIN{
  ls=1;
  lc=0;
}
{
  if(($1!="") && (substr($1,1,1)!="#")){
    # line counter, equal to NR if no empty lines are found in file
    lc++;

    # we are in line number one of the five data lines blocks
    if(lc == ls){
      #
      # the events name
      #First line: Hypocenter line
      #[1-4]   Hypocenter reference catalog (e.g., PDE for USGS location, ISC for
      #        ISC catalog, SWE for surface-wave location, [Ekstrom, BSSA, 2006])
      #[6-15]  Date of reference event
      #[17-26] Time of reference event
      #[28-33] Latitude
      #[35-41] Longitude
      #[43-47] Depth
      #[49-55] Reported magnitudes, usually mb and MS
      #[57-80] Geographical location (24 characters)


      catalog=substr($0,1,4);

      date=substr($0,6,10);split(date,dstring,"/");
      month=dstring[2];day=dstring[3];year=dstring[1];

      time=substr($0,17,10);split(time,tstring,":");
      hour=tstring[1];minute=tstring[2];second=tstring[3];

      # convert to seconds since epoch
      the_time=sprintf("%i %i %i %i %i %i",year,month,day,hour,minute,int(second+0.5));
      secs = mktime(the_time);


      lat=sprintf("%lf",substr($0,28,6));

      lon=sprintf("%lf",substr($0,35,7));

      if(lon<0.0)
	lon+=360.0;
    
      dep=sprintf("%lf",substr($0,43,5));

      mb=sprintf("%lf",substr($0,49,3)); # Mb
      
      # location string
      loc_string=substr($0,57,24);

      #print(catalog,month,day,year,hour,minute,second,lat,lon,dep,mb,loc_string);
      
    }else if(lc == ls+1){
#
# in line number two
# Second line: CMT info (1)
#
#
#[1-16]  CMT event name. This string is a unique CMT-event identifier. Older
#        events have 8-character names, current ones have 14-character names.
#        See note (1) below for the naming conventions used.
#[18-61] Data used in the CMT inversion. Three data types may be used: 
#        Long-period body waves (B), Intermediate-period surface waves (S),
#        and long-period mantle waves (M). For each data type, three values
#        are given: the number of stations used, the number of components 
#        used, and the shortest period used.
#[63-68] Type of source inverted for: "CMT: 0" - general moment tensor; 
#        "CMT: 1" - moment tensor with constraint of zero trace (standard); 
#        "CMT: 2" - double-couple source.
#[70-80] Type and duration of moment-rate function assumed in the inversion. 
#        "TRIHD" indicates a triangular moment-rate function, "BOXHD" indicates
#        a boxcar moment-rate function. The value given is half the duration
#        of the moment-rate function. This value is assumed in the inversion,
#        following a standard scaling relationship (see note (2) below),
#        and is not derived from the analysis.

      eventid=substr($0,1,17);
#      print(eventid)
    }else if(lc == ls+2){
      #
      # line number three
      #

#Third line: CMT info (2)
#[1-58]  Centroid parameters determined in the inversion. Centroid time, given
#        with respect to the reference time, centroid latitude, centroid
#        longitude, and centroid depth. The value of each variable is followed
#        by its estimated standard error. See note (3) below for cases in
#        which the hypocentral coordinates are held fixed.
#[60-63] Type of depth. "FREE" indicates that the depth was a result of the
#        inversion; "FIX " that the depth was fixed and not inverted for;
#        "BDY " that the depth was fixed based on modeling of broad-band 
#        P waveforms.
#[65-80] Timestamp. This 16-character string identifies the type of analysis that
#        led to the given CMT results and, for recent events, the date and 
#        time of the analysis. This is useful to distinguish Quick CMTs ("Q-"), 
#        calculated within hours of an event, from Standard CMTs ("S-"), which 
#        are calculated later. The format for this string should not be 
#        considered fixed.

#  1             2   3    4     5     6      7     8     9   10 
#CENTROID:     13.8 0.2 -29.25 0.02 -176.96 0.01  47.8  0.6 FREE O-00000000000000

      latc = $4;
      lonc = $6;
      if(lonc < 0)lonc += 360;
      depc=$8;
 
#     print(latc,lonc,depc);

    }else if(lc == ls+3){
      #
      # fourth line
      #
#Fourth line: CMT info (3)
#[1-2]   The exponent for all following moment values. For example, if the
#        exponent is given as 24, the moment values that follow, expressed in 
#        dyne-cm, should be multiplied by 10**24.
#[3-80]  The six moment-tensor elements: Mrr, Mtt, Mpp, Mrt, Mrp, Mtp, where r
#        is up, t is south, and p is east. See Aki and Richards for conversions
#        to other coordinate systems. The value of each moment-tensor
#	  element is followed by its estimated standard error. See note (4)
#	  below for cases in which some elements are constrained in the inversion.
    
      exponent=$1;
      for(i=1;i<=6;i++){
	m[i]=  $(2+(i-1)*2)	# don't multiply with 10**exponent here
	msd[i]=$(3+(i-1)*2);
      }
#      print(exponent,m[1],m[2],m[3],m[4],m[5],m[6]);
      
    }else if(lc == ls+4){
      # 
      # fifth line
      # 
#Fifth line: CMT info (4)
#[1-3]   Version code. This three-character string is used to track the version 
#        of the program that generates the "ndk" file.
#[4-48]  Moment tensor expressed in its principal-axis system: eigenvalue, 
#        plunge, and azimuth of the three eigenvectors. The eigenvalue should be
#        multiplied by 10**(exponent) as given on line four.
#[50-56] Scalar moment, to be multiplied by 10**(exponent) as given on line four.
#[58-80] Strike, dip, and rake for first nodal plane of the best-double-couple 
#        mechanism, repeated for the second nodal plane. The angles are defined
#        as in Aki and Richards.
#
# 1       2    3  4      5    6  7    8     9   10     11   12 13   14  15 16   17
#
# V10   8.940 75 283   1.260  2  19 -10.190 15 110   9.560 202 30   93  18 60   88
      for(i=1;i <= 3;i++){# eigenvectors
	e_val[i]=   $(2+(i-1)*3); # don't multiply with 10**exponent
	e_plunge[i]=$(3+(i-1)*3);
	e_strike[i]=$(4+(i-1)*3);
      }
      # best double couple
      scalar_moment=$11;# in units of 10**24
      for(i=1;i <= 2;i++){
	strike[i]=$(12+(i-1)*3);# first and second nodal planes
	dip[i]=$(13+(i-1)*3);
	rake[i]=$(14+(i-1)*3);
      }
#      print(e_val[1],e_plunge[1],e_strike[1],e_val[2],e_plunge[2],e_strike[2],e_val[3],e_plunge[3],e_strike[3]);
#      print(scalar_moment,strike[1],dip[1],rake[1],strike[2],dip[2],rake[2]);


      mw = 2/3*(log(scalar_moment*10**(exponent))*0.4342944819032518)-10.7;

      #
      #
      # OUTPUT OF EVENT 
      #
      # do we want this event?
      if(select > 0){
	use=0;
	if(select==1){
	  if((e_plunge[1] >= 45)&&(e_plunge[2] < 45))# thrust
	    use=1;
	}else if(select == 2){# strike slip, counts also those with 
	  if(e_plunge[2] >= 45) # both >= 45
	    use=1;
	}else if(select == 3){# normal
	  if((e_plunge[1]<45)&&(e_plunge[2]<45))
	    use=1;
	}else if(select==4){
	  if(event==""){
	  print("error, for select mode 4, have to define variable event") > \
	    "/dev/stderr";
	  nextfile;
	  }
	  if((eventid==event) || (substr(eventid,2)==event))
	    use=1
	}else{
	  print("error, select mode ",select," undefined") > "/dev/stderr";
	  nextfile;
	}
      }else{
	use=1;
      }
      if(use){
	if(format==0){
	  #
	  # use centroid location 
	  #
	  #
	  # (c) Focal mechanisms in CMT convention for best double couple
	  #             lonc, latc, strike1, dip1, rake1, strike2, dip2, rake2, moment, newX, newY, event_title
	  #
	  #             with moment in 2 columns : mantiss and exponent corresponding to seismic moment in dynes-cm
	  #
	  print(lonc,latc,depc,
		strike[1],dip[1], rake[1], strike[2], dip[2], rake[2], 
		scalar_moment, exponent, lonc, latc,
		eventid);
	}else if(format==1){
	  #
	  #
	  # (m) Moment tensor, components are not scaled up with 10**exponent
	  #
	  #	lonc, latc, depth, mrr, mtt, mff, mrt, mrf, mtf, exp, newX, newY, event_title
	  #
	  printf("%8.3f %8.3f %8.2f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %2i %8.3f %8.3f %s\n",
		 lonc,latc,depc,m[1],m[2],m[3],m[4],m[5],m[6],exponent,lonc,latc,eventid);
	#
	# lonc latc depth moment [Nm, 1 dyne cm = 10^{-5} N 10^{-2} m = 10^{-7} Nm] time lon lat
	#
	}else if(format==2){
	  printf("%8.3f %8.3f %8.3f %12.4e %02i/%02i/%04i:%02i:%02i:%02i %g %g\n",
		 lonc,latc,depc,scalar_moment*10**(exponent-7),
		 month,day,year,hour,minute,second,lon,lat);
	}else if(format==3){
	  #
	  # two plunges
	  #
	  print(e_plunge[1],e_plunge[2]);
	}else if(format==4){
	  #
	  # 4: lon lat z mw year time_sec_since_epoch
	  #
	  printf("%8.2f %8.2f %6.1f %5.1f %4i %12i\n",
		 lonc,latc,depc,mw,year,secs);
	}else if(format == 5){
	    #
	    # 5: (m) type moment tensor with time in sec since epoch
	    # lonc, latc, depth, mrr, mtt, mff, mrt, mrf, mtf, exp, newX, newY, time
	    #
	    printf("%8.3f %8.3f %8.2f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %2i %8.3f %8.3f %s\n",\
		 lonc,latc,depc,m[1],m[2],m[3],m[4],m[5],m[6],exponent,lonc,latc,
		  secs);
	}else if(format==6){
#                                       date         time       elat    elon   edep  Mag

	    printf("%4i-%02i-%02i %02i:%02i:%07.4f %8.3f %8.3f %8.1f %4.1f\n",
		   year,month,day,hour,minute,second,latc,lonc,depc,mw);


	}else if(format == 7){
	    #
	    # 5: (m) type moment tensor with time in sec since epoch
	    # lonc, latc, depth, mrr, mtt, mff, mrt, mrf, mtf, exp, newX, newY, Mw-date-time
	    #
	    printf("%8.3f %8.3f %8.2f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %2i %8.3f %8.3f %.1f-%s-%s-%s-%s-%s\n",\
		   lonc,latc,depc,m[1],m[2],m[3],m[4],m[5],m[6],exponent,lonc,latc,
		   mw,year,month,day,hour,minute);
	}else {
	  print("error, format",format,"undefined") > "/dev/stderr";
	  nextfile;
	}
      } # end use branch

      
      #
      # reset the current data block line counter
      #
      ls = lc+1;

    } # end fifth line branch

  } # non-zero line
}
END{
  
}
