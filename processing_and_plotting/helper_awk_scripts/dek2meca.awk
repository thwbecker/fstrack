#
# awk script to read Harvard CMT moment tensor solutions
# in the "dek" format
#
# is select is set to
# 0: (default) use all events
# 1: select thrust mechanisms only
# 2: select strike-slip mechanisms only
# 3: select normal mechanism only
# 4: select eventid `event'
#
#
# format selects the output:
#
# 0: CMT best double couple 
# 1: moment tensor 
# 2: lon lat z moment[Nm] 
# 3: e_plunge[1] e_plunge[2] 
# 4: lon lat z mw year time_sec_since_epoch
#
# $Id: dek2meca.awk,v 1.5 2010/07/06 22:11:31 becker Exp $
#
BEGIN{
  ls=1;
}
{
# we are in line number one of the four data lines blocks
  if(NR==ls){
#
#B010177C  1/ 1/77 11:33:41.6  30.66  137.06 476.05.20.0SOUTH OF HONSHU, JAPAN
#^        ^        ^          ^     ^       ^     ^     ^
#1        10       19         30    36      44    50    56
# the events name
    eventid=substr($0,1,8);
    date=substr($0,10,8);split(date,dstring,"/");
    month=dstring[1];day=dstring[2];year=dstring[3];
    if(year<70)year += 2000;else year+=1900;
    time=substr($0,19,10);split(time,tstring,":");
    hour=tstring[1];minute=tstring[2];second=tstring[3];
# convert to seconds since epoch
    the_time=sprintf("%i %i %i %i %i %i",year,month,day,hour,minute,second);
    secs = mktime(the_time);

    lat=sprintf("%lf",substr($0,30,6));
    lon=sprintf("%lf",substr($0,37,7));

    if(lon<0.0)lon+=360.0;
    dep=sprintf("%lf",substr($0,45,6));
    mb=sprintf("%lf",substr($0,50,3)); # Mb
  }else if(NR==ls+1){# in line number two
#MLI BW: 5 14  45 MW: 0  0   0 DT=   4.3 0.7  30.62 0.07  136.80 0.10 476.5  4.8
#^                             ^             ^           ^            ^     ^   
#1	                       31            45          57           70    76
    if(substr($0,31,3)!="DT="){
      print("error line two formatting",substr($0,31,3));
      nextline;
    }
    latc=sprintf("%lf",substr($0,45,6));# centroid latitude
    lonc=sprintf("%lf",substr($0,57,7));# centroid longitude
    depc=sprintf("%lf",substr($0,70,5));# centroid depth
    if(lonc<0)lonc+=360.0;
  }else if(NR==ls+2){
#
# line number three
#
# moment tensor
    if(substr($0,10,2)!="EX"){
      printf("line mismatch, 10 and 11 are not EX\n") > "/dev/stderr";
    }
    exponent=sprintf("%i",substr($0,13,2));
# Mrr, Mss, Mee, Mrs, Mre, Mse or
# mrr, mtt, mpp, mrt, mrp, mtp
# DUR 1.0 EX 23  4.25 1.05  5.74 1.04-10.00 1.89  4.60 1.18  1.15 1.09 -9.28 1.01
#              ^--> 15th field
#                    ^-->21th field
#                         ^-->26th field
#                               ^-->32nd field 
#              1234561234512345612345
    sf=15;
    for(i=1;i<=6;i++){
      m[i]=  sprintf("%f",substr($0,sf,6));sf+=6;
      msd[i]=sprintf("%f",substr($0,sf,5));sf+=5;
    }
  }else{
# fourth line
# best double couple
    if(NF != 16){
      print("fourth line doesn't have 16 fields") > "/dev/stderr";
      nextfile;
    }
    for(i=1;i<=3;i++){# eigenvectors
      e_val[i]=$(1+(i-1)*3);
      e_plunge[i]=$(2+(i-1)*3);
      e_strike[i]=$(3+(i-1)*3);
    }
    # best double couple
    scalar_moment=$10;# in units of 10**24
    for(i=1;i<=2;i++){
      strike[i]=$(11+(i-1)*3);# first and second nodal planes
      dip[i]=$(12+(i-1)*3);
      rake[i]=$(13+(i-1)*3);
    }

    mw = 2/3*(log(scalar_moment*10**(exponent))*0.4342944819032518)-10.7);

#
#
# OUTPUT OF EVENT 
#
# do we want this event?
    if(select > 0){
      use=0;
      if(select==1){
	if((e_plunge[1]>=45)&&(e_plunge[2]<45))# thrust
	  use=1;
      }else if(select==2){# strike slip, counts also those with 
	if(e_plunge[2]>=45) # both >= 45
	  use=1;
      }else if(select==3){# normal
	if((e_plunge[1]<45)&&(e_plunge[2]<45))
	  use=1;
      }else if(select==4){
	if(event==""){
	  print("error, for select mode 4, have to define variable event") > \
	    "/dev/stderr";
	  nextfile;
	}
	if((eventid==event)||(substr(eventid,2)==event))
	  use=1
      }else{
	print("error, select mode ",select," undefined") > "/dev/stderr";
	nextfile;
      }
    }else{
      use=1;
    }
    if(use){
      if(format==0)
#
# use centroid location 
#
#
# (c) Focal mechanisms in CMT convention for best double couple
#             lonc, latc, strike1, dip1, rake1, strike2, dip2, rake2, moment, newX, newY, event_title
#             with moment in 2 columns : mantiss and exponent corresponding to seismic moment in dynes-cm
	print(lonc,latc,depc,strike[1],dip[1], rake[1], 
	      strike[2], dip[2], rake[2], scalar_moment, exponent, lonc, latc,
	      eventid);
      else if(format==1)
# (m) Moment tensor
#	lonc, latc, depth, mrr, mtt, mff, mrt, mrf, mtf, exp, newX, newY, event_title
	printf("%8.3f %8.3f %8.2f %8.5f %8.5f %8.5f %8.5f %8.5f %8.5f %2i %8.3f %8.3f %s\n",
	      lonc,latc,depc,m[1],m[2],m[3],m[4],m[5],m[6],exponent,lonc,latc,
	       eventid);
#
# lonc latc depth moment [Nm, 1 dyne cm = 10^{-7} Nm] time lon lat
#
      else if(format==2){
	printf("%8.3f %8.3f %8.3f %12.4e %02i/%02i/%04i:%02i:%02i:%02i %g %g\n",
	       lonc,latc,depc,scalar_moment*10**(exponent-7),
	       month,day,year,hour,minute,second,lon,lat);
      }
# e_plunge[1] e_plunge[2]
      else if(format==3)
	print(e_plunge[1],e_plunge[2]);
      else if(format==4)
# 4: lon lat z mw year time_sec_since_epoch 
	printf("%8.2f %8.2f %6.1f %5.1f %4i %12i\n",
	       lonc,latc,depc,mw,year,secs);

      else {
	print("error, format",format,"undefined") > "/dev/stderr";
	nextfile;
      }
    }
#
#reset the current data block line counter
#
    ls=NR+1;
  }
  
}
END{
  
}
