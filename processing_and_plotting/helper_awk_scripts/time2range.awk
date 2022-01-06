#
# given input time in 
#
# year month day hours minutes dt_left dt_right
#
# format, print out 
#
# year1 month1 day1 hours1 minutes1 year2 month2 day2 hours2 minutes2
#
# where the first and second time are dt_left minutes and dt_right minutes away from the input
# date and time
#
BEGIN{

}
{
    datespec=sprintf("%04i %02i %02i %02i %02i 00",$1,$2,$3,$4,$5);
    time_stamp0 = mktime(datespec); # time in seconds since epoch
    time_stamp1 = time_stamp0 + $6*60;
    time_stamp2 = time_stamp0 + $7*60;

    print(strftime("%Y %m %d %H %M",time_stamp1),
	  strftime("%Y %m %d %H %M",time_stamp2))

}
