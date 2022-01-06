# convert unix timestamp seconds to decimal year
#
BEGIN{


}
{
    unixt=$1;

    year = strftime("%Y",unixt);
    tday_per_year = (year%4==0)?(366):(365);
    day_of_year = (strftime("%j",unixt)); # 1...365 or 366
    hour = strftime("%H",unixt);	  # 0...23
    min = strftime("%M",unixt);		  # 0...59
    sec = strftime("%S",unixt);		  # 0..60
    
    yfrac = ((day_of_year-1)+(hour+min/60+sec/360)/24)/tday_per_year;
    
    
    printf("%.12f\n",year+yfrac);

}
