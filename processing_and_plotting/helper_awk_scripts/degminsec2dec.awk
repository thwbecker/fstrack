# lon lat
# 39°40'32" 6°1'42"
BEGIN{

}
{
    if(NF>=2){
	lon=$1;lat=$2;

	split(lon,lona,"°");split(lona[2],lonaa,"'");split(lonaa[2],lonaaa,"\"");
	split(lat,lata,"°");split(lata[2],lataa,"'");split(lataa[2],lataaa,"\"");

	lon = lona[1]+lonaa[1]/60+lonaaa[1]/60/60;
	lat = lata[1]+lataa[1]/60+lataaa[1]/60/60;

	printf("%20.10e %20.10e ",lon,lat);
 
	for(i=3;i<=NF;i++)
	    printf("%s ",$i);
	printf("\n");
    }

}