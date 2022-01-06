#
# convert a direction [0...360] to orientation of a 180deg period quantity
#
#
BEGIN{
}
{
    if(substr($1,1,1)!="#"){
	if((tolower($1)=="nan")||(tolower($1)=="inf")||(tolower($1)=="-inf")){
	    printf("%s ",$1);
	}else{
	    azi=$1;
	    while(azi < 0)
		azi += 360;
	    while(azi >= 360)
		azi -= 360;
	    if(azi > 180)
		azi = azi-180;
	    printf("%g ",azi);
	}

	for(i=2;i<=NF;i++)
	    printf("%s ",$i);
	printf("\n");
    }
}
