# convert year month day string to unix second
BEGIN{
}
{
    if(substr($1,1,1)!="#" && NF>=3){
	datespec=sprintf("%04i %02i %02i 0 0 0",$1,$2,$3);
	time_stamp = mktime(datespec);
	printf("%i ",time_stamp);
	for(i=4;i<=NF;i++)
	    printf("%s ",$i);
	printf("\n");
    }
}
