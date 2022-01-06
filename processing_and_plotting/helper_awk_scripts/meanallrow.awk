#
# compute mean values for a row
#
BEGIN{
}
{
    if(substr($1,1,1)!="#"){
	for(i=1;i<=skip;i++)
	    printf("%s ",$i);
	sum=0;j=0;
	for(i=skip+1;i <= NF;i++){
	    if(tolower($i) != "nan"){
		sum += $i;
		j++;
	    }
	}
	printf("%g\n",sum/j);
    }
}
END{
}
