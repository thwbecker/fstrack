#
# filter out even or odd lines
#
BEGIN{
    if(even=="")
	even=-1;
    count = 0;
}
{
    if((substr($1,1,1)!="#") && ($1!="")){
	count++;
	if((even==-1)||((even==1)&&(count%2==0))||((even==0)&&(count%2!=0)))
	    print($0)
    }
}