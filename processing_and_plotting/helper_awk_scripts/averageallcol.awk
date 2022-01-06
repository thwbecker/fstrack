#
# calculate the average of all columns without a "#" at the 
# beginning
#
BEGIN{

}
{
  if((substr($1,1,1)!="#")&&($1!="")){
    if(NF > n)
      n = NF;
    m = (n > NF)?(NF):(n);
    for(i=1;i<=m;i++)
      if(tolower($i) != "nan"){
	x[i] += $i;
	ic[i] ++;
      }
  }
}
END{
    for(i=1;i<=n;i++){
	if(ic[i]!=0)
	    printf("%22.16e ",x[i]/ic[i]);
	else
	    printf("NaN ");
    }
    printf("\n");
}
