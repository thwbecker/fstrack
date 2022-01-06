#
# calculate the standard deviation of all columns, 
# fast and inaccurate
#
BEGIN {

}
{
  if(substr($1,1,1)!="#"){
    if(NF>cmax)cmax=NF;
    for(i=1;i<=NF;i++){
      if(tolower($i)!="nan"){
	mean[i] += $i;
	sum2[i] += $i * $i;
	n[i]++;
      }
    }
  }
}
END {
    for(i=1;i<=cmax;i++){
      if(n[i] > 1){
	std = sqrt ((n[i] * sum2[i] - mean[i] * mean[i]) / ((n[i]*(n[i]-1))));
      }else{
	std=0;
      }
      printf("%g ",std);
    }
    printf("\n");
    
}
