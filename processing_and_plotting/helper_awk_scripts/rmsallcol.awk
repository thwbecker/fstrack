#
# calculate the sqrt((\sum_i x_i^2)/n) RMS
#
BEGIN{
  mc=0;
}
{
  if((NF>=1)&&(substr($1,1,1)!="#")){
    if(NF>mc){
      mc = NF;
    }
    for(i=1;i<=mc;i++){
      if(tolower($i) != "nan"){
	n[i]++;
	x[i] += $(i) * $(i);
      }
    }
  }
}
END{
  for(i=1;i<=mc;i++){
    if(n[i])
      printf("%lg ",sqrt(x[i]/n[i]));
    else
      printf("nan");
  }
  printf("\n");

}
