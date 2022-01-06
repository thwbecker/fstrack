#
# calculate the sqrt((\sum_i x_i^2)/n) RMS
# there's also vrms which takes out the mean first
# 
BEGIN{
  n=0;
  if(col==0)
    col=1;
  rms=0;
}
{
  if((NF>=col) && 
     (substr($1,1,1)!="#")&&
     (tolower($col)!="nan")&&
     ($col!="")){
      n++;
      rms += ($(col)) * ($(col));
  }
}
END{
  if(n){
      printf("%.15e\n",sqrt(rms/n));
  }else{
    print("nan");
  }
}
