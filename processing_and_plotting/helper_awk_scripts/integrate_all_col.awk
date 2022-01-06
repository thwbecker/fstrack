#
#
# uses mid point method to integrate all columns
#
# cannot deal with NaN!
#
BEGIN{

  oldf=-1;
  n=0;
}
{
  if((substr($1,1,1)!="#") && (NF>1)){
    if(oldf == -1){
      oldf = NF;
      nval= oldf-1;
    }
    if(NF != oldf)print("error, change in col number") > "/dev/stderr";
    n++;
    x[n]=$1;
    for(i=1;i<=nval;i++)
      y[n*nval+i] = $(i+1);
  }
}
END{
# simple sum
  if(n>1){
    for(i=1;i<=nval;i++){
      sum = 0.0;
      for(j=2;j<=n;j++){
	sum += (x[j]-x[j-1]) * (y[j*nval+i]+y[(j-1)*nval+i])/2.;
      }
      printf("%11g ",sum);
    }
    printf("\n");
  }
    


}
