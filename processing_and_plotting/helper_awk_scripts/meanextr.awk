BEGIN{
  sum=0.0;
  max=-9e20;
  min=9e20;
  n=0;
  if(col==0)
    col=1;
  printf("col %i",col);
}

{
  x[1]=$1;
  x[2]=$2;
  x[3]=$3;
  x[4]=$4;
  x[5]=$5;
  x[6]=$6;
  x[7]=$7;
  x[8]=$8;
  if(x[col]!=""){
    sum += x[col];
    if(x[col]<min)min=x[col];
    if(x[col]>max)max=x[col];
    n++;
  }
}
END{
  if(n!=0){
    sum/=n;
    printf(" n %i min %g mean %g max %g\n",n,min,sum,max);
  }else
    print("no data");
}
