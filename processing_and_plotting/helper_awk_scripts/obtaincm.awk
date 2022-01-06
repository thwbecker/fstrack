BEGIN{
  xsum[1]=0.0;
  xsum[2]=0.0;
  xsum[3]=0.0;
  n=0;
}
{
  if((substr($1,1,1)!="#") && $1 !=""){
    x[1]=$1;
    x[2]=$2;
    x[3]=$3;
    dx[1]=$4;
    dx[2]=$5;
    dx[3]=$6;
    
    for(i=1;i<=3;i++)
      xsum[i] += dx[i];
    n++;
  }
}
END{
  print(xsum[1]/n,xsum[2]/n,xsum[3]/n);
}
