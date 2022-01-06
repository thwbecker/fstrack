BEGIN{
  n=0;
}
{
  n++;
  x[n]=$1;y[n]=$2;
}
END{
  if(gm){
    meanx=1.0;meany=1.0;
    for(i=1;i<=n;i++){
      meanx*=x[i];
      meany*=y[i];
    }
    meanx=(meanx^(1./n));
    meany=(meany^(1./n));
  }else{ 
    meanx=0.0;meany=0.0;
    for(i=1;i<=n;i++){
      meanx+=x[i];
      meany+=y[i];
    }
    meanx/=n;
    meany/=n;
  }
  sdx=0.0;sdy=0.0;
  for(i=1;i<=n;i++){
    sdx += (x[i]-meanx)*(x[i]-meanx);
    sdy += (y[i]-meany)*(y[i]-meany);
  }
  sd =sdx+sdy;
  sd/=2*n-2;
  sd*=(2/n);
  sd=sqrt(sd);
  t=(meanx-meany)/sd;
  print(t);
}
