BEGIN{

  n=0;
}
{
  n++;
  x[n]=$1;y[n]=$2;
}
END{
  meanx=meany=0.0;
  for(i=1;i<=n;i++){
    meanx+=x[i];
    meany+=y[i];
  }
  meanx/=n;
  meany/=n;
  sum1=sum2=sum3=0.0;
  for(i=1;i<=n;i++){
    xm=x[i]-meanx;
    ym=y[i]-meany;
    sum1+=xm*ym;
    sum2+=xm*xm;
    sum3+=ym*ym;
  }
  r= sum1/sqrt(sum2*sum3);
  nu=n-2;
  t=r*sqrt(nu/(1-r*r));
  printf("r: %g t: %g n: %g nu: %g\n",
	 r,t,n,nu);
}
