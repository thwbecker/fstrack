BEGIN{
  n=0;
  xmin=1e20;
  xmax=-1e20;
  xm=0.0;
  sd=0.0; 
# number of bins
  if(nb=="")
    nb=n/3;
}

{
  if($1 != "" && (substr($1,1,1)!="#")){
    n++;
    x[n]=$1;
    xm+=x[n];
    if(x[n]>xmax)xmax=x[n];
    if(x[n]<xmin)xmin=x[n];
  }
}
END{
  xm/=n;
  xstddev=0.0;
  for(i=1;i<=n;i++){
    xstddev+=(x[i]-xm)**2;
  }
  xstddev/=(n-1);
  xstddev=sqrt(xstddev);
  xr=xmax-xmin;
  dx=xr/nb;
  j=0.5;
  for(i=1;i<=nb;i++){
    bin[i]=0;
    binx[i]=xmin+j*dx;
    j+=1.0;
  }
  for(i=1;i<=n;i++){
    xl=(x[i]-xmin)/xr;
    if(xl==1)
      j=nb;
    else
      j=1+int(xl*nb);
#    print(x[i],xl,j);
    bin[j]++;
  }
  
  printf("# n: %i bins: %i xmin: %g xmax: %g xmean: %g xstd_dev: %g\n",
	 n,nb,xmin,xmax,xm,xstddev);
  print("# x_bin n_bin n_bin/n");
  check=0;
  for(i=1;i<=nb;i++){
    check += bin[i];
    printf("%12g %12g %12g\n",
	   binx[i],bin[i],bin[i]/n);
  }
  if(check!=n)
    print("error, don't add up");
}
