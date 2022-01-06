#
# remove a grid mean integrated over y
#
BEGIN{
  i=1;
  j=1;
}
{
  z[j*nx+i] = $3;
  if(j == 1)
    x[i] = $1;
  if(i == 1)
    y[j] = $2;
  i++;
  if(i > nx){
    i = 1;
    j++;
  }
}
END{
  if(mode == 1){
    # compute mean with regard to y 
    if(mean){
      for(j=1;j<=nx;j++) {
	sum[j] = 0.0;
	for(i=1;i<=ny;i++)
	  sum[j] += z[i*nx+j];
	sum[j] /= ny;
      }
    }else{
      for(i=1;i<=nx;i++)
	sum[i] = z[1*nx+i];

    }

    for(i=1;i<=ny;i++)
      for(j=1;j<=nx;j++)
	printf("%g %g %g\n",x[j],y[i],z[i*nx+j]-sum[j]);
  }else{
    # compute mean with regard to x
    for(j=1;j<=ny;j++)
      sum[j] = 0.0;
    
    for(j=1;j<=nx;j++)
      for(i=1;i<=ny;i++)
	sum[i] += z[i*nx+j];
    
    for(j=1;j<=ny;j++)
      sum[j] /= nx;

    for(i=1;i<=ny;i++)
      for(j=1;j<=nx;j++)
	printf("%g %g %g\n",x[j],y[i],z[i*nx+j]-sum[i]);




  }
}