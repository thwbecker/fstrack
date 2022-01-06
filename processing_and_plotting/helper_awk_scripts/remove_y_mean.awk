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
  # compute mean with regard to y 
  for(i=1;i<=nx;i++)
    sum[i] = 0.0;
  
  for(i=1;i<=ny;i++)
    for(j=1;j<=nx;j++)
      sum[i] += z[i*nx+j];
  for(i=1;i<=nx;i++)
    sum[i] /= ny;
  for(i=1;i<=ny;i++)
    for(j=1;j<=nx;j++)
      printf("%g %g %g\n",x[j],y[i],z[i*nx+j]-sum[j]);
}