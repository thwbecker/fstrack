#
# read in a square matrix in ASCII format and print the diagonal 
# to stdout
#
BEGIN{
  nn=0;
}
{
  if((substr($1,1,1)!="#")&&(NF>0)){
    for(i=1;i<=NF;i++)
      a[++nn] = $i;
  }
}
END{
  n = sqrt(nn);
  for(i=0;i<n;i++)
    print(a[1+i*n+i]);
}
