#
# prints the norm of all columns
#
BEGIN{
  
}
{
  if(substr($1,1,1)!="#"){
    for(i=1;i<=NF;i++){
      p[i]++;
      sum[i] += ($i)*($i);
    }
    if(NF>maxnf)maxnf=NF;
  }
}
END{
  for(i=1;i<=maxnf;i++)
    printf("%g ",sum[i]/p[i]);
  printf("\n");
}
