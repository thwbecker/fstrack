#
# compute the symmetric part of a 3x3 tensor and print upper right hand side 
#
BEGIN{
}
{
  if(substr($1,1,1)!="#"){
    if(NF >= 9){
      k=1;
      for(i=1;i<=3;i++)
	for(j=1;j<=3;j++){
	  a[(i-1)*3+j] = $k;
	  k++;
	}
      for(i=1;i<=3;i++)
	for(j=i;j<=3;j++)
	  printf("%8.4e ",.5*(a[(i-1)*3+j]+a[(j-1)*3+i]));
      for(i=10;i<=NF;i++)
	printf("%s ",$i);
      printf("\n");
    }
  }
}
