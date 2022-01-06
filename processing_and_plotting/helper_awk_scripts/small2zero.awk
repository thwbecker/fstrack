BEGIN{

}
{
  if(substr($1,1,1)!="#"){
    for(i=1;i<=NF;i++){
      x=$i;if(x<0)x=-x;
      if(x < 1e-14)
	printf("0 ");
      else
	printf("%lg ",$i);
    }
    printf("\n");
  }
}
