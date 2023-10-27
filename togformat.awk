BEGIN{

}
{
  if(substr($1,1,1)!="#"){
    for(i=1;i<=NF;i++)
      printf("%lg ",$i);
    printf("\n");
  }
}
