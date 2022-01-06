{
if(substr($1,1,1)!="#"){
  for(i=1;i<=NF;i++){
    x=$i;
    if(x<0)
      x= -x;
    printf("%.15e ",x);
  }
  printf("\n");
}

}
