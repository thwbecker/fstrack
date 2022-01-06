BEGIN{
}
{
  if((substr($1,1,1)!="#") && $1!=""){
    dim=NF;
    for(i=1;i<=NF;i++){
      s[i]+=$i;
    }
    n++;
  }
}
END{
  for(i=1;i<=dim;i++)
    printf("%g ",s[i]);
  printf("\n");
}
