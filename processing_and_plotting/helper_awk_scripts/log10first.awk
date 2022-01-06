# print decadic log of first column, and leave rest
BEGIN{
#    f=1.0/ln(10.0);
    f=0.4342944819032518;
}
{
  if((substr($1,1,1)!="#")&&(NF>=1)){
      printf("%20.15e ",log($1)*f);
      for(i=2;i<=NF;i++)
	  printf("%s ",$i);
      printf("\n");
  }
}
