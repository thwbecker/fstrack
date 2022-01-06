#
# print decadic log of input
# 
#
BEGIN{
#    f=1.0/ln(10.0);
    f=0.4342944819032518;
}
{
  if((substr($1,1,1)!="#")&&(NF>=1)){
      if(only_first){		# only take log of first
	  printf("%20.15e ",log($1)*f);
	  for(i=2;i<=NF;i++)
	      printf("%s ",$i);
      }else{			# all
	  for(i=1;i<=NF;i++)
	      printf("%20.15e ",log($i)*f);

      }
    printf("\n");
  }
}
