#
# convert r t p to lon lat z
#
BEGIN{
  pif = 57.295779513082320876798154814105;
  
}
{
  if($1==">")
    print($0);
  else{
    if((substr($1,1,1)!="#") && (NF>=3)){
      printf("%20.15e  %20.15e  %20.15e ",
	     $3*pif,90.-pif*$2,(1-$1)*6371);
      for(i=4;i<=NF;i++)
	printf("%s ",$i);
      printf("\n");
    }
  }
}
