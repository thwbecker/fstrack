# find max in all rows
BEGIN{
}
{
  if((substr($1,1,1)!="#")&& ($1 != ">")){
    for(i=1;i<=NF;i++){
      if(NF>maxnf)maxnf=NF;
      if(tolower($i) != "nan"){
	c[i]++;
	if(c[i]==1)
	  max[i]=$i;
	else if($i > max[i])
	  max[i]=$i;
      }
    }
  }
}
END{
  for(i=1;i<=maxnf;i++)
    printf("%s ",max[i]);
  printf("\n");
}
