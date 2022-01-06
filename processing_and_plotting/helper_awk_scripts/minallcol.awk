BEGIN{
}
{
  if((substr($1,1,1)!="#")&& ($1 != ">")){
    for(i=1;i<=NF;i++){
      if(NF>minnf)minnf=NF;
      if(tolower($i) != "nan"){
	c[i]++;
	if(c[i]==1)
	  min[i]=$i;
	else if($i < min[i])
	  min[i]=$i;
      }
    }
  }
}
END{
  for(i=1;i<=minnf;i++)
    printf("%s ",min[i]);
  printf("\n");
}
