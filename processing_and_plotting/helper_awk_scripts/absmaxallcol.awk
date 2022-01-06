BEGIN{
  max=-9e20;
  if(col==0)
    col=1;
}

{
  if((substr($1,1,1)!="#")&&($1!="")){
    for(i=1;i<=NF;i++){
      if(tolower($i)!="nan"){
	c=($i>=0)?($i):(-$i);
	if(c > max)max=c;
      }
    }
  }
}
END{
  print(max);
}
