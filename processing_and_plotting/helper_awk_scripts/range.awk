BEGIN{
  min=9e20;
  max=-9e20;
  if(col==0)
    col=1;
}

{
  if((substr($1,1,1)!="#") && $1 != ""){
    if($col!=""){
      if($col < min)min=$col;
      if($col > max)max=$col;
    }
  }
}
END{
  print(sqrt((min-max)*(min-max)));
}
