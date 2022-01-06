BEGIN{
  max=-9e20;
  if(col==0)
    col=1;
}

{
  if((substr($1,1,1)!="#") && (NF>=col) && (tolower($col)!="nan")){
    c=($col>=0)?($col):(-$col);
    if(c > max)max=c;
  }
}
END{
  print(max);
}
