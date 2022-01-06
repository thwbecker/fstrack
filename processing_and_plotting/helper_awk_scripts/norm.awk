#
# prints the norm of the first column
#
BEGIN{
  sum=0.0;
  if(col==0)
    col=1;
}
{
  if($1!="" && (substr($1,1,1)!="#")){
    sum += ($col)*($col);
  }
}
END{
  print(sqrt(sum));
}
