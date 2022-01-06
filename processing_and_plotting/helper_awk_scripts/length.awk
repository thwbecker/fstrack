BEGIN{

}
{
if(($1 !="") && (substr($1,1,1) !="#")){
  x+=$1*$1;
}
}
END{
  if(x<=0.0)
    print(0);
  else
    print(sqrt(x))
} 
