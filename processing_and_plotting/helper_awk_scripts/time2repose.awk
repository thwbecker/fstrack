BEGIN{
  b=0;
}
{
  if($1!="" && (substr($1,1,1)!="#")){
    n++;
    if(n==1){
      xo=$1;
    }else{
      x=$1;
      print(x-xo);
      xo=x;
    }
  }
  
}
