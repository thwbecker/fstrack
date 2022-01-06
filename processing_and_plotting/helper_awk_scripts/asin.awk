BEGIN{


}
{
  if((NF>=1)&&(substr($1,1,1)!="#")){
    for(i=1;i<=NF;i++)
      printf("%20.15e ",asin($i));
    printf("\n");
  }else{
    print($0);
  }

}

function asin( x ) {
  tmp=atan2(x,sqrt(1.0-x*x+1.0e-17));
  return(tmp);
}
