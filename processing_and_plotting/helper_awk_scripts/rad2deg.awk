#
# convert radians to degree
#
BEGIN{
  pif = 57.2957795130823;
}
{
  if((NF>=1)&&(substr($1,1,1)!="#")){
    for(i=1;i<=NF;i++)
      printf("%g ",$i*pif);
    printf("\n");
  }else{
    print($0);
  }
}
