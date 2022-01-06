#
# normalize vector that is given in row form
#
BEGIN{
}
{
  if((NF>=3) && (substr($1,1,1)!="#")){
    len=$1*$1 + $2*$2 + $3*$3;
    if(len > 0.0)
      len=sqrt(len);
    else
      len=1.0;
    printf("%12g %12g %12g ",$1/len,$2/len,$3/len);
    for(i=4;i<=NF;i++)
      printf("%s ",$i);
    printf("\n");
  }
}
END{

}
