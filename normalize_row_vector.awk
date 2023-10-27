#
# normalize vector that is given in row form
#
BEGIN{
}
{
  if(substr($1,1,1)!="#"){
    len=0.0;
    for(i=1;i<=NF;i++)
      len += ($(i))*($(i));
    if(len >= 0.0)
      len=sqrt(len);
    else
      len=1.0;
    for(i=1;i<=NF;i++)
      printf("%.10e ",($(i))/len);
    printf("\n");
  }
}
END{

}
