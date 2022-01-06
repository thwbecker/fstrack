#
# convert a 6x6 elasticity matrix and scale it if fac is set
#
BEGIN{
  if(fac=="")
    fac = 1.0;
}
{
  j=NR;
  if(NF==6)
    for(k=1;k <= 6;k++)
      sav[(j-1) * 6 + k] = $(k) * fac;
}
END{
  printf("0 0 0\t");
  for(j=1;j <= 6;j++)
    for(k=j;k <= 6;k++){
      printf("%11g ",sav[(j-1) * 6 + k]);
    }
  printf("\n");


}
