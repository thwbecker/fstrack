#
# convert vera's Cij format 
# i j cij
# ....
#
# to mine
# x y z c11 c12 c13 c14 .... c66
#
BEGIN{


}
{
  if((NF==3)&&(substr($1,1,1)!="#")){
    j=$1;k=$2;
    sav[(j-1) * 6 + k] = $3;
    sav[(k-1) * 6 + j] = $3;
  }
}
END{
  printf("0 0 0\t");
  for(j=1;j <= 6;j++)
    for(k=j;k <= 6;k++){
      printf("%11g ",sav[(j-1) * 6 + k]);
    }
  printf("\n");
}
