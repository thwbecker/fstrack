#
# convert Sav type stiffness tensor file to 6x6
#
BEGIN{
}
{
  if((substr($1,1,1)!="#")&&(NF>=24)){
    lon=$1;lat=$2;d=$3;
    os=4;
    for(j=1;j <= 6;j++)
      for(k=j;k <= 6;k++){
	sav[(j-1) * 6 + k] = $os;
	sav[(k-1) * 6 + j] = $os;
	os++;
      }
    for(i=1;i<=6;i++){
      for(j=1;j<=6;j++)
	printf("%11g ",sav[(i-1) * 6 + j]);
      printf("\n");
    }
  }
}
END{

}

function vera_out()
{
#  print("# Cij at lon ",lon," lat ",lat, " z ",d);
#  print("i j Cij");
#  print("");
  printf("%2i %i %7.2f\n",1,1,sav[(1-1)*6+1]);
  printf("%2i %i %7.2f\n",2,2,sav[(2-1)*6+2]);
  printf("%2i %i %7.2f\n",3,3,sav[(3-1)*6+3]);
  printf("%2i %i %7.2f\n",4,4,sav[(4-1)*6+4]);
  printf("%2i %i %7.2f\n",5,5,sav[(5-1)*6+5]);
  printf("%2i %i %7.2f\n",6,6,sav[(6-1)*6+6]);
  printf("%2i %i %7.2f\n",1,2,sav[(1-1)*6+2]);
  printf("%2i %i %7.2f\n",1,3,sav[(1-1)*6+3]);
  printf("%2i %i %7.2f\n",1,4,sav[(1-1)*6+4]);
  printf("%2i %i %7.2f\n",1,5,sav[(1-1)*6+5]);
  printf("%2i %i %7.2f\n",1,6,sav[(1-1)*6+6]);
  printf("%2i %i %7.2f\n",2,3,sav[(2-1)*6+3]);
  printf("%2i %i %7.2f\n",2,4,sav[(2-1)*6+4]);
  printf("%2i %i %7.2f\n",2,5,sav[(2-1)*6+5]);
  printf("%2i %i %7.2f\n",2,6,sav[(2-1)*6+6]);
  printf("%2i %i %7.2f\n",3,4,sav[(3-1)*6+4]);
  printf("%2i %i %7.2f\n",3,5,sav[(3-1)*6+5]);
  printf("%2i %i %7.2f\n",3,6,sav[(3-1)*6+6]);
  printf("%2i %i %7.2f\n",4,5,sav[(4-1)*6+5]);
  printf("%2i %i %7.2f\n",4,6,sav[(4-1)*6+6]);
  printf("%2i %i %7.2f\n",5,6,sav[(5-1)*6+6]);
}
