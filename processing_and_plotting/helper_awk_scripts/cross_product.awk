#
# calculates a \cross b
#
# input is in format:
# a_1 a_2 a_3 b_1 b_2 b_3
#
{
if((substr($1,1,1)!="#") && $1!=""){
  a[1]=$1;a[2]=$2;a[3]=$3;
  b[1]=$4;b[2]=$5;b[3]=$6;
  printf("%g %g %g ",
	 ( b[3]* a[2] - b[2]* a[3]),
	 ( b[1]* a[3] - b[3]* a[1]),
	 ( b[2]* a[1] - b[1]* a[2]));
  for(i=7;i<=NF;i++)
    printf("%s ",$i);
  printf("\n");
}
}
