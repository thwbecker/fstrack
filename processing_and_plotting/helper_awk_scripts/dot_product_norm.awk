#
# calculates (a \dot b)/|b|^2 for 3-D vectors
#
# input is in format:
#
# a_1 a_2 a_3 b_1 b_2 b_3 __whatever__
#
# output is in format
# a_\dot_b/|b|^2 __whatever__
#
{
if((substr($1,1,1)!="#") && $1!=""){
  bn = $4*$4+$5*$5+$6*$6;
  printf("%g ",($1*$4+$2*$5+$3*$6)/bn);
  for(i=7;i<=NF;i++)
    printf("%s ",$i);
  printf("\n");
}
}
