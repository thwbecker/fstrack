#
# reverse the order of a column
#
BEGIN{
}
{
  for(i=NF;i>=1;i--)
    printf("%s ",$i);
  printf("\n");
}
