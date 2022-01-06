#
#
# given two files pasted into this script with the same number of columns
# compute the ratio of all columns
#
{
  n=NF/2;
  for(i=1;i<=n;i++){
    if($(i+n)==0)
      printf("nan        ");
    else
      printf("%11g ",$i/$(i+n));
  }
  printf("\n");
}
