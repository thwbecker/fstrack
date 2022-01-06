#
# reverse the order of a file
#
{
    n++;
    x[n]=$0;
}
END{
  for(i=n;i>0;i--)
      print(x[i]);
}
