#
# convert a format of 4:13:20 to seconds
#
{
  n=split($1,a,":");
  if(n>3)
    print("error") > "/dev/stderr";
  if(n==3)
    print(a[1]*3600+a[2]*60+a[3]);
  else if(n==2)
    print(a[1]*60+a[2]);
  else
    print(a[1]);
  

}
