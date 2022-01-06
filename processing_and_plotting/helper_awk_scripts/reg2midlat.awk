# find the mean latitude of a GMT region
BEGIN{
    if(frac=="")
	frac = 0.5;
}
{
  n=split($1,a,"/");
  print(a[3]+(a[4]-a[3])*frac);
}
