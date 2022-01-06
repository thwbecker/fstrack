# find the mean latitude of a GMT region and write -Acos(mlat)
BEGIN{
    if(frac=="")
	frac = 0.5;
}
{
  n=split($1,a,"/");
  mlat = (a[3]+a[4])*frac;
  printf("-A%g",cos(mlat/57.295779513));
}
