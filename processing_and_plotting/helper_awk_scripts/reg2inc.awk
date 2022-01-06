BEGIN{
  if(dx=="")
    dx=10;
  if(dy=="")
    dy=dx;
}{

  n=split(substr($1,3),a,"/");
  xr=a[2]-a[1];
  yr=a[4]-a[3];
  printf("-I%0g/%0g\n",xr/dx,yr/dy);

}
