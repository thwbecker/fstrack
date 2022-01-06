#
# converts input in format
# lon lat v_x v_y ...
# to 
# lon lat azimuth(clockwise in degrees from north) length ....
#
BEGIN{
  f=57.295779513082321;
}
{
  if((substr($1,1,1)!="#") && $1!="" && NF>=4){
    vx=$3;
    vy=$4;
    azi=atan2(vx,vy)*f;
    if(azi<0)
      azi+=360.;
    v=sqrt(vx*vx+vy*vy);
    printf("%g %g %g %g ",$1,$2,azi,v);
    for(i=5;i<=NF;i++)
      printf("%s ",$i);
    printf("\n");
  }
}
