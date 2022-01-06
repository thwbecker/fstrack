#
# converts input in format
# lon lat v_r v_theta v_phi ....
# to 
# lon lat azimuth(clockwise in degrees from north) length ....
#
# (v_r is discarded)
#
BEGIN{
  f=57.295779513082321;
}
{
  if(((substr($1,1,1)!="#")) && ($1!="") && (NF>=5)){
    vx= $5;
    vy=-$4;
    azi=atan2(vx,vy)*f;
    if(azi<0)
      azi+=360;
# check for numerical inaccuracies
    tmp=vx*vx+vy*vy;
    if(tmp<=0)
      v=0.0;
    else
      v=sqrt(tmp);
    printf("%g %g %g %g ",$1,$2,azi,v);
    for(i=6;i<=NF;i++)
      printf("%s ",$i);
    printf("\n");
  }
}
