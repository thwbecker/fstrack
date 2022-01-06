#
# convert box boundaries in r,theta,phi to GMT file with 
#
#
# assumes rmin rmax theta_min theta_max phi_min phi_max input
#
BEGIN{
  pif=57.29577951308232087679;
}
{
  rmin=$1;rmax=$2;
  ymax=90-$3*pif;ymin=90-$4*pif;
  xmin=$5*pif;xmax=$6*pif;
 
  for(r=rmin;r<=rmax+0.1;r+=(rmax-rmin)){
    y=ymin;for(x=xmin;x<=xmax;x+=0.1)print(x,y,r);print(">");
    y=ymax;for(x=xmin;x<=xmax;x+=0.1)print(x,y,r);print(">");
    x=xmin;for(y=ymin;y<=ymax;y+=0.1)print(x,y,r);print(">");
    x=xmax;for(y=ymin;y<=ymax;y+=0.1)print(x,y,r);print(">");
  }
  
  print(xmin,ymin,rmin);
  print(xmin,ymin,rmax);
  print(">");
  print(xmin,ymax,rmin);
  print(xmin,ymax,rmax);
  print(">");
  print(xmax,ymin,rmin);
  print(xmax,ymin,rmax);
  print(">");
  print(xmax,ymax,rmin);
  print(xmax,ymax,rmax);
  print(">");


 
}
