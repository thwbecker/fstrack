BEGIN{
  atpi=0;
  got=0;
  checkmax=0;
  n=0;m=0;p=0;o=0;
}
{
  t=$1;
  x=$2;
  y=$3;
  z=$4;
  
  coord_to_check=y;
  poincarecoord=x;

  if(NR>startat){
  if((NR-startat)==1)
    {lastcoord_to_check=coord_to_check;lastpoincarecoord=poincarecoord;lastpiy=coord_to_check;slope=1.0;}
  else{
    if(slope != 0 && ((coord_to_check-lastcoord_to_check)/slope)>=0.0)
      { slope=(coord_to_check-lastcoord_to_check);}
    else 
      {
	if(slope>0){
	  maxtime[n]=t;
	  max[n]=coord_to_check;
	  n++;}
	else{
	  mintime[m]=t;
	  min[m]=coord_to_check;
	  m++;}
	slope=(coord_to_check-lastcoord_to_check);
      }
    lastcoord_to_check=coord_to_check;
    if((lastpoincarecoord <= 0.0 && poincarecoord > 0.0)||
       (lastpoincarecoord >= 0.0 && poincarecoord < 0.0))
      {
	poincarecoordzero[p]=coord_to_check;
	poincarecoordtime[p]=t;
	lastpoincarecoord=poincarecoord;
	p++;
      }
  }
  }
}
END{
  if(pmax)
    for(i=1;i<n;i++)
       print(maxtime[i],max[i]);
  else if(pmin)
    for(i=1;i<m;i++)
      print(mintime[i],min[i]);
  else if(pmindelay)
    for(i=2;i<m;i++)
      print(min[i-1],mintime[i-1],min[i],mintime[i]);
  else if(pmaxdelay)
    for(i=2;i<n;i++)
      print(max[i-1],maxtime[i-1],max[i],maxtime[i]);
  else if(ptimedelay)
     for(i=2;i<n-1;i++)
      print(maxtime[i]-maxtime[i-1],maxtime[i+1]-maxtime[i]);
  else if(ppoincare)
     for(i=1;i<p-1;i++)
      print(poincarecoordzero[i],poincarecoordtime[i]);
}
