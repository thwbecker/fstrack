BEGIN{
  atpi=0;
  got=0;
  checkmax=0;
}
{
  #sign=$7/sqrt($7*$7);
  #phi=$7;
  t=$1;
  x=$2;
  y=$3;

  
  if(NR==1)
    {lasty=y;slope=1.0;}
  else{
    if(slope != 0 && ((y-lasty)/slope)>=0.0)
      { slope=(y-lasty);}
    else 
      {
	if(slope>0){
	  maxtime[n]=t;
	  max[n]=y;
	  n++;}
	else{
	  mintime[m]=t;
	  min[m]=y;
	  m++;}
	slope=(y-lasty);
      }
    lasty=y;
  }
}
END{
  if(printmin)
    for(i=0;i<m;i++)
      print(mintime[i],min[i]);
  else
    for(i=0;i<n;i++)
      print(maxtime[i],max[i]);
  


}
