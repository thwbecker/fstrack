BEGIN{
  if(x0==0.0)print("x0=0, error?");
  if(x1==0.0)print("x1=0, error?");
  if(y0==0.0)print("y0=0, error?");
  if(y1==0.0)print("y1=0, error?");
  uy=y1-y0;
  ux=x1-x0;
}
{
  x=$1;
  y=$2;
  vy=y-y0;
  vx=x-x0;
  if((ux*vy-uy*vx)>=0.0)
    print(1,$0);
  else 
    print(0,$0);
}
