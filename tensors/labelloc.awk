BEGIN{
  x0=-0.05;dx=0.2;
  y0=1.9;dy=0.25;
  n=6;
}
{

  if($1 <= n){
    x=x0 + ($1-1) * dx;
    y = y0;
  }else if($1 <= 2*n){
    x=x0 + ($1-n-1) * dx;
    y = y0-dy;
  }else{
    x=x0 + ($1-2*n-1) * dx;
    y = y0-2*dy;


  }
  print(x,y);

}
