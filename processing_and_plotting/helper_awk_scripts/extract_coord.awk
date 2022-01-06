{
if(NR==1){
  nx=$1;ny=$2;nz=$3;
}else{
  i=NR-1;
  if(i>=1 && i<=nx)
    if(px)
      print($1>=0?$1:360+$1);
  if(i>=nx+1 && i<=nx+ny)
    if(py)
      print($1);
  if(i>=nx+ny+1 && i<=nx+ny+nz)
    if(pz)
      print($1);
}
}
