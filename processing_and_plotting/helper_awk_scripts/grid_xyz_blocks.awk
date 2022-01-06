BEGIN{
  n=0;
  printf("nx: %i ny: %i reg: %s inc: %s\n",nx,ny,reg,inc) > "/dev/stderr";
  if(nx==0)
    exit;
}
{
  if(NF>=3 && (substr($1,1,1)!="#")){
    n++;
    x[n]=$1;
    y[n]=$2;
    z[n]=$3;
    if(n==nx*ny
  }


}
END{


}
