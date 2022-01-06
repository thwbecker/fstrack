#
# converts r,theta,phi to Cartesian
#
BEGIN{
}
{
  if($1!="" && (substr($1,1,1)!="#")){
    r=$1;
    t=$2;
    p=$3;
    
    tmp=sin(t)*r;
    x=tmp * cos(p);
    y=tmp * sin(p);
    z=cos(t)*r;

    printf("%20.16e %20.16e %20.16e ",x,y,z);
    for(i=4;i<=NF;i++)
      printf("%s ",$(i));
    printf("\n");
  }
}
END{
}
