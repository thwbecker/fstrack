BEGIN{

}
{
  if($1!="" && (substr($1,1,1)!="#")){
    x=$1;
    y=$2;
    z=$3;
    tx=$4;
    ty=$5;
    tz=$6;
    
    printf("%20.16e %20.16e %20.16e\n", y*tz-z*ty,
	   z*tx-x*tz,
	   x*ty-y*tx);
  }
}
