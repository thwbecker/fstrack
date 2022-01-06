BEGIN{
}
{
  if((substr($1,1,1)!="#") && (NF>=3)){
    x=$1;
    y=$2;
    z=$3;
    tmp1=x*x + y*y;
    tmp2=tmp1 + z*z;
    r=sqrt(tmp2);
    theta=atan2(sqrt(tmp1),z);
    phi=atan2(y,x);
    printf("%20.16e %20.16e %20.16e ",
	   r,theta,phi);
    for(i=4;i<=NF;i++)
      printf("%s ",$i);
    printf("\n");
  }
}
END{
}
