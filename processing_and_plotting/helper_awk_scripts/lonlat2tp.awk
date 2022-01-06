#
#
#
BEGIN{
  pif  = 57.295779513082320876798154814105;

}
{
  if((substr($1,1,1)!="#") && (NF >= 2)){
    theta = (90-$2)/pif;
    phi = $1/pif;
    printf("%.15e %.15e ",theta,phi);
    for(i=3;i<=NF;i++)
      printf("%s ",$i);
    printf("\n");
    

  }

}
END{


}