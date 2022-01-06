{
if(($3 <  taperedas)&&($3 != ymin)&&((($2 == xmax)&&(rightfree != 1))||(($2 == xmin)&&(leftfree != 1))))
{
  print($1,0,fliflo,0.0);
};

if((leftfree != 1)&&($2 == xmin)&&($3 < lithosphere)&&($3 >  taperedas))
{
  print($1,0,(lithosphere-$3)/(lithosphere-taperedas)*fliflo,0.0);
};
     
}
