BEGIN{

}
{
  x=$1;
  y=$2;
#
  epp=$3;etp=$4;ett=$5;
  e1=$6;e2=$7;eazi=$8;
#
  err=$9;ert=$10;erp=$11;

  fie = ett + err + epp;
 
  sie  = ett*epp; 
  sie += ett*err;
  sie += epp*err;
  sie -= etp*etp;
  sie -= ert*ert;
  sie -= erp*erp;


  print(x,y,sie);
}
