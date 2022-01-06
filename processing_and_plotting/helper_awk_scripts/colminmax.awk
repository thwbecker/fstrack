# script to convert the data range of colorscale files for GMT
# by rescaling
#
# part of the iGMT distribution
# 
# Thorsten Becker, 04/05/99
#
BEGIN{
  if(min==0)min=-1.;
  if(max==0)max=1.;
  mean=(max+min)/2.;
  range=(max-min)/2.;
  printevery=25;
}
{
    xl[NR]=$1;rl[NR]=$2;gl[NR]=$3;bl[NR]=$4;xr[NR]=$5;rr[NR]=$6;gr[NR]=$7;br[NR]=$8;
}
END{
  for(i=1;i<=255;i++){
    if(reverse)
      j=256-i;
    else
      j=i;
    if(printevery - NR > 0)
      print(xl[i]*range+mean,rl[j],bl[j],bl[j],xr[i]*range+mean,rr[j],br[j],gr[j]);
    else {
      print(xl[i]*range+mean,rl[j],bl[j],bl[j],xr[i]*range+mean,rr[j],br[j],gr[j],"L");
      printevery=NR+printevery;
    }
  }

}
