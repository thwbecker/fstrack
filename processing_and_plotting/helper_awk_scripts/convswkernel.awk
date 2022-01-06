#
# 
# convert a SW A F L kernel
# 
#
BEGIN{
# these are the elastic constants from the monatgner paper
#
#> cijtoacfnl < cij.dat
#       Gc,s/L  Cc,s/N  Cc,s/A  Bc,s/A  Hc,s/F
#         0.028  -0.005  -0.002   0.034   0.003
#        -0.008   0.003   0.001  -0.002   0.006
# 
  bca=0.0342138; hcf=0.00284428; gcl=0.0277003;
  bsa=-0.00174672; hsf=0.00568852; gsl=-0.00752514;
  n=0;max=-1000;
}
{
  if(NF >= 4){
    n++;
    z[n] = $1;A = $2;F = $3; L = $4;
# 2 phi sensitivity
    sens[n] = sqrt((A * bca + F *hcf + L *gcl)**2 + (A * bsa + F *hsf + L * gsl)**2);
    if(sens[n] > max)
      max = sens[n];
  }
}
END{
  if(normalize)
    norm = max;
  else
    norm = 1.0;
  for(i=1;i<=n;i++)
    print(z[i],sens[i]/norm);

}
