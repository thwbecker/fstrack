#
# computes the geometric mean x = ( Prod_(i=1..n) x_i )^(1/n)
#
#
BEGIN{
  uselog=1;			# use the log method to actually compute the mean
  sum = 0.0;
  prod = 1.0;
  n=0;
}
{
  if((substr($1,1,1)!="#")&&(NF>=1)){
    x = $1;
    if(x < 0)
      print("error, x shouldn't be negative for geometric mean") > "/dev/stderr";
    if(uselog)
      sum += log(x);
    else
      prod *= x;
    n++;
  }
}
END{
  if(uselog)
    print(exp(1./n*sum));
  else
    print(prod**(1/n));

}
