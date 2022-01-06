#
# given a floating point number, return mantiss * 10 ** exponent
#
BEGIN{

}
{
  if($1 < 0){
    tval = -$1;
    sign = -1.0;
  }else{
    sign = 1.0;
    tval = $1;
  }
  exponent = (int(log(tval)/2.30258509299405+1e-7));
  val = 10**exponent;
  mant = sign * tval/val;
    
  printf("%20.15f %5i\n",mant,exponent)
}
