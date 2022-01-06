function mod (y, x){
  tmp_mod = y - x * int(y/x);
  if ( tmp_mod < 0) 
    tmp_mod += x;
  return tmp_mod;
}
