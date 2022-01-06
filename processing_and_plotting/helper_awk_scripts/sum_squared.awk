#
# sum the squares
# 
BEGIN{
  n=0;
  if(col==0)
    col=1;
  sum =0.;
}
{
  if((NF>=col)&&(substr($1,1,1)!="#")&&
     (tolower($col)!="nan")&&($col!="")){
    n++;
    sum += $(col) * $(col);
  }
}
END{
    printf("%.15e\n",sum);

}
