#
# calculate sum of column col
# if col not set, use col=1
#
# $Id: sum.awk,v 1.1 2015/03/12 00:18:23 becker Exp becker $
#
BEGIN{
  s=0.0;
  n=0;
  if(col==0)
      col=1;
}
{
    if((substr($1,1,1)!="#") && ($(col)!="") && (tolower($col)!="nan")){
	s += $(col);
	n++;
    }
}
END{
    if(n)
	printf("%20.15e\n",s);
    else
	print("NaN")
}
