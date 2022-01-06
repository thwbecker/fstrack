#
# calculate arithmetic mean of col, if col is not set, use col=1
#
# there's also avg, which can do different kinds of means
#
# $Id: mean.awk,v 1.1 2004/10/28 20:40:49 becker Exp twb $
#
BEGIN{
    sum=0.0;
  n=0;
  if(col==0)
    col=1;
}

{
    if((substr($1,1,1)!="#") && ($col!="") && (!match(tolower($col),"nan"))){
	sum += $col;
	n++;
    }
}
END{
  if(n!=0)
      printf("%20.17e\n",sum/n);
  else
      print("NaN");
}
