#
#
# compute the quartiles of col col, 0.5 is median
#
# if reject_zero is set to unity, will not use zero values
#
BEGIN{
  if(col=="")
    col=1;
  if(quartile=="")
      quartile=0.5
  n=0;
}
{
  if((substr($1,1,1)!="#") && (NF>=col)){
    if((!reject_zero)||($(col) != 0)){
	if(tolower($(col))!="nan"){
	    n++;
	    x[n] = $(col);
	}
    }
  }
}
END{
  if(n){
      # sort the values
      asort(x);
      #
      nuse = int(n*quartile+.5);
      if(nuse<1)
	  nuse = 1;
      if(nuse>n)
	  nuse = n;
      print(x[nuse]);
  }
}