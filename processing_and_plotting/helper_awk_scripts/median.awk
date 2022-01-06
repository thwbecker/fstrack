#
#
# compute the median of col col
# if reject_zero is set to unity, will not use zero values
#
BEGIN{
  if(col=="")
    col=1;
  n=0;
}
{
    if((substr($1,1,1)!="#")&&(NF>=col)&&(tolower($(col))!="nan")){
	if((!reject_zero)||($(col) != 0)){
	    n++;
	    x[n] = $(col);
	}
    }
}
END{
  if(n){
      # sort the values
      asort(x);
      #
      if(n%2 !=0 )			# odd
	  print(x[(n+1)/2]);
      else				# even
	  print((x[n/2]+x[n/2+1])/2.);
  }
}
