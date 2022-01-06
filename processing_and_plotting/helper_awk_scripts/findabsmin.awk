#
# find min of the absolute value of the  column  and print the whole line
#
BEGIN{
  if(col==0)
    col=2;
  min = 1e20;
}
{
  if((substr($1,1,1)!="#")&&(NF>=col)){
      if($(col) < 0)
	  use_val = - ($(col));
      else
	  use_val = $(col);
      if(use_val < min){
	  min = $(col);
	  line = $0;
      }
  }
}
END{
    print(line);

}
