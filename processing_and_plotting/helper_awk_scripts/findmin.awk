#
# find min in column  and print the whole line
#
BEGIN{
  if(col==0)
    col=2;
  min = 1e20;
}
{
  if((substr($1,1,1)!="#")&&(NF>=col)){
      
      if($(col) < min){
	  min = $(col);
	  line = $0;
      }
  }
}
END{
    print(line);

}
