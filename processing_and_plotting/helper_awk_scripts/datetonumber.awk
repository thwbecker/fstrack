#
# convert day month year strings to UNIX seconds 
#

BEGIN{
  i=1;
}
{
if((substr($1,1,1)!="#"))
{
  year=$3;
  month=$2;
  day=$1;
  hour = minute = second = 0;
  the_time=sprintf("%i %i %i %i %i %i",
		   year,month,day,
		   hour,minute,second);
  secs = mktime(the_time);

  print(NR-i,secs);
}
else
  {i++;};



}
