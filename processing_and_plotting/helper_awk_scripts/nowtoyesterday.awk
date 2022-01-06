{
  day=$1;month=$2;year=$3;
  if(day == 1)
    {
      if((month==2)||(month==4)||(month==6)||(month==8)||(month==9)||(month==11))
	day=31;
      else if(month == 3)
	if(int(abs(2000.0-year)/4.0) == (abs(2000.0-year)/4.0))
	  day=29;
	else 
	  day=28;
      else if(month == 1)
	{day=31;year--;}
      else
	day=30;
      month--;
    }
  else
    day--;
  printf("%02d%02d%02dtnsseismogram.gif\n",day,month,year);
}


