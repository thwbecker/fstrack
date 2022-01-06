{
# boundary conditions
# pin corners
if((($2 == xmin)&&(($3 == ymin)||($3 == ymax)))||
   (($2 == xmax)&&(($3 == ymin)||($3 == ymax))))
{
  print($1,$1,1,1,1);
}
else
{
# free slip on top 
  if($3 == ymax)
    {print($1,$1,1,0,1);}
# free slip and outflow for bottom
  else if(($3 == ymin)&&(lowerfree != 1))
    print($1,$1,1,0,1);
# left and right bcs
  else if($2 == xmax)
    {
      if($3 >= taperedas)
	print($1,$1,1,1,1);
      else
	{
	  if(rightfree != 1)
	    print($1,$1,1,1,0);
#	  else
#	    print($1,$1,1,0,1);
	};
    }
  else if($2 == xmin)
    {
      if($3 >= lithosphere)
	print($1,$1,1,1,1);
      else 
	{
	  if($3 >= taperedas)
	    print($1,$1,1,1,0);
	  else
	    {
	      if(leftfree != 1)
		print($1,$1,1,1,0);
#	      else
#		print($1,$1,1,0,1);
	    };
	};
    };
};
}
