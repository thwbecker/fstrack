#
# find the maximum in an x y series and print the x_max and y_max values
#
#
BEGIN{
  if(col=="")			# which column to use for y
      col=2;
  if(interpolate=="")		# interpolate the x and y values?
      interpolate = 0;
  if(plotrow == "")		# print out the whole row with the max?
      plotrow = 0;
  maxinit=-1e30;
  max=maxinit;
  n=0;
}
{
  if((substr($1,1,1)!="#")&&(NF>=col)){
    if(tolower($(col))!="nan"){
	n++;
	x[n] = $1;
	y[n] = $(col);
	row[n] = $0;
    }
  }
}
END{
    for(i=1;i<=n;i++){
	if(y[i] > max){
	    max = y[i];
	    maxtime = x[i];
	    maxrow = row[i];
	}
    }
    if(max == maxinit){
	print("no data found") > "/dev/stderr";
	print("nan");
    }else{
	if(plotrow)
	    print(maxrow)
	else
	    print(maxtime,max);
    }
}
