#
# determine if point is in a list of polygons, need to call with pfile defined
#
BEGIN{
  nnode=0;
  npol=0;
  ncount=0;
  np=0;
}
{
  if(FILENAME == pfile){	# reading the list of polygons
    if(!nnode)
      print("error, need to read node file first") > "/dev/stderr";
    if((substr($1,1,1)!="#")&&($1!="")){
      ncount ++;
      if(ncount == 1){		# read number of polygons
	ntpol = $1;
      }else if(ncount == 2){
	npn = $1;			# number of nodes in this polygon
	npol=0;
      }else{
# read polygon node
	npol++;
	xp[npol] = $1;
	yp[npol] = $2;
	if(npol == npn){
	  np++;			# this polygon is read in completely

	  sum = 0;
# loop through nodes and count	 
	  for(ii=1;ii<=nnode;ii++)
	    sum +=  pnpoly(npol, xp, yp, x[ii], y[ii]);
#
# print (1+n)/N
#
	  print(np,sum,nnode/ntpol);
# reset polygon
	  ncount = 1;
	}
      }
#    print(pnpoly(npol,xp,yp,$1,$2));
    }
  }else{
# read in ndoe list
    if((substr($1,1,1)!="#") && (NF >= 2)){
      nnode++;
      x[nnode] = $1;
      y[nnode] = $2;
    }
    


  }
  
}


function  pnpoly(npol, xp, yp, x, y)
{

  

  c = 0;
  i =0;
  j = npol-1;
  for (; i < npol;  ) {
    if ((((yp[i+1] <= y) && (y < yp[j+1])) ||
	 ((yp[j+1] <= y) && (y < yp[i+1]))) &&
	(x < (xp[j+1] - xp[i+1]) * (y - yp[i+1]) / (yp[j+1] - yp[i+1]) + xp[i+1])){
      c = !c;
    }
    j = i++;
  }
  return c;
}
