#
# given a list of y values, interpolate to y1
# 
BEGIN{
  n=0;
}
{
    if(substr($1,1,1)!="#"){
	n++;
	y[n]=$1;
	if(n>1)
	    if(y[n]<y[n-1])
		print("error, assuming sorted ascendingly") > "/dev/stderr";
    }
}
END{
    
    j=1;
    while( (j<n) && (y[j] < y1))
	j++;
    if(j==1)
	j=2;
    i=j-1;
    
    fac=(y1 - y[i])/(y[j]-y[i]);
    fac2=1.0-fac;
    print(fac,j,y[j],fac2,i,y[i])
}

