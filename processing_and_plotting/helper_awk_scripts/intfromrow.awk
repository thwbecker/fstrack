#
# provide interpolation weights from a row of ascending values
#
BEGIN{
    if(z=="")print("error, z should be defined") > "/dev/stderr"
}
{
    if(NR==1){
	for(i=1;i<=NF;i++){
	    x[i] = $i;
	    if(i>1){
		if(x[i] <= x[i-1])
		    print("error, row values should increase") > "/dev/stderr"
		
	    }
	}
	n=NF;
    }
}
END{
    j=1;
    while( (j<n) && (x[j] < z))
	j++;
    if(j==1)
	j=2;
    i=j-1;
    
    fac=(z - x[i])/(x[j]-x[i]);
    fac2=1.0-fac;
#
# output of weight_1 val_1 weight_2 val_2
#
    print(fac2,x[i],fac,x[j]);
    
    
}