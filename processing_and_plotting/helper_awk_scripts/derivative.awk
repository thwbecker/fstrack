#
# given x y data, compute first derivative by central derivatives
#
#
BEGIN{
    n=1;
    if(order=="")
	order = 1;		# default is first derivative
    x[n]=0;y[n]=0;
}
{
    if((substr($1,1,1)!="#") && (NF >= 2)){
	n++;
	x[n] = $1;y[n]=$2;
	if(n>3)
	    if(x[n]<x[n-1])
		print("ordering error") > "/dev/stderr"
    }
}
END{
    # sort 

    
    if(n > 2){
    # make up end points
	h = x[3] - x[2];dy = (y[3]-y[2])/h;
	x[1] = x[2] - h;y[1] = y[2] - h * dy;
	h = x[n]-x[n-1];dy = (y[n]-y[n-1])/h;
	x[n+1] = x[n] + h;y[n+1] = y[n] + h * dy;
    }


    if(order == 1){
	i1 = 2;i2 = n;
	for(i=i1;i <= i2;i++){
	    dx1 = (y[i+1]-y[i])/(x[i+1]-x[i]); # forward 
	    dx2 = (y[i]-y[i-1])/(x[i]-x[i-1]); # backward
	    d[i] = (dx1+dx2)/2;
	}
    }else if(order == 2){			# second derivative
	i1 = 2;i2 = n;
	for(i=i1;i <= i2;i++){		       # compute central differences around point in center
	    dx1 = (y[i+1]-y[i])/(x[i+1]-x[i]); # forward 
	    dx2 = (y[i]-y[i-1])/(x[i]-x[i-1]); # backward
	    h = (x[i+1] + x[i])/2 - (x[i] + x[i-1])/2
	    d[i] = (dx1-dx2)/h;
	}
    }else if(order == 0){	# for testing
	i1 = 1;i2 = n+1;
	for(i=i1;i<=i2;i++)
	    d[i] = y[i];
    }else{
	i1=2;i2=1;
    }
    for(i=i1;i <= i2;i++){
	printf("%.15e %.15e\n",x[i],d[i]);
    }
}