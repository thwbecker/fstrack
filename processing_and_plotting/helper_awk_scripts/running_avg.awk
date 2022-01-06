#
# compute running average and std given dx spacing
#
{
    if(substr($1,1,1)!="#"){
	n++;
	x[n]=$1;
	y[n]=$2;
    }
    
}
END{
    if(dx==""){
	dx = (x[n]-x[1])/100;	# default sampling
    }


    xlast = x[1];
    jlast = 1;
    for(i=2;i<=n;i++){
	if((x[i]-xlast > dx)||(i==n)){
	    npoints=i-jlast;
	    if(npoints){
		xavg=yavg=0.;std=0.1;
		for(j=jlast;j <= i;j++){
		    xavg += x[j];
		    yavg += y[j];
		}
		xavg /= npoints;
		yavg /= npoints;
		for(j=jlast;j <= i;j++){
		    tmp = y[j]-yavg;
		    std += tmp*tmp;
		}
		std = sqrt(std/(npoints-1));
		print(xavg,yavg,std);
	    }
	    jlast = i;
	    xlast = x[i];
	}
    }
}
