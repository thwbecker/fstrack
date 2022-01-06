#
# convert lon1 lat1 lon2 lat2 input to meanlon meanlat
#
BEGIN{
    f=0.017453292519943296;
    twopi = 6.28318530717958647;
}
{
    if((substr($1,1,1)!="#") && (substr($1,1,1)!=">") && (NF>=4)){
	
	lambda[1]=$2*f;
	phi[1]   =$1*f;

	lambda[2]=$4*f;
	phi[2]   =$3*f;
	
	for(i=1;i<=2;i++){
	    tmp = cos(lambda[i]);
	    x[i]=tmp * cos(phi[i]);
	    y[i]=tmp * sin(phi[i]);
	    z[i]=sin(lambda[i]);
	}
	xm = (x[1] + x[2])/2;
	ym = (y[1] + y[2])/2;
	zm = (z[1] + z[2])/2;

	tmp = xm*xm + ym*ym;

	thetam=atan2(sqrt(tmp),zm);
	phim=atan2(ym,xm);
	if(phim < 0)
	    phim += twopi;
	if(phim >= twopi)
	    phim -= twopi;
	meanlon = phim/f;
	meanlat = 90.-thetam/f;

	printf("%.15e %.15e ",meanlon,meanlat)
	
	for(i=5;i<=NF;i++)
	    printf("%s ",$i);
	printf("\n");
    }
}
