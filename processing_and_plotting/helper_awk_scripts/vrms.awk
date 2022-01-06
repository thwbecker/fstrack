#
# calculate the x-x_mean RMS
#
BEGIN{
  n=0;
  if(col==0)
    col=1;
  xm=0.0;
}
{
  if((NF>=col) && 
     (substr($1,1,1)!="#") && 
     (tolower($col)!="nan") && 
     ($col!="")){
    n++;
    x[n] = sprintf("%lf",$(col));
    xm += x[n];
  }
}
END{
    if(n){
# get mean
	xm /= n;
	rms=0.0;
	for(i=1;i <= n;i++){
	    tmp = x[i] - xm;
	    rms += (tmp * tmp);
	}
	printf("%.15e\n",sqrt(rms/n));
    }else{
	print("nan");
    }
}
