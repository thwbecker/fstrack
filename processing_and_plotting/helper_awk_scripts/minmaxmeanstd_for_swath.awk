BEGIN{
    n=0;
}
{
    code=$4;
    l=$3;
    e=$5;

    if(old_code != code){
	print_values(old_l,n,x);
	n=0;
    }
    n++;
    x[n]=e;
    old_code=code;
    old_l=l;
}
END{
    
    print_values(l,n,x);
}

function print_values(l,n,sum,sum2,min,max)
{
    if(n){
	asort(x);
	#
	if(n%2 !=0 )			# odd
	    median = x[(n+1)/2];
	else				# even
	    median = (x[n/2]+x[n/2+1])/2.;

	mean=0;
	min=1e60;
	max=-1e60;
	for(i=1;i<=n;i++){
	    mean += x[i];
	    if(x[i]>max)
		max=x[i];
	    if(x[i]<min)
		min=x[i];
	}
	mean /= n;

	if(n>1){
	    sum2 = 0.0;
	    for(i=1;i<=n;i++){
		x[i] -= mean;
		sum2 += x[i]*x[i];
	    }
	    std = sqrt(sum2/(n-1));
	}else{
	    std="nan";
	}
	print(min,max,mean,median,std,n);
    }
}
