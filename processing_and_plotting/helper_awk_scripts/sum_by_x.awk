BEGIN{
    
    if(avg=="")
	avg = 0;
}
{
    n=int(($1-xmin)/dx)+1;
    if(n > 0){
	b[n] += $2;
	num[n] ++;
	if(n > nmax)
	    nmax = n;
	sum += $2;
    }

}
END{
    if(avg==0){
	for(i=1;i <= nmax;i++){
	    print(dx/2+xmin+(i-1)*dx,b[i],b[i]/sum,num[i])
	    sum2 += b[i];
	}
	print(sum,sum2) > "/dev/stderr"
    }else{
	for(i=1;i <= nmax;i++){
	    if(num[i])
		print(dx/2+xmin+(i-1)*dx,b[i]/num[i],b[i]/num[i])
	}

    }
}
