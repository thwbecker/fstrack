# query number of rows m and columns n
BEGIN{
    m=0;
}
{
    if(substr($1,1,1)!="#"){
	m++;
	n=NF;
	if(m>1){
	    if(n!=nold)
		print("error, length of line mismatch") > "/dev/stderr"
	}
	nold = n;
    }
}
END{
    print(m,n);
}
