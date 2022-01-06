# read matrix and transpose
BEGIN{
    m=0;
}
{
    if(substr($1,1,1)!="#"){
	m++;
	n=NF;
	if(m>1){
	    if(n!=nold)
		print("error") > "/dev/stderr"
	}
	nold = n;
	for(i=1;i<=n;i++)
	    a[m*n+i] = $i
    }
}
END{

    for(i=1;i<=n;i++){
	for(j=1;j<=m;j++)
	    printf("%s ",a[j*n+i])
	printf("\n");
    }
	
    
}
