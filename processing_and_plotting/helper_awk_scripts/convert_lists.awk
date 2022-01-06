BEGIN{
    n=0;
    m=9;
}
{
    if($1!=""){
	gsub(",",".");
	n++;
	x[n]=$0;
	if(n == m){
	    for(i=1;i<=m;i++)
		printf("%s ; ",x[i])
	    printf("\n");
	    n=0;
	}

    }
}