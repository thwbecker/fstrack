# convert Mw to M_0 moment in Nm
{
    if(substr($1,1,1)!="#"){
	m0 = 10**(3./2.*$1 + 9.05);
	printf("%.6e ",m0);
	for(i=2;i<=NF;i++)
	    printf("%s ",$i);
	printf("\n");
    }

}

