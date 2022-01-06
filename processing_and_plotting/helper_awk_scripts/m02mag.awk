# convert moment in Nm to magnitude

{
    if(substr($1,1,1)!="#"){
	mag = 2./3. * (0.4342944819032518*log($1)-9.1);
	printf("%.5f ",mag);
	for(i=2;i<=NF;i++)
	    printf("%s ",$i);
	printf("\n");

    }
}