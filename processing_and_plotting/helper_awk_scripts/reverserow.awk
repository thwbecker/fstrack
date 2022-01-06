BEGIN{
}
{
    for(i=1;i<=NF;i++)
	a[i]=$i;
    for(i=NF;i>=1;i--)
	printf("%s ",a[i]);
    printf("\n");
    
}
