BEGIN{
    if(col=="")
	col=1
}
{
    if(substr($1,1,1)!="#"){
	printf("%s ",$(col));
    }
}
END{
    printf("\n");
}
