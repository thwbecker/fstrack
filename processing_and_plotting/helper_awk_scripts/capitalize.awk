BEGIN{
}{
    for(i=1;i<=NF;i++){
	a = $i;
	printf("%s%s",toupper(substr(a,1,1)),substr(a,2));
    }
    printf("\n");
}
