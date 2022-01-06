{
    for(i=1;i<=NF;i++){
	if(substr($i,length($i),1)=="/")
	    printf("%s ",substr($i,1,length($i)-1));
	else
	    printf("%s ",$i);
    }
    printf("\n");
}