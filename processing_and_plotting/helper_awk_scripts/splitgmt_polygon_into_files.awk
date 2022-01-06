BEGIN{
    n=m=0;
}
{
    n++;x[n]=$1;y[n]=$2;
    if($1==">"){
	m++;
	filename=sprintf("%s.%i",FILENAME,m);
	print(x[1],y[1]) > filename;
	for(i=2;i<n;i++)
	    print(x[i],y[i]) >> filename;
	n=0;
    }
}