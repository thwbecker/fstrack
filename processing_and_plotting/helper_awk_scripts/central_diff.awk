BEGIN{

}
{
    if(substr($1,1,1)!="#"){
	n++;
	t[n]=$1;
	x[n]=$2;
    }
}
END{
    for(i=2;i<n;i++){
	d1=(x[i]-x[i-1])/(t[i]-t[i-1]);
	d2=(x[i+1]-x[i])/(t[i+1]-t[i]);
	d=(d1+d2)/2;
	print(t[i],d);


    }



}
