BEGIN{

    sum=0;
    n=0;
    if(col=="")
	col=2;
}
{
    if((substr($1,1,1)!="#")&&(NF>=col)){
	if(tolower($(col))!="nan"){
	    n++;
	    x[n]=$1;
	    y[n]=$(col);
	    sum += y[n];
	}
    }
}
END{
    mean = sum/n;
    
    for(i=1;i<=n;i++)
	print(x[i],y[i]-mean);

}
