#
# compute centroid of x y data
#
BEGIN{
    x=y=0;
}
{
    if((substr($1,1,1)!="#")&&(NF>=2)){
	n++;
	x += $1;
	y += $2;
    }
}
END{
    if(n)
	print(x/n,y/n);
}