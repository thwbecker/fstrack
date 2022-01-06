# find the mean longitude of a GMT region
# and make a grid of even sampled points
BEGIN{
    if(dx=="")
	dx=1;
    if(dy=="")
	dy = dx;
    if(offset=="")
	off=0;
    else
	off=(dx+dy)/4;
    if(add=="")
	add=0;
}
{
    split($1,a,"/");
    w = sprintf("%f",substr(a[1],3));
    e=sprintf("%f",a[2]);
    s=sprintf("%f",a[3]);
    n=sprintf("%f",a[4]);

    for(x=w+off-add;x <= e+1e-12+add;x+=dx)
	for(y=s+off-add;y <= n+1e-12+add;y+=dy)
	    printf("%g %g\n",x,y);

}
END{}
