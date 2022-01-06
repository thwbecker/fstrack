BEGIN{
    if(e<0)e+=180;
    if(w<0)w+=180;
}
{
    x=$1;if(x<0)x+=180;
    y=$2;
    if((x>=w) && (x<=e) && (y>=s) && (y<=n))
	print($0,1);
    else 
	print($0,0);
}