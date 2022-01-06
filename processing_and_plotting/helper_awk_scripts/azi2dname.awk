BEGIN{
    w=360/8;
}
{
    if(substr($1,1,1)!="#"){
	azi = $1;
	while(azi<0)
	    azi += 360;
	while(azi>360)
	    azi -= 360;
	if(azi < 360-7.5*w)
	    name="N";
	else if(azi < 360-6.5*w)
	    name="NE";
	else if(azi < 360-5.5*w)
	    name="E";
	else if(azi < 360-4.5*w)
	    name="SE";
	else if(azi < 360-3.5*w)
	    name="S";
	else if(azi < 360-2.5*w)
	    name="SW";
	else if(azi < 360-1.5*w)
	    name="W";
	else if(azi < 360-0.5*w)
	    name="NW";
	else
	    name="N";
	printf("%s ",name);
	for(i=2;i<=NF;i++)
	    printf("%s ",$i);
	printf("\n");
    }
}
