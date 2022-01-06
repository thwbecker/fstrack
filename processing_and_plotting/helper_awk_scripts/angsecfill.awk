#
# count how many sectors of 0...360 or 180 are filled
#
BEGIN{
    if(orient=="")		# orient will only use 0....180
	orient=0;
    if(da=="")
	da=30;
    if(orient)
	nbin=180/da;
    else
	nbin=360/da;
    for(i=1;i<=nbin;i++)
	bin[i]=0;
}
{
    if(substr($1,1,1)!="#" && tolower($1)!="nan"){
	a = $1;
	while(a<0)
	    a+=360;
	while(a>360)
	    a-=360;
	if(orient){
	    if(a>180)
		a-=180;
	}
	i = int(a/da);
	bin[i]++;
    }
}
END{
    if(phist){			# print histogram
	for(i=1;i<=nbin;i++){
	    print((i-1)*da,i*da,bin[i]);
	}
    }
    if(pfilled){			# how many have any entries?
	j=0;
	for(i=1;i<=nbin;i++){
	    if(bin[i])
		j++;
	}
	print(j/nbin);		# fractional
    }


}
