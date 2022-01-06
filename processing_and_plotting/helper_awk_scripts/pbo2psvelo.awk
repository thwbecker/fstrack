BEGIN{
    #
    # original velocities are in m/yr
    #
    if(mode == "")
	mode=1;
    if((mode == 1)||(mode == 2)) #  horizontal, convert to cm/yr
	scale = 100;
    else			# mode 3, verticals, convert to mm/yr
	scale = 1000;
}
{
    #
    # change the sign of the up velocity to be positive upward
    #
    if((NF==30)&&($1 != "*Dot#")){

	lon=$9;lat=$8;
	height=$10;		# 

	vn=$20*scale;		# north
	ve=$21*scale;		# east
	vu=$22*scale;		# up

	svn=$23*scale;		# errors
	sve=$24*scale;
	svu=$25*scale;
	corrne=$26;		# N-E correlation
	code=$1;		# station code
	
	if((mode == 1)||(mode == 3))
	    printf("%.2f %.2f %g %g %g %g %g %s\n",lon,lat,ve,vn,sve,svn,corrne,code);
	else if(mode == 2){
	    if($25 != "9.99999")
		print(lon,lat,vu,svu)
	}

    }
}
