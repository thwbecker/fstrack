#
# convert my plate velocity expansions to Bernhard's format
#
BEGIN{
    if(lmax=="")
	lmax_lim = -1;
    else
	lmax_lim = lmax;
    l=m=0;
    pi = atan2(1,0)*2.;
    pre_factor = 4.*pi*pi;
}
{
    if(NR==1){
	lmax = $1;
	if(lmax<lmax_lim){
	    print("error ",lmax,lmax_lim) > "/dev/stderr";
	    exit
	}
	if(lmax_lim < 0)
	    lmax_lim = lmax;
	
	print("velocity field spherical harmonic expansion")
	printf("%i =degree, cspol, cstor: \n",lmax_lim);

    }else{
	pa=$1;pb=$2;ta=$3;tb=$4;
	if(l <= lmax_lim){
	    #print(l,m,"A",pa,ta);print(l,m,"B",pb,tb);
	    if(m%2)
		fac=-pre_factor;
	    else
		fac=pre_factor;
	    
	    printf("%22.15e %22.15e\n",pa/fac,ta/fac);
	    printf("%22.15e %22.15e\n",pb/fac,tb/fac);
	    
	    m++;
	    if(m>l){
		l++;m=0;
	    }
	}
    }

}
END{


}
