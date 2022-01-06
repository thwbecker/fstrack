#
# convert Bernhard's plate velocity expansions to my format
# this is similar to Rick's format, but does not quite seem the same (?)
# 
BEGIN{
    if(lmax=="")
	lmax_lim=-1;
    else
	lmax_lim = lmax;
    l=m=0;
    isb=0;
    pi = atan2(1,0)*2.;
    pre_factor = 4.*pi*pi;
}
{
    if(NR==2){
	lmax = $1;
	if(lmax_lim<0)
	    lmax_lim = lmax;
	print(lmax_lim);
    }else if(NR>2){
	if(m%2)
	    fac=-pre_factor;
	else
	    fac=pre_factor;

	if(isb){
	    pb = $1;
	    tb = $2;
	    isb=0;
	    if(l<=lmax_lim)
		printf("%22.15e %22.15e %22.15e %22.15e\n",pa*fac,pb*fac,ta*fac,tb*fac);
	    m++;
	}else{
	    pa = $1;
	    ta = $2;
	    isb=1;
	}
	if(m>l){
	    l++;m=0;
	}
    }

}
END{


}
