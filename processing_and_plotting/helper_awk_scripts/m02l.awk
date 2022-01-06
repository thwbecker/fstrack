#
# convert moment in Nm to rupture length and mean slip
#
BEGIN{
    mu = 3e10;			# shear modulus
    alpha = 3.5e-5;
    cs = mu*alpha;			# constant times stress drop, from slip = alpha * length scaling, this needs to be 
                                        # alpha * mu
    w0 = 15e3;			# width for transition
    mc= mag(w0**3 * cs);	# transition magnitude
    print("alpha ",alpha, "stress drop times constant ",cs/1e6, " Mpa, transition mag: ",mc) > "/dev/stderr"

    
}
{
    if(substr($1,1,1)!="#"){
	m0 = $1;
	if(mag(m0) < mc){
	    # small event
	    l = (m0/(cs))**(1./3.);
	}else{
	    l = sqrt(m0/(cs*w0));
	}
	printf("%.5e %.5e ",l,alpha*l);
	for(i=2;i<=NF;i++)
	    printf("%s ",$i);
	printf("\n");
    }


}
function mag(m0x){
    magx = 2./3. * (0.4342944819032518*log(m0x)-9.1);
    return magx;
}