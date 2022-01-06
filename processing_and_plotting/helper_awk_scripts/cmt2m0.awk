#
# convert cmt solution format to m0 [Nm]
#
BEGIN{


}
{
    # read
    #
    # lon lat depth mrr mtt mpp mrt mrp mtp  exp
    #
    x=$1;
    y=$2;
    z=$3
    m[1*3+1]=$4;		# rr
    m[2*3+2]=$5;		# tt
    m[3*3+3]=$6;		# pp
    m[1*3+2]=m[2*3+1]=$7;	# rt
    m[1*3+3]=m[3*3+1]=$8;	# rp
    m[2*3+3]=m[3*3+2]=$9;	# tp
    
    
    mexp=$10 - 7;		# go from dyn.cm to Nm
    #
    # tensor norm
    #
    scalar_mom = 0;
    for(i=1;i <= 3;i++)
	for(j=1;j <= 3;j++)
	    scalar_mom += m[i*3+j] * m[i*3+j];
    scalar_mom =  sqrt(scalar_mom);

    scalar_mom /= sqrt(2);


    m0 = 10**mexp * scalar_mom;

    printf("%g %g %g %e\n",x,y,z,m0);
}