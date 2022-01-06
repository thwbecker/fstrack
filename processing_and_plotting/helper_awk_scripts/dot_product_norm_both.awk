#
# calculates (a \dot b)/(|a||b|) for 3-D vectors
#
# input is in format:
#
# a_1 a_2 a_3 b_1 b_2 b_3 __whatever__
#
# output is in format
# a_\dot_b/(|a||b|) __whatever__
#
{
    if((substr($1,1,1)!="#") && ($1!="") && (NF>=6)){
	an = sqrt($1*$1 + $2*$2 + $3*$3);
	bn = sqrt($4*$4 + $5*$5 + $6*$6);
	
	printf("%.8e ",($1*$4 + $2*$5 + $3*$6)/(an*bn));
	for(i=7;i<=NF;i++)
	    printf("%s ",$i);
	printf("\n");
    }
}
