#
# normalized tensor dot product
#
# reads:
#
# axx axy axz ayy ayz azz   bxx bxy bxz byy byz bzz
#
# and prints the normalized dot product
#
# a.b/|a||b| = cos(theta)
#
BEGIN{

}
{
    if(NF>=12 && substr($1,1,1)!="#"){
	a[1][1]=$1;
	a[1][2]=a[2][1]=$2;
	a[1][3]=a[3][1]=$3;
	a[2][2]=$4;
	a[2][3]=a[3][2]=$5;
	a[3][3]=$6;    

	b[1][1]=$7;
	b[1][2]=b[2][1]=$8;
	b[1][3]=b[3][1]=$9;
	b[2][2]=$10;
	b[2][3]=b[3][2]=$11;
	b[3][3]=$12;

	dot = 0;
	for(i=1;i<=3;i++)
	    for(j=1;j<=3;j++)
		dot += a[i][j] * b[i][j];
	ndot = dot/tnorm(a)/tnorm(b);
	
	printf("%.8e ",ndot);
	for(i=13;i<=NF;i++)
	    printf("%s ",$i);
	printf("\n");
    }
}

function tnorm(a) {

    norm = 0;
    for(i=1;i<=3;i++)
	for(j=1;j<=3;j++)
	    norm += a[i][j] * a[i][j];
    norm = sqrt(norm);
    
    return norm;

}



