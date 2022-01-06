#
# read a 3d symmetric tensor in upper right triangle format
# 1   2  3  4  5 6
# xx xy xz yy yz zz
#
# and print second invariant
#
BEGIN{
    if(incompressible=="")
	incompressible = 1;
}
{
    if(NF>=6 && substr($1,1,1)!="#"){
	a[1][1]=$1;
	a[1][2]=a[2][1]=$2;
	a[1][3]=a[3][1]=$3;
	a[2][2]=$4;
	a[2][3]=a[3][2]=$5;
	a[3][3]=$6;    
	
	printf("%.8e ",sec_inv(a,incompressible));
	for(i=7;i<=NF;i++)
	    printf("%s ",$i);
	printf("\n");
    }
}

function sec_inv(a, incompressible) {
    if(incompressible){
	s2 = 0;
	for(i=1;i<=3;i++)
	    for(j=1;j<=3;j++)
		s2 += a[i][j] * a[i][j];
	s2 = sqrt(s2/2);

    }else{
	s2  = a[1][1]*a[2][2] + a[1][1]*a[3][3] + a[2][2]*a[3][3];
	s2 -= a[1][2]**2 + a[1][3]**2 + a[2][3]**2;
	if(s2<0)
	    s2 = sqrt(-s2);
	else
	    s2 = sqrt(s2);
    }
    return s2;

}



