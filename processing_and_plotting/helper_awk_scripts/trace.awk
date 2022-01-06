#
# read a 3d symmetric tensor in upper right triangle format
# 1   2  3  4  5 6
# xx xy xz yy yz zz
#
# and print trace
#
BEGIN{

}
{
    if(NF>=6 && substr($1,1,1)!="#"){
	a[1][1]=$1;
	a[1][2]=a[2][1]=$2;
	a[1][3]=a[3][1]=$3;
	a[2][2]=$4;
	a[2][3]=a[3][2]=$5;
	a[3][3]=$6;    
	
	printf("%.8e ",trace(a));
	for(i=7;i<=NF;i++)
	    printf("%s ",$i);
	printf("\n");
    }
}

function trace(a) {
    return a[1][1]+a[2][2]+a[3][3];

}



