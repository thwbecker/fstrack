#
# read a 3d symmetric tensor in upper right triangle format
# 1   2  3  4  5 6
# xx xy xz yy yz zz
#
# and print some stats
#
BEGIN{

}
{
    for(i=1;i<=6;i++)
	s[i]=$i;

    trace=(s[1]+s[4]+s[6])

    print(trace)

}


