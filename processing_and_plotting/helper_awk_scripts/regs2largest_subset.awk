#
#
# given two -R style gmt regions, print the largest subset
BEGIN{}
{
    split($1,r1,"/");
    split($2,r2,"/");
    w[1]=substr(r1[1],3);e[1]=r1[2];s[1]=r1[3];n[1]=r1[4];
    w[2]=substr(r2[1],3);e[2]=r2[2];s[2]=r2[3];n[2]=r2[4];
    if(w[1]<w[2])wuse=

    print(w[1],e[1],s[1],n[1],w[2],e[2],s[2],n[2])
}
