# convert 3x3 matrix to symmetric matrix in six component upper right triangle form
BEGIN{


}
{
    if((substr($1,1,1)!="#")&&(NF>=9)){
	x[1*3+1] = $1;x[1*3+2] = $2;x[1*3+3] = $3;
	x[2*3+1] = $4;x[2*3+2] = $5;x[2*3+3] = $6;
	x[3*3+1] = $7;x[3*3+2] = $8;x[3*3+3] = $9;
	for(i=1;i<=3;i++)
	    for(j=1;j<=3;j++)
		s[i*3+j] = 0.5 * (x[i*3+j] - x[j*3+i]);
	print(s[1*3+1],s[1*3+2],s[1*3+3],s[2*3+2],s[2*3+3],s[3*3+3]);
    }
}