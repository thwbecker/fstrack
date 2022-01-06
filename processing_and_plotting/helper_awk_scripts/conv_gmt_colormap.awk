# new to old GMT colormap conversion
BEGIN{
}
{
    if((substr($1,1,1)!="#")&&(NF>=4)){
	n++;
	v1=$1;
	split($2,col1,"/");
	v2=$3;
	split($4,col2,"/");
	if(n==1){
	    first_col[1] = col1[1];
	    first_col[2] = col1[2];
	    first_col[3] = col1[3];
	}
	printf("%g %i %i %i %g %i %i %i\n",
	       v1,col1[1],col1[2],col1[3],
	       v2,col2[1],col2[2],col2[3]);
    }

}
END{
    last_col[1] = col2[1];
    last_col[2] = col2[2];
    last_col[3] = col2[3];

    printf("B %i %i %i\n",first_col[1],first_col[2],first_col[3]);
    printf("F %i %i %i\n",last_col[1],last_col[2],last_col[3]);
    print("N 128 128 128");

}
