#
# calculate weighted mean, input is in format
# 
# weight x_1 x_2 ....
#
# if col is set, will select x_col (not the actual col number)
#
BEGIN{
    maxnf=1;
}
{
    if((substr($1,1,1)!="#") && (NF >=2) ){
	w = $1;
	for(i=2;i <= NF;i++){
	    if(!match(tolower($(i)),"nan")){
		sumw[i-1] += w;
		sum[i-1]  += w * $(i);
	    }
	}
	if(NF-1 > maxnf)
	    maxnf=NF-1;
    }
}
END{
    if(col == ""){
	for(i=1;i<=maxnf;i++){
	    if(sumw[i] != 0)
		printf("%.15e ",sum[i]/sumw[i]);
	    else
		printf("NaN ");
	}
	printf("\n");
    }else{
	if((col > maxnf) || (col < 1))
	    print("NaN");
	else
	  printf("%.15e\n",sum[col]/sumw[col]);  

    }
}
