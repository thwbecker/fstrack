#
# mean and fast, inaccurate stddev for all columns
#
# $Id: meanandstddevallcol.awk,v 1.2 2005/07/16 00:14:14 becker Exp becker $
#
BEGIN{
}
{
    if(substr($1,1,1)!="#"){
	if(NF > nfmax)
	    nfmax = NF;
	for(i=1;i <= NF;i++){
	    if(tolower($i) != "nan"){
		sum[i] += $i;
		sum2[i] += $i * $i;
		n[i]++;
	    }
	}
    }
}
END{
    for(i=1;i <= nfmax;i++){
	if(n[i] < 2)
	    printf("NaN ");
	else{
	    stddev = sqrt (1e-16+ ((n[i] * sum2[i] - sum[i] * sum[i]) / ((n[i]*(n[i]-1.0)))  ));
	    printf("%g %g\t\t",sum[i]/n[i],stddev);
	}
    }
    printf("\n");
}
