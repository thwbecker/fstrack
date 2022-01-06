#
# reads k l wave_length[km] A1 B1 A2 B2
# and computes correlation as a function of wavelength
#
BEGIN{
  n=0.0;
  old_freq=-1;
}
{
    if(($3 != old_freq)&&(NR!=1)){
	print(old_freq,corr(x,y,n));
	n=0.0;
    }
    n += 1.0;
    x[n]=$4;
    y[n]=$6;

    n += 1.0;
    x[n]=$5;
    y[n]=$7;


    old_freq = $3;
}
END{
    print(old_freq,corr(x,y,n));
}

#
# given x and y, compute linear correlation
#
function corr(x,y,n)
{
    if(n >= 2){
# calculate means
	mx=0.0;
	my=0.0;
	for(i=1;i<=n;i++){
	    mx += x[i];
	    my += y[i];
	}
	mx /= n;
	my /= n;
# calculate correlation
	s1=0.0;s2=0.0;s3=0.0;
	for(i=1;i<=n;i++){
	    dx=x[i]-mx;
	    dy=y[i]-my;
	    s1 += (dx*dy);
	    s2 += (dx*dx);
	    s3 += (dy*dy);
	}
	tmp=sqrt(s2)*sqrt(s3);
	if(tmp != 0.0)
	    loc_r = s1/tmp;
	else
	    loc_r = "NaN"
    }else{
	loc_r = "NaN";
    }
    return loc_r;

}
