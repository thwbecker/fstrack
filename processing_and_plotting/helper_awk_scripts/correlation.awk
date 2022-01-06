#
# calculate the linear correlation coefficient
# input is two columns, x_i and y_i
#
#
# output is r_corr t_student N-3
#
BEGIN{
  n=0.0;
}
{
  if((NF>=2) && (substr($1,1,1)!="#") && (tolower($1)!="nan") && (tolower($2)!="nan")){
    n += 1.0;
    x[n]=$1;
    y[n]=$2;
  }
}
END{
    r = corr(x,y,n);
# stat quantity
    tmp = 1.0 - r*r;
    if(tmp > 0){
	t=r*sqrt((n-2.)/tmp);
    }else{
	if(r > (1.0-1e-14))
	    t=1e30;
	else
	    t=0.0;
    }
    print(r,t,n-2);
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
