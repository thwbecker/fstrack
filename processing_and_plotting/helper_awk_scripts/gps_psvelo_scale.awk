#
# apply a non-linear scaling for GPS velocities in psvelo format
#
# need to define scaling for GPS (don't auto adapt for scales!)
#
#
BEGIN{
    n = 0;
    min_amp = 1e20;
    max_amp = -1e20;
    if(min_scale == "")		# for scaling
	min_scale = 1;
    if(max_scale == "")
	max_scale = 10;		# 
}
{
    if(substr($1,1,1)!="#"){
	n++;
	lon[n] = $1;
	lat[n] = $2;
	vx[n]  = $3;
	vy[n]  = $4;
	amp[n] = sqrt(vx[n]**2+vy[n]**2);
	if(amp[n] > max_amp)
	    max_amp = amp[n];
	if(amp[n] < min_amp)
	    min_amp = amp[n];
	svx[n]  = $5;
	svy[n]  = $6;
	corr[n] = $7;
	stat[n] = $8;
	for(i=9;i<=NF;i++)
	    stat[n]=sprintf("%s %s",stat[n],$i);
    }
}
END{
    # don't auto adapt to allow plotting of scales!
    print("determined range of ",min_amp," to ",max_amp,", compared to scaling from ",min_scale,
	  " to ",max_scale) > "/dev/stderr"
    if(max_amp > max_scale){
	print("error, increase max_scale") > "/dev/stderr"
    }
    else if(min_amp < min_scale){
	print("error, decrease min_scale") > "/dev/stderr"

    }else{
	range = max_scale-min_scale;
	
	for(i=1;i<=n;i++){
	    factor = (amp[i]-min_scale)/range;
	    factor = atan2(factor*7.5,1)/3.14159265358*2+0.025;
	    factor /= amp[i];
	    print(lon[i],lat[i],vx[i]*factor,vy[i]*factor,"NaN","NaN",corr[i],stat[i]);
	    
	    
	}
    }

}