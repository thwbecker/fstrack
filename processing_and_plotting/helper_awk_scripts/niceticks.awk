#
# compute nicely spaced xticks given 
#
# min max steps
BEGIN{
#    f=1.0/ln(10.0);
    f=0.4342944819032518;
}
{
    tempStep = ($2-$1)/$3;
    log_s = log(tempStep)*f;	# log10
    # floor
    if(log_s > 0)
	mag = int(log_s);
    else
	mag = int(log_s-.5);
    # 
    magPow = 10**mag;


    magMsd = int(tempStep/magPow + 0.5);

    
    if (magMsd > 5.0)
	magMsd = 10.0;
    else if (magMsd > 2.0)
	magMsd = 5.0;
    else if (magMsd > 1.0)
	magMsd = 2.0;
    else
	magMsd = 1;

    print magMsd*magPow;
}

