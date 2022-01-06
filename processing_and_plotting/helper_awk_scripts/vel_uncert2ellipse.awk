#
# convert sig_x sig_y correlation from velocity uncertainties to 
# major minor angle formulation for ellipses
#
BEGIN{
    pif = 57.2957795131;
}
{
    if((substr($1,1,1)!="#") && (NF >= 3)){
	sigx = $1;
	sigy = $2;
	rho = $3;

	#conrad = sqrt( -2.0 * log(1.0 - confid)); 
	conrad = 1.0;
	
	a = sigy*sigy - sigx*sigx;a *= a;
	b = rho*sigx*sigy;b *= 4 * b;

	c = sigx*sigx + sigy*sigy;

	# minimum eigenvector (semi-minor axis) 
	eigen1 = conrad * sqrt((c - sqrt(a + b))/2.);

	# maximu eigenvector (semi-major axis) */
	eigen2 = conrad * sqrt((c + sqrt(a + b))/2.);

	d = 2. * rho * sigx * sigy;
	e = sigx*sigx - sigy*sigy;

	ang = atan2(d,e)/2.;
	
	# output is e1 e2 angle (to e1, CCW from x in deg)
	# (for plotting with psxy -Se)
	print(eigen1,eigen2,ang*pif-90)
    }
}