#
# convert Aki type focal mechanisms into moment
#
#     input format 
#
#     lon, lat, depth, strike, dip, rake, mag
#
BEGIN{
  f = 57.2957795130823;
  log10_fac=0.4342944819032518;
}
{
  if((substr($1,1)!="#")&&(NF>=7)){
    lon = $1;
    if(lon<0)
	lon+=360.0;
    lat = $2;
    depth = $3;
# angles
    phi = $4/f;			        # strike
    delta = $5/f;			# dip
    lamda = $6/f;			# rake
# magnitude
    ml = $7;
# moment
    M0 = 10**(3./2.*ml + 9.1);	# Nm
# trigonometry
    sin_phi   = sin(phi);  
    cos_phi   = cos(phi);  
    sin_2phi  = sin(2.0*phi);  
    cos_2phi  = cos(2.0*phi);
    sin2_phi  = sin_phi**2;
    cos2_phi  = cos_phi**2;
    sin_delta = sin(delta);
    cos_delta = cos(delta);
    sin_2delta=sin(2.0*delta);
    cos_2delta=cos(2.0*delta);
    sin_lamda   = sin(lamda);
    cos_lamda = cos(lamda);
#
# moment components
#  +x as North, +y as East, and +z as down. Sinc
    Mxx = -(sin_delta  * cos_lamda * sin_2phi +     sin_2delta * sin_lamda * sin2_phi);
    
    Mxy =  (sin_delta  * cos_lamda * cos_2phi + 0.5*sin_2delta * sin_lamda * sin_2phi);
    
    Mxz = -(cos_delta  * cos_lamda * cos_phi  +     cos_2delta * sin_lamda * sin_phi);
    
    Myy =  (sin_delta  * cos_lamda * sin_2phi -     sin_2delta * sin_lamda * cos2_phi);
    
    Myz = -(cos_delta  * cos_lamda * sin_phi  -     cos_2delta * sin_lamda * cos_phi);
    
    Mzz =  (sin_2delta * sin_lamda);
#
# convert to regular spherical r,theta,phi system
#
    Mrr =  Mzz;
    Mrt =  Mxz;
    Mrp = -Myz;
    Mtt =  Mxx;
    Mtp = -Mxy;
    Mpp =  Myy;

    if(cmt_style){
	# might have to multiply with sqrt(2) 
	# take log10 and 
	mexp = int(log(M0)*log10_fac+0.5);
	fac = 10**mexp; 
	M0 /= fac;
	# add 7 for going from Nm to dyn.cm
	mexp += 7;
	# X, Y, depth, mrr, mtt, mff, mrt, mrf, mtf, exp, newX, newY, event_title
	printf("%11g %11g %11g\t %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %i\t",
	       lon,lat,depth, Mrr*M0,Mtt*M0,Mpp*M0,Mrt*M0,Mrp*M0,Mtp*M0,mexp);
    }else{
#
# output: lon lat depth Mrr Mrt Mrp Mtt Mtp Mpp 
#
	printf("%11g %11g %11g\t %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\t",
	       lon,lat,depth,Mrr*M0,Mrp*M0,Mrt*M0,Mtt*M0,Mtp*M0,Mpp*M0);
    }
    for(i=8;i<=NF;i++)
	printf("%s ",$i);
    printf("\n");
  }

}
