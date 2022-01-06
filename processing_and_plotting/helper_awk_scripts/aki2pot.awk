#
# convert Aki type focal mechanisms into potency 
#
# using Yehuda's formula
#     input format 
#
#     X, Y, depth, strike, dip, rake, mag
#
# flags: 
#         normalize  0: use potencies
#                    1: use unity potency
#         weights:   0: no weighting
#                    1: multiply potency with 1-misfit/max_misfit
BEGIN{
  f = 57.2957795130823;
  max_misfit = 0.56;		# max misfit for weighting
#  low_cutoff = 0.9;		# lower mag cutoff, needs to be < 3.5
  low_cutoff = -1;		# lower mag cutoff, needs to be < 3.5
  high_cutoff = 5.0;		# high cutoff, needs to be > 3.5
}
{
  if((substr($1,1)!="#")&&(NF>=7)){
    lon = $1;if(lon<0)lon+=360.0;
    lat = $2;
    depth = $3;
# angles
    phi = $4/f;			        # strike
    delta = $5/f;			# dip
    lamda = $6/f;			# rake
# magnitude
    ml = $7;
# misfit
    if(NF>=10)
      misfit = $10;
    else
      misfit = 0.0;
# scalar potency ([m^3]?)
    if(ml < low_cutoff){
      print("skipping, magnitude too small:",ml) > "/dev/stderr";
      P0 = 0.0;
    }else if(ml < 3.5){		# lower scaling
      P0 = (normalize ? 1.0 : 10**(1.45 * ml - 5.69));
    }else if(ml <= high_cutoff){ # upper scaling
      P0 = (normalize ? 1.0 : 10**(1.08 * ml - 4.87));
    }else{
      print("skipping, magnitude too large:",ml) > "/dev/stderr";
      P0 = 0.0;
    }
    if(weights){
      if(misfit > max_misfit)
	print("misfit of ",misfit," is too large") > "/dev/stderr";
      weight = 1.0 - misfit/max_misfit;
    }else{
      weight = 1.0;
    }
    if(P0 != 0.0){
# trigonometry
      sin_phi     = sin(phi);  cos_phi   = cos(phi);  sin_2phi  =sin(2.0*phi);  cos_2phi  =cos(2.0*phi);
      sin2_phi=sin_phi**2;cos2_phi=cos_phi**2;
      sin_delta   = sin(delta);cos_delta = cos(delta);sin_2delta=sin(2.0*delta);cos_2delta=cos(2.0*delta);
      sin_lamda   = sin(lamda);cos_lamda = cos(lamda);
# potency components
#  +x as North, +y as East, and +z as down. Sinc
      Pxx = -P0 * (sin_delta  * cos_lamda * sin_2phi +     sin_2delta * sin_lamda * sin2_phi);

      Pxy =  P0 * (sin_delta  * cos_lamda * cos_2phi + 0.5*sin_2delta * sin_lamda * sin_2phi);

      Pxz = -P0 * (cos_delta  * cos_lamda * cos_phi  +     cos_2delta * sin_lamda * sin_phi);

      Pyy =  P0 * (sin_delta  * cos_lamda * sin_2phi -     sin_2delta * sin_lamda * cos2_phi);

      Pyz = -P0 * (cos_delta  * cos_lamda * sin_phi  -     cos_2delta * sin_lamda * cos_phi);

      Pzz =  P0 * (sin_2delta * sin_lamda);
#
# convert to regular spherical r,theta,phi system
#
      Prr =  Pzz;
      Prt =  Pxz;
      Prp = -Pyz;
      Ptt =  Pxx;
      Ptp = -Pxy;
      Ppp =  Pyy;
      
#    print(Pxx+Pxy+Pxz+Pyy+Pyz+Pzz) > "/dev/stderr";
#
# output: lon lat depth Prr Prt Prp Ptt Ptp Ppp 
#
      printf("%11g %11g %11g\t %.4e %.4e %.4e %.4e %.4e %.4e\t%g ",
	     lon,lat,depth,Prr,Prp,Prt,Ptt,Ptp,Ppp,ml);
      if(weights)
	printf("%.3f ",weight);
      for(i=8;i<=NF;i++)
	printf("%s ",$i);
      printf("\n");
    }
  }

}
