#
# converts data in lon lat z v_x v_y v_z
# format from cartesian vector components to polar components
# input is 
#
# lon lat z v_x v_y v_z ....
#
# output is 
#
# lon lat z v_r v_theta v_phi ....
#
# there's also a C version of this in progs/src/misc
#
# $Id: lonlatcv2pv.awk,v 1.2 2005/02/23 01:42:59 becker Exp $
#
BEGIN{
   f=0.017453292519943296;
}
{
  if((substr($1,1,1)!="#") && (NF>=5)){
# coords
    phi=$1*f;
    theta=(90.0-$2)*f;
# polar components
    vx=$4;
    vy=$5;
    vz=$6;
    
# base vecs
    ct=cos(theta);cp=cos(phi);
    st=sin(theta);sp=sin(phi);
    # r base vec
    polar_base_x[1]= st * cp;
    polar_base_y[1]= st * sp;
    polar_base_z[1]= ct;
    # theta base vec
    polar_base_x[2]= ct * cp;
    polar_base_y[2]= ct * sp;
    polar_base_z[2]= -st;
    # phi base vec
    polar_base_x[3]= -sp;
    polar_base_y[3]=  cp;
    polar_base_z[3]= 0.0;
# convert vector
    for(i=1;i<=3;i++){
      polar_vec[i]  = polar_base_x[i] * vx;
      polar_vec[i] += polar_base_y[i] * vy;
      polar_vec[i] += polar_base_z[i] * vz;
    }
    printf("%20.15e %20.15e %20.15e %20.15e %20.15e %20.15e ",
	   $1,$2,$3,polar_vec[1],polar_vec[2],polar_vec[3]);
    for(i=7;i<=NF;i++)
      printf("%s ",$i);
    printf("\n");
  }
}
