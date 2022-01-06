#
# converts data in lon lat v_r v_theta v_phi [ ... ]
#
# THERE IS ALSO  A BINARY TO DO THAT
#
#
# format from polar vector components to cartesian components
#
# input is lon lat v_r v_theta v_phi [ ... ]
# 
# output is 
#
# lon lat v_x v_y v_z [ ... ]
#
# $Id: lonlatpv2cv.awk,v 1.2 2001/10/17 20:38:24 becker Exp $ 
#
BEGIN{
   f=0.017453292519943296;
}
{
  if((substr($1,1,1)!="#") && (NF>=5)){
# coords
    theta=(90.0-$2)*f;
    phi=$1*f;
# polar components
    vr=$3;
    vtheta=$4;
    vphi=$5;
#
# base vecs
    ct=cos(theta);
    cp=cos(phi);
    st=sin(theta);
    sp=sin(phi);
#
    polar_base_r[1]= st * cp;
    polar_base_r[2]= st * sp;
    polar_base_r[3]= ct;
#
    polar_base_theta[1]= ct * cp;
    polar_base_theta[2]= ct * sp;
    polar_base_theta[3]= -st;
#
    polar_base_phi[1]= -sp;
    polar_base_phi[2]=  cp;
    polar_base_phi[3]= 0.0;
# convert vector
    for(i=1;i<=3;i++){
      cart_vec[i]  = polar_base_r[i]    * vr;
      cart_vec[i] += polar_base_theta[i]* vtheta;
      cart_vec[i] += polar_base_phi[i]  * vphi;
    }
    printf("%20.15e %20.15e %20.15e %20.15e %20.15e ",
	   $1,$2,cart_vec[1],cart_vec[2],cart_vec[3]);
    for(i=6;i<=NF;i++)
      printf("%s ",$i);
    printf("\n");
  }

}
