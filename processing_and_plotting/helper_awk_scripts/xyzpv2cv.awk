#
# converts data in x y z v_r v_theta v_phi [ ... ]
#
# format from polar vector components to cartesian components
#
# input is x y z v_r v_theta v_phi [ ... ]
# 
# output is 
#
# x y z v_x v_y v_z [ ... ]
#
BEGIN{
}
{
  if((substr($1,1,1)!="#") && (NF>=5)){
    x=$1;
    y=$2;
    z=$3;
# convert cartesian to spherical coordinates    
    tmp1  = x*x  + y*y;
    tmp2  = tmp1 + z*z;
    r = sqrt(tmp2);
    theta = atan2(sqrt(tmp1),z);
    phi   = atan2(y,x);
# polar components
    vr     = $4;
    vtheta = $5;
    vphi   = $6;
# spherical base vectors in cartesian frame
    ct=cos(theta);cp=cos(phi);
    st=sin(theta);sp=sin(phi);
    polar_base_r[1]= st * cp;
    polar_base_r[2]= st * sp;
    polar_base_r[3]= ct;
    polar_base_theta[1]= ct * cp;
    polar_base_theta[2]= ct * sp;
    polar_base_theta[3]= -st;
    polar_base_phi[1]= -sp;
    polar_base_phi[2]=  cp;
    polar_base_phi[3]= 0.0;
# convert vector
    for(i=1;i<=3;i++){
      cart_vec[i]  = polar_base_r[i] * vr ;
      cart_vec[i] += polar_base_theta[i] * vtheta ;
      cart_vec[i] += polar_base_phi[i]   * vphi;
    }
#
# output of x y z vx vy vz
    printf("%20.15e %20.15e %20.15e %20.15e %20.15e  %20.15e ",
	   x,y,z,cart_vec[1],cart_vec[2],cart_vec[3]);
    for(i=7;i<=NF;i++)
      printf("%s ",$i);
    printf("\n");
  }

}
