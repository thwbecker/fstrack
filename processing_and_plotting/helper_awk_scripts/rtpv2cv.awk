#
# converts data in r theta phi v_r v_theta v_phi 
#
# format from polar vector components to cartesian components
# using Earth values
#
BEGIN{
}
{
# coords
  theta=$2;
  phi=$3;

# polar components
  vr=$4;
  vtheta=$5;
  vphi=$6;

# base vecs
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
    cart_vec[i]  = polar_base_r[i]    * vr;
    cart_vec[i] += polar_base_theta[i]* vtheta;
    cart_vec[i] += polar_base_phi[i]  * vphi;
  }
  print(cart_vec[1],cart_vec[2],cart_vec[3]);

}
