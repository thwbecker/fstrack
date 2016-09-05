/*
  
output routines

$Id: output.c,v 1.33 2007/02/28 19:12:44 becker Exp becker $
  
*/


#include "fstrack_flow.h"

//
//  writes ASCII format data from array x at layer layer (0..n[FSTRACK_R]-1) to file
//
void write_ascii_layer(COMP_PRECISION *x,int layer,
		       struct mod *model,
		       char *arrname,char *filename,
		       COMP_PRECISION scalefac)
{
  FILE *out;
  int i,j,os;
  COMP_PRECISION theta,phi;
  if((layer < 0) || (layer >= MDP->n[FSTRACK_R])){
    fprintf(stderr,"write_ascii_layer: error, array %s has only %i layers, %i was given\n",
	    arrname,MDP->n[FSTRACK_R],layer);
    exit(-1);
  }
  out=myopen(filename,"w","write_ascii_layer");
  fprintf(stderr,"write_ascii_layer: writing layer %i(%i) z: %g  of %s to file %s, scaling by %g\n",
	  layer+1,MDP->n[FSTRACK_R],ZDEPTH(MDP->rlevels[layer]),arrname,filename,scalefac);
  for(theta=MDP->dtheta*0.5,os=MDP->n[TPPROD]*layer,i=0;i<MDP->n[FSTRACK_THETA];
      i++,theta += MDP->dtheta,os += MDP->n[FSTRACK_PHI])
    for(phi=0.0,j=0;j<MDP->n[FSTRACK_PHI];j++,phi+=MDP->dphi)
      fprintf(out,"%12g %12g %12g\n",PHI2LONGITUDE(phi),THETA2LATITUDE(theta),
	      x[os+j]*scalefac);
  fclose(out);
}
void write_ascii_layer_v(VPREC *x,int layer,
			 struct mod *model,
			 char *arrname,char *filename,
			 COMP_PRECISION scalefac)
{
  FILE *out;
  int i,j,os;
  COMP_PRECISION theta,phi;
  if((layer < 0) || (layer >= MDP->n[FSTRACK_R])){
    fprintf(stderr,"write_ascii_layer: error, array %s has only %i layers, %i was given\n",
	    arrname,MDP->n[FSTRACK_R],layer);
    exit(-1);
  }
  out=myopen(filename,"w","write_ascii_layer");
  fprintf(stderr,"write_ascii_layer: writing layer %i(%i) z: %g  of %s to file %s, scaling by %g\n",
	  layer+1,MDP->n[FSTRACK_R],ZDEPTH(MDP->rlevels[layer]),arrname,filename,scalefac);
  for(theta=MDP->dtheta*0.5,os=MDP->n[TPPROD]*layer,i=0;i<MDP->n[FSTRACK_THETA];
      i++,theta += MDP->dtheta,os += MDP->n[FSTRACK_PHI])
    for(phi=0.0,j=0;j<MDP->n[FSTRACK_PHI];j++,phi+=MDP->dphi)
      fprintf(out,"%12g %12g %12g\n",PHI2LONGITUDE(phi),THETA2LATITUDE(theta),
	      x[os+j]*scalefac);
  fclose(out);
}

void write_tracer_history(struct trc *tracer, char *filename, 
			  int mode,my_boolean verbose,
			  struct mod *model)
{
  FILE *out;
  if(tracer->discarded){
    fprintf(stderr,"write_tracer_history: this tracer was discarded\n");
    return;
  }
  out=myopen(filename,"w","write_tracer_history");
  print_tracer_data(tracer,mode,out,verbose,filename,model);
  fclose(out);
}

void write_all_tracer_history(struct mod *model, char *filename,
			      my_boolean gmt_style)
{
  FILE *out;
  int i;
  out=myopen(filename,"w","write_all_tracer_history");
  for(i=0;i<model->ntracer;i++){
    if(model->tracer[i].discarded){
      fprintf(stderr,"write_all_tracer_history: tracer %i was discarded\n",i);
    }
    if(gmt_style){
      print_tracer_data((model->tracer+i),LOCATION_LLZ,out,FALSE,filename,model);
      fprintf(out,">\n");
    }else{
      print_tracer_data((model->tracer+i),LOCATION_XYZ,out,FALSE,filename,model);
      fprintf(out,"\n");
    }
  }
  fprintf(stderr,"write_all_tracer_history: written all tracer histories to %s in ",
	  filename);
  if(gmt_style)
    fprintf(stderr,"gmt style\n");
  else
    fprintf(stderr,"cartesian coordinates style\n");
  fclose(out);
}
/*

output of quantities for individual tracer at all states

*/
void print_tracer_data(struct trc *tracer, int mode, 
		       FILE *out, my_boolean verbose,
		       char *filename,struct mod *model)
{
  int i,j,k,m,n,os1,os2,nwork=12,nxny;
  COMP_PRECISION eval[3],evec[9],ecvec[9],s[9],xc[3],xp[3],Fc[3][3],
    dx[12],l[9],l2[9],Savr[36],afactor[4],bfactor[4],period,strain_rate,vorticity,
    azi[2],amp[2],fac_l,fac_n,kernel[5],ba[2],hf[2],gl[2],ca[2],cn[2],
    fac_a,fac_c,fac_f    ;
  my_boolean frozen;
  char secfn[300];
  FILE *secout;
  switch(mode){
  case LOCATION_TRTP:/* location in  
			
		     time r theta phi  
		     
		     format */
    for(i=0;i<tracer->nstate;i++)
      fprintf(out,"%11g %13.6e %13.6e %13.6e\n",tracer->state[i].t,
	      tracer->state[i].x[FSTRACK_R],tracer->state[i].x[FSTRACK_THETA],
	      tracer->state[i].x[FSTRACK_PHI]);
    if(verbose)
      fprintf(stderr,"print_tracer_data: t r theta phi (%i steps) written to %s\n",
	      tracer->nstate,filename);
    break;
  case LOCATION_XYZ:/* location in 
			
		    x y z
			
		    (cartesian coordinates) format */
    for(i=0;i<tracer->nstate;i++){
      rtp2xyz(tracer->state[i].x,xc); /* convert to cartesian coordinates */
      fprintf(out,"%11g %11g %11g\n",xc[FSTRACK_X],xc[FSTRACK_Y],xc[FSTRACK_Z]);
    }
    if(verbose)
      fprintf(stderr,"print_tracer_data: x y z (%i steps) written to %s\n",
	      tracer->nstate,filename);
    break;
  case LOCATION_LLZ:/* location in 
			
		    lon lat z  
			
		    format */
    for(i=0;i<tracer->nstate;i++){
      print_llz((tracer->state+i),out);
      fprintf(out,"\n");
    }
    if(verbose)
      fprintf(stderr,"print_tracer_data: lon lat z (%i steps) written to %s\n",
	      tracer->nstate,filename);
    break;
  case LOCATION_TLLZ:/* location in 
			
		     t lon lat z  
		     
		     format */
    for(i=0;i<tracer->nstate;i++)
      print_tllz((tracer->state+i),out);
    if(verbose)
      fprintf(stderr,"print_tracer_data: t lon lat z (%i steps) written to %s\n",
	      tracer->nstate,filename);
    break;
  case STRAIN_COMP:/* 

		   left stretch tensor strain in 
		   
		   time L upper right half matrix format 
		   
		   ie.
		   
		   time L_rr L_rt L_rp L_tt L_tp L_pp
		   
		   format, spherical system
		   
		   */
    for(i=0;i < tracer->nstate;i++){
      /* 
	 convert the strain matrix from Cartesian to spherical 
      */
      cart_to_polar_mat_at_r((tracer->state[i].x+3),s,tracer->state[i].x);
      calc_l2_left_stretch(s,l2);calc_sqrt_sym_matrix(l2,l); 
      fprintf(out,"%11g %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e\n",
	      tracer->state[i].t,l[RR],l[RT],l[RP],l[TT],l[TP],l[PP]);
    }
    if(verbose)
      fprintf(stderr,"print_tracer_data: t L_ij strains (%i steps) written to %s\n",
	      tracer->nstate,filename);
    break;
  case XYZ_STRAIN_CART_COMP:/* 
			       tensorial strain from L (left-stretch)
			       converted to cartesian time L upper
			       right half matrix format
			       
			       x y z time L_xx L_xy L_xz L_yy L_yz L_zz
		      
			       format, Cartesian system
			    */
    for(i=0;i < tracer->nstate;i++){
      /* get L = sqrt(L^2), l2 is already in Cartesian system */
      calc_sqrt_sym_matrix(tracer->state[i].left_stretch,l); 
      symmat9to3x3(l,Fc);
      rtp2xyz(tracer->state[i].x,xc); /* location vector in cartesian */
      fprintf(out,"%15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n",
	      xc[FSTRACK_X],xc[FSTRACK_Y],xc[FSTRACK_Z],tracer->state[i].t,
	      Fc[FSTRACK_X][FSTRACK_X],Fc[FSTRACK_X][FSTRACK_Y],Fc[FSTRACK_X][FSTRACK_Z],
	      Fc[FSTRACK_Y][FSTRACK_Y],Fc[FSTRACK_Y][FSTRACK_Z],Fc[FSTRACK_Z][FSTRACK_Z]);
    }
    if(verbose)
      fprintf(stderr,"print_tracer_data: x y z t L_xx L_xy L_xz L_yy L_yz L_zz (%i steps)  written to %s\n",
	      tracer->nstate,filename);
    break;
  case STRAIN_EVAL:/* 
		      tensorial strain from L^2 in
		       
		      time sqrt(eigenvalue(s)) eigenvectors format, ie.
		       
		      time lon lat z e1 e2 e3 e1r e1t e1p e2r e2t e2p e3r e3t e3p
		      
		      where e1 > e2 > e3

		      spherical system
		   */
    for(i=0;i < tracer->nstate;i++){
      /* 
	 compute eigensystem of strain which is given in Cartesian 
      */
      calc_eigensystem_sym(tracer->state[i].left_stretch,eval,ecvec,TRUE);
      if(!MDP->rotate_grain_coord_sys){
	if(i == 0)
	  fprintf(stderr,"output: WARNING: seval is not rotated with xp, global Cartesian\n");
	for(j=0;j<3;j++)
	  a_equals_b_vector3d((evec+j*3),(ecvec+j*3));
      }else{
	/* 
	   convert all three vectors from Cartesian to spherical system 
	*/
	for(j=0;j<3;j++)
	  cart_vec2polar_vec_at_xp((ecvec+j*3),(evec+j*3),tracer->state[i].x);
      }
      /* 
	 output of eigenvectors in spherical system
      */
      print_tllz((tracer->state+i),out);
      fprintf(out,"%13.6e %13.6e %13.6e %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",
	      save_sqrt(eval[2]),save_sqrt(eval[1]),save_sqrt(eval[0]),
	      evec[6+FSTRACK_R],evec[6+FSTRACK_THETA],evec[6+FSTRACK_PHI],
	      evec[3+FSTRACK_R],evec[3+FSTRACK_THETA],evec[3+FSTRACK_PHI],
	      evec[  FSTRACK_R],evec[  FSTRACK_THETA],evec[  FSTRACK_PHI]);
    }
    if(verbose)
      fprintf(stderr,"print_tracer_data: t lon lat z sqrt(EV(L2)) e1 e2 e3 (%i steps, rtp)  written to %s\n",
	      tracer->nstate,filename);
    break;
  case STRAIN_EVAL_CART:/* 

			tensorial strain from L^2
		      
			time sqrt(eigenvalue(s)) eigenvectors format,
		       
			time x y z e1 e2 e3 e1x e1y e1z e2x e2y e2z e3x e3y e3z
		      
			converted to cartesian system

			where e1 > e2 > e3
			   
			*/
    for(i=0;i<tracer->nstate;i++){
      /* compute the eigensystem of the strain, which is in Cartesian */
      calc_eigensystem_sym(tracer->state[i].left_stretch,eval,ecvec,TRUE);
      /* convert location to cartesian */
      rtp2xyz(tracer->state[i].x,xc);// location in cartesian 
      /* 
	 output of eigensystem in Cartesian
      */
      fprintf(out,"%11g %11g %11g %11g %13.6e %13.6e %13.6e %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",
	      tracer->state[i].t,xc[FSTRACK_X],xc[FSTRACK_Y],xc[FSTRACK_Z],
	      save_sqrt(eval[2]),save_sqrt(eval[1]),save_sqrt(eval[0]),
	      ecvec[6+FSTRACK_X],ecvec[6+FSTRACK_Y],ecvec[6+FSTRACK_Z],
	      ecvec[3+FSTRACK_X],ecvec[3+FSTRACK_Y],ecvec[3+FSTRACK_Z],
	      ecvec[  FSTRACK_X],ecvec[  FSTRACK_Y],ecvec[  FSTRACK_Z]);
    }
    if(verbose)
      fprintf(stderr,"print_tracer_data: t lon lat z sqrt(EV(L2)) e1 e2 e3 (%i steps, cartesian)  written to %s\n",
	      tracer->nstate,filename);
    break;
  case CAUCHY_STRAIN_EVAL:/* 
			     Cauchy tensor strain in

			     time 1/sqrt(eigenvalue(B^1)) eigenvectors format, ie.
			     
			     time lon lat z e3 e2 e1 e3r e3t e3p e2r e2t e2p e1r e1t e1p
			     
			     where e1 > e2 > e3,

			     spherical system
			  */
    for(i=0;i < tracer->nstate;i++){
      calc_cauchy_strain((tracer->state[i].x+3),s);
      /* compute eigenvectors of Cauchy strain */
      calc_eigensystem_sym(s,eval,ecvec,TRUE);
      /* convert to spherical system */
      for(j=0;j<3;j++)
	cart_vec2polar_vec_at_xp((ecvec+j*3),(evec+j*3),tracer->state[i].x);
      print_tllz((tracer->state+i),out);
      /* output of eigensystem in spherical coordinates */
      fprintf(out,"%13.6e %13.6e %13.6e %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",
	      1./sqrt(eval[0]),1./sqrt(eval[1]),1./sqrt(eval[2]),
	      evec[  FSTRACK_R],evec[  FSTRACK_THETA],evec[  FSTRACK_PHI],
	      evec[3+FSTRACK_R],evec[3+FSTRACK_THETA],evec[3+FSTRACK_PHI],
	      evec[6+FSTRACK_R],evec[6+FSTRACK_THETA],evec[6+FSTRACK_PHI]);
    }
    if(verbose)
      fprintf(stderr,"print_tracer_data: t lon lat z 1/sqrt(EV(B-1)) e3 e2 e1 (%i steps)  written to %s\n",
	      tracer->nstate,filename);
    break;
  case DEFORMATION_TENSOR:
    /*
      
    output of F deformation tensor (non-symmetric) in 
    
    time F_11 F_12 F_13 ... 
    
    format, spherical system
    
    */
    for(i=0;i<tracer->nstate;i++){
      fprintf(out,"%11g ",tracer->state[i].t);
      /* 
	 convert the deformation tensor from Cartesian to spherical
	 system
      */
      cart_to_polar_mat_at_r((tracer->state[i].x+3),s,tracer->state[i].x);
      for(j=0;j<9;j++)
	fprintf(out,"%11.4e ",s[j]);
      fprintf(out,"\n");
    }
    if(verbose)
      fprintf(stderr,"print_tracer_data: t F_ij deformation (%i steps) for tracer written to %s\n",
	      tracer->nstate,filename);
    break;
  case ISA_AXES:
    /* 
       
    compute the ISA axis for this tracer at each timestep format:

    time ISA_r ISA_t ISA_p ISA_exists PI_par

    */
    for(i=0;i<tracer->nstate;i++){
      fprintf(out,"%11g ",tracer->state[i].t);
      /* convert the ISA from Cartesian to spherical */
      cart_vec2polar_vec_at_xp(tracer->state[i].isa,evec,tracer->state[i].x);
      fprintf(out,"%11.4e %11.4e %11.4e %1i %11g\n",
	      evec[FSTRACK_R],evec[FSTRACK_THETA],evec[FSTRACK_PHI],tracer->state[i].isa_exists,
	      tracer->state[i].pipar);
    }
    if(verbose)
      fprintf(stderr,"print_tracer_data: t ISA_r ISA_t ISA_p Pi_par (%i steps) for tracer written to %s\n",
	      tracer->nstate,filename);

    break;
  case LOCATION_VGM_RTPTVGM:	/* 
				   output of r theta phi t VGM where
				   VGM is the velocity gradient matrix
				   as given in C format by
				   fse_ellderivs in spherical system
				   
				*/
    for(i=0;i < tracer->nstate;i++){
      fprintf(out,"%16.8e %16.8e %16.8e %16.8e ",
	      tracer->state[i].x[FSTRACK_R],
	      tracer->state[i].x[FSTRACK_THETA],
	      tracer->state[i].x[FSTRACK_PHI],
	      tracer->state[i].t);
      /* 
	 obtain the velocity gradient matrix in spherical system
      */
      fse_derivs_wrapper(tracer->state[i].t,tracer->state[i].x,
			 dx,nwork,MDP,MDP->v_loc,MDP->vgm,FALSE,TRUE,
			 &frozen,MDP->strain_fraction_from_gamma); 
      if(frozen)
	fprintf(stderr,"print_tracer_data: WARNING: fse_derivs returned frozen\n");
      if(MDP->strain_fraction_from_gamma)
	fprintf(stderr,"print_tracer_data: WARNING: fse_deriv using only alpha strain fraction\n");
      /*     
	     we want output to be sorted like

	     d_r(vr) d_t(vr) d_p(vr) ...
	     d_r(vt) d_t(vt) d_p(vt) ...
	     d_r(vp) d_t(vp) d_p(vp) 
	     
      */
      fprintf(out,"%16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n",
	      MDP->vgm[RR],MDP->vgm[RT],MDP->vgm[RP],
	      MDP->vgm[TR],MDP->vgm[TT],MDP->vgm[TP],
	      MDP->vgm[PR],MDP->vgm[PT],MDP->vgm[PP]);
    }
    if(verbose)
      fprintf(stderr,"print_tracer_data: r theta phi t VGM_ij  (%i steps) for tracer written to %s\n",
	      tracer->nstate,filename);
    break;
  case LOCATION_VGM_RTPTVGM_CART:	/* output of r theta phi t VGM where
					   VGM is the velocity gradient matrix
					   as given in C format by
					   fse_ellderivs in a Cartesian system
				   
					*/
    for(i=0;i < tracer->nstate;i++){
      fprintf(out,"%16.8e %16.8e %16.8e %16.8e ",
	      tracer->state[i].x[FSTRACK_R],
	      tracer->state[i].x[FSTRACK_THETA],
	      tracer->state[i].x[FSTRACK_PHI],
	      tracer->state[i].t);
      /* 
	 obtain the velocity gradient matrix in Cartesian system
      */
      fse_derivs_wrapper(tracer->state[i].t,tracer->state[i].x,
			 dx,nwork,MDP,MDP->v_loc,MDP->vgm,TRUE,TRUE,
			 &frozen,MDP->strain_fraction_from_gamma); 
      if(MDP->strain_fraction_from_gamma)
	fprintf(stderr,"print_tracer_data: WARNING: fse_deriv using only alpha strain fraction\n");
      if(frozen)
	fprintf(stderr,"print_tracer_data: WARNING: fse_derivs returned frozen\n");
      /*     
	     we want output to be sorted like

	     d_x(vx) d_y(vx) d_z(vx) ...
	     d_x(vy) d_y(vy) d_z(vy) ...
	     d_x(vz) d_y(vz) d_z(vz) 
	     
      */
      fprintf(out,"%16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n",
	      MDP->vgm[XX],MDP->vgm[XY],MDP->vgm[XZ],
	      MDP->vgm[YX],MDP->vgm[YY],MDP->vgm[YZ],
	      MDP->vgm[ZX],MDP->vgm[ZY],MDP->vgm[ZZ]);
    }
    if(verbose)
      fprintf(stderr,"print_tracer_data: r theta phi t VGM_ij  (%i steps) for tracer written to %s\n",
	      tracer->nstate,filename);
    break;
  case STRAINRATE_VORTICITY:	/* output of strain-rate and vorticity from Cartesian 
				   velocity matrix
				*/
    for(i=0;i < tracer->nstate;i++){
      /* 
	 obtain the velocity gradient matrix in Cartesian system
      */
      fse_derivs_wrapper(tracer->state[i].t,tracer->state[i].x,
			 dx,nwork,MDP,MDP->v_loc,MDP->vgm,TRUE,TRUE,
			 &frozen,MDP->strain_fraction_from_gamma); 
      if(MDP->strain_fraction_from_gamma)
	fprintf(stderr,"print_tracer_data: WARNING: fse_deriv using only alpha strain fraction\n");
      if(frozen)
	fprintf(stderr,"print_tracer_data: WARNING: fse_derivs returned frozen\n");
      calc_strainrate_and_vorticity(MDP->vgm,&strain_rate,&vorticity,TRUE);
      fprintf(out,"%16.8e %16.8e\n",strain_rate,vorticity);
    }
    if(verbose)
      fprintf(stderr,"print_tracer_data: strain rate and vorticity (%i steps) for tracer written to %s\n",
	      tracer->nstate,filename);
    break;
  case LOCATION_VGM_XYZTVGM_CART:	/* output of x y z t VGM where
					   VGM is the velocity
					   gradient matrix as given in
					   C format by fse_ellderivs
					   in a Cartesian system
					   
					*/
    for(i=0;i < tracer->nstate;i++){
      rtp2xyz(tracer->state[i].x,xc); /* convert to cartesian */
      fprintf(out,"%16.8e %16.8e %16.8e %16.8e ",
	      xc[FSTRACK_X],xc[FSTRACK_Y],xc[FSTRACK_Z],tracer->state[i].t);
      /* 
	 obtain the velocity gradient matrix in Cartesian system
      */
      fse_derivs_wrapper(tracer->state[i].t,tracer->state[i].x,
			 dx,nwork,MDP,MDP->v_loc,MDP->vgm,TRUE,TRUE,
			 &frozen,MDP->strain_fraction_from_gamma); 
      if(frozen)
	fprintf(stderr,"print_tracer_data: WARNING: fse_derivs returned frozen\n");
      if(MDP->strain_fraction_from_gamma)
	fprintf(stderr,"print_tracer_data: WARNING: fse_deriv using only alpha strain fraction\n");

      /*     
	     we want output to be sorted like
	     
	     d_x(vx) d_y(vx) d_z(vx) ...
	     d_x(vy) d_y(vy) d_z(vy) ...
	     d_x(vz) d_y(vz) d_z(vz) 
	     
      */
      fprintf(out,"%16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e\n",
	      MDP->vgm[XX],MDP->vgm[XY],MDP->vgm[XZ],
	      MDP->vgm[YX],MDP->vgm[YY],MDP->vgm[YZ],
	      MDP->vgm[ZX],MDP->vgm[ZY],MDP->vgm[ZZ]);
    }
    if(verbose)
      fprintf(stderr,"print_tracer_data: x y z t VGM_ij  (%i steps) for tracer written to %s\n",
	      tracer->nstate,filename);
    break;
  case TEXTURE_SAV:
    /*
      
    output of texture derived symmetrix stiffness matrix with 36
    components in UPPER TRIANGULAR FORMAT
		   
    lon lat depth SAV_11 SAV_12 SAV_13 ... SAV_16 \
    SAV_22 SAV_23 ... SAV_26 \
    ... 
    SAV_66
                             
    
    format. 

    THE MATRIX WILL BE ROTATED, such that locally x points south, y
    points east, and z ppoints up
    
    */
    for(i=0;i< tracer->nstate;i++){

      /* 
	 rotate the stiffness matrix from the gloval Caretsian to a
	 local Caretsian frame, where z is up, x is south, and y is
	 East
      */
      print_llz((tracer->state+i),out);
      if(!MDP->rotate_grain_coord_sys){
	if(i == 0)
	  fprintf(stderr,"output: WARNING: SAV is not rotated with xp, global Cartesian\n");
	print_sym_6by6(tracer->state[i].Sav,out);
      }else{
	/* 
	   rotate into geographic local frame, x = S, y = E, z = UP, the regular cartesian system
	*/
	rotate_mat_to_reg_cart_system(tracer->state[i].Sav, tracer->state[i].x,Savr);
	/* print */
	print_sym_6by6(Savr,out);
      }
    }
    if(verbose)
      fprintf(stderr,"print_tracer_data: t Sav_ij stiffness (%i steps) in %s system for tracer written to %s\n",
	      tracer->nstate,"Cartesian",filename);
    break;
  case TEXTURE_TI:
    /*
      
    output of texture derived best-fit hexagoal anisotropy axis and 
    anisotropy percent
		   
    time t_r t_t t_phi vti_ani_perc iso_frac hex_frac tet_frac orth_frac mon_fra tri_frac
    
    format (spherical system)
    
    */
    for(i=0;i< tracer->nstate;i++){
      fprintf(out,"%11g\t",tracer->state[i].t);
      process_ti_tens_frac(tracer->state[i].x,
			   tracer->state[i].Sav,
			   MDP->rotate_grain_coord_sys,out);
      fprintf(out,"\n");
      /* 
	 note that this would give the same result: in the rotated
	 reference frame z is up, x is south and y is east
      */
      /* coord_convetion = DREX_REG_CART_CONVENTION; */
      /*  fprintf(stderr,"%g %g %g\n",evec[FSTRACK_R],evec[FSTRACK_THETA],evec[FSTRACK_PHI]); */
      /*       drex_calc_euler_rad_from_xp(tracer->state[i].x,&alpha,&beta,&gamma,&coord_convention); */
      /*       drex_rotate_6x6_rad_ftrn(tracer->state[i].Sav,Savr,&alpha,&beta,&gamma); */
      /*       drex_decsym(Savr,amp,......ecvec); */
      /*       fprintf(stderr,"%g %g %g\n",ecvec[FSTRACK_Z],ecvec[FSTRACK_X],ecvec[FSTRACK_Y]); */
    }
    if(verbose)
      fprintf(stderr,"print_tracer_data: t t_i perc best fit TI (%i steps) for tracer written to %s\n",
	      tracer->nstate,filename);
    break;
  case TEXTURE_ODF:
    /* 
       output of the orientation density functions for olivine and
       enstatite. the major file will only hold the times, size, nx, ny, and
       the name of the subordinate files with the data

    */
    nxny = DREX->np[0] * DREX->np[1];
    system("mkdir poles/ 2> /dev/null");
    for(i=0;i< tracer->nstate;i++){
      if(tracer->state[i].npole){ 
	/* only those states that have saved pole figure densitites */
	/* 
	   filename for all binary output for this state 
	*/
	sprintf(secfn,"poles/%s.%i",filename,i);
	/* information file:
	   
	time isize nx ny lon lat z grains_in_ENU_system name
	
	if grains_in_ENU_system is TRUE, ODF will have been rotated
	such that everything is in a geographic system with East,
	North, and Up the local x, y, and z axis. Else, ODF is in 
	a global Cartesian system

	*/
	lonlatz_from_xp(tracer->state[i].x,xp);
	fprintf(out,"%11g\t%i\t%i %i\t%12g %12g %12g\t%i\t%s\n",
		tracer->state[i].t,DREX->size,
		DREX->np[0],DREX->np[1],xp[FSTRACK_X],xp[FSTRACK_Y],xp[FSTRACK_Z],
		model->dp->rotate_grain_coord_sys,secfn);
	if(tracer->state[i].npole != nxny * 3 * 2){
	  fprintf(stderr,"print_tracer_data: error: pole n: %i for state %i, should be %i\n",
		  tracer->state[i].npole,i+1,nxny * 3 * 2);
	  exit(-1);
	}
	for(j=os1=0;j < 2;j++,os1 += 3*nxny){		
	  /* olivine / enstatite loop */
	  for(k=os2=0;k < 3;k++,os2 += nxny){ 
	    /* axis loop */
	    if(j==0)		/* olivine */
	      sprintf(secfn,"poles/%s.%i.o.%1i",filename,i,k+1);
	    else			/* enstatite */
	      sprintf(secfn,"poles/%s.%i.e.%1i",filename,i,k+1);
	    secout = myopen(secfn,"w","print_tracer_data");
	    for(m=0;m  < DREX->np[1];m++){/* lat loop */
	      for(n=0;n < DREX->np[0];n++){ /* lon loop */
		//		fprintf(secout,"%.3e\n",
		//tracer->state[i].pdens[os1+os2+m*DREX->np[0]+n]);
		/* binary output */
		fwrite((tracer->state[i].pdens+os1+os2+m*DREX->np[0]+n),
		       sizeof(COMP_PRECISION),1,secout);
	      }
	    }
	    fclose(secout);
	  } /* end axis loop */
	}	/* end oli/ens loop */
      }
    }
    if(verbose)
      fprintf(stderr,"print_tracer_data: ODF parameters (%i steps) in %s system for tracer written to %s\n",
	      tracer->nstate,"Cartesian",filename);
    break;
  case TEXTURE_STATS:
    /* 
       
    print statistics on the olivine axes orientations, format: unweighted and weighted

    lon lat z time      a_r a_t a_p 1-<a.<a>>   aw_r aw_t aw_p 1-<aw.<aw>>  

    i.e. those are reordered from there internal E-N-U format

    */
    for(i=0;i< tracer->nstate;i++){
      print_llzt((tracer->state+i),out);
      for(j=0;j < 2;j++)		/* loop through olivine [100]
					   unweighted and weighted
					   mean axes and mean
					   deviations from it */ 
	fprintf(out,"%10.6f %10.6f %10.6f\t%12.5e\t",
		tracer->state[i].olivine_axes_stats[j*3+FSTRACK_X],
		tracer->state[i].olivine_axes_stats[j*3+FSTRACK_Y],
		tracer->state[i].olivine_axes_stats[j*3+FSTRACK_Z],
		tracer->state[i].olivine_axes_stats[6+j]);
      fprintf(out,"\t%12.5e\n",
	      tracer->state[i].olivine_axes_stats[8]); /* J index */
    }
    if(verbose)
      fprintf(stderr,"print_tracer_data: olivine [100] stats (%i steps) for tracer written to %s\n",
	      tracer->nstate,filename);
    break;
  case TEXTURE_RAYLEIGH_RPHI:
    if(!model->sw_sens_init){
      fprintf(stderr,"print_tracer_data: error: sensitivity kernels not initialized\n");
      exit(-1);
    }  
    if(!MDP->rotate_grain_coord_sys){
      fprintf(stderr,"output: rphi terms will not be printed, since no rotation was selected\n");
    }else{
      /* 
	 
      
      compute the 2phi and 4phi terms for Rayleigh waves of certain periods
      format:

      lon lat z time l n period azi amp period azi amp period azi amp
      
      where l and n are the vsv vsh factors l=\rho v_sv^2 and n=\rho
      v_sh^2
      
      */
      for(i=0;i< tracer->nstate;i++){
	print_llzt((tracer->state+i),out);
	for(j=0;j < N_SW_SENS;j++){ /* loop through periods */
	  period = model->swsens[j].p;
	  /* 
	     
	  this computes the 2phi and 4phi factors at the tracer depth
	  in the rotated coordinate frame (z=up, x=down, y=east) using
	  sensitivity kernels. phi is azimuth clockwise from north
	  a[0] = 2 phi cos;a[1] = 2 phi sin
	  a[2] = 4 phi cos;a[3] = 4 phi sin
	  
	  call with DREX_SW_CART_CONVENTION to rotate directly into
	  the SW system

	  */
	  compute_phi_terms_from_Sav(tracer->state[i].Sav,period,
				     tracer->state[i].x,
				     afactor,bfactor,&fac_l,&fac_n,
				     &fac_a,&fac_c,&fac_f,
				     kernel,(kernel+1),(kernel+3),(kernel+4),(kernel+5),
				     ba,hf,gl,ca,cn,
				     model,DREX_SW_CART_CONVENTION);
	  /* compute the azimuths (degrees CW from North) and
	     aplitudes */
	  compute_azideg_amp_from_afactors(afactor,azi,amp);
	  if(j==0)/* 
		     vsv/vsh factors
		     
		     l=\rho v_sv^2
		     n=\rho v_sh^2
		     
		  */
	    fprintf(out,"%.4e %.4e\t",fac_l,fac_n);
	  fprintf(out,"%5.1f\t%7.2f %7.2e %7.2f %7.2e\t",
		  period,azi[0],amp[0],azi[1],amp[1]);
	}
	fprintf(out,"\n");
      }
      if(verbose)
	fprintf(stderr,"print_tracer_data: Rayleigh 2/4phi (%i steps) for tracer written to %s\n",
		tracer->nstate,filename);
    }
    break;
  default:
    fprintf(stderr,"print_tracer_data: mode %i not defined\n",
	    mode);
    exit(-1);
    break;
  }
}
/*
  
write the states of tracers first_tracer through last_tracer to file
first_tracer and last_tracer run from 0 .. model->ntracer-1
  
*/
void write_tracer_field(struct mod *model,int first_tracer,int last_tracer, 
			int nstate,char *filename,int mode,my_boolean verbose)
{
  FILE *out;
  int i,j,usestate;
  COMP_PRECISION eval[3],evec[9],ecvec[9],l[9],xp[3],fac_n,fac_l,
    fac_a,fac_c,fac_f,bfactor[4],
    s[9],l2[9],period,afactor[4],azi[2],amp[2],kernel[5],ba[2],hf[2],gl[2],ca[2],Savr[36],cn[2];
  char message[200];
  
  if((first_tracer < 0) || (last_tracer > model->ntracer-1)){
    fprintf(stderr,"write_tracer_field: index problem: %i through %i but only %i tracers in model\n",
	    first_tracer,last_tracer,model->ntracer);
    PEE("write_tracer_field: HINT: use 0 .. ntracer-1 convention");
  }
  if(nstate == -1)
    strcpy(message, "last state");
  else
    sprintf(message,"state  %3i",nstate);
  /* 
     open output file
  */
  out=myopen(filename,"w","write_tracer_field");
  if(mode == ALL_TRACER_LLZAT)// check if we have read in attributes
    if(model->tinitmode != SPOTTED_3D_WITH_ATTR){// if not, use normal LLZ output
      if(verbose)
	fprintf(stderr,"write_tracer_field: WARNING: no attributes read, reverting to simple LLZ output\n");
      mode = ALL_TRACER_LLZT;
    }
  if(verbose)
    switch(mode){// write message
    case ALL_TRACER_RTPT:
      fprintf(stderr,"write_tracer_field: locations from %i to %i at %s in r theta phi t  to %s\n",
	      first_tracer,last_tracer,message,filename);
      break;
    case ALL_TRACER_LLZT:
      fprintf(stderr,"write_tracer_field: locations from %i to %3i at %s in lon lat z t  to %s\n",
	      first_tracer,last_tracer,message,filename);
      break;
    case ALL_TRACER_LLZAT:
      fprintf(stderr,"write_tracer_field: pos. + attr from %i to %3i at %s in lon lat z a_1 ... t  to %s\n",
	      first_tracer,last_tracer,message,filename);
      break;
    case ALL_TRACER_STRAIN_COMP:
      fprintf(stderr,"write_tracer_field: strain from %i to %i at %s in lon lat z L_ij Dt dz dx to %s\n",
	      first_tracer,last_tracer,message,filename);
      break;
    case ALL_TRACER_ISA:
      fprintf(stderr,"write_tracer_field: ISA from %i to %i at %s in lon lat z ISA_r ISA_t ISA_p isa_exists pi_parameter %s\n",
	      first_tracer,last_tracer,message,filename);
      break;
    case ALL_TRACER_STRAIN_EIGEN:
      fprintf(stderr,"write_tracer_field: strain from %i to %i at %s in lon lat z sqrt(EV(L2)) e1 e2 e3 e1r e1t e1p ... age  to %s\n",
	      first_tracer,last_tracer,message,filename);
      break;
    case ALL_TRACER_STRAIN_EIGEN_ANGLE:
      fprintf(stderr,"write_tracer_field: strain from %i to %i at %s in lon lat z sqrt(EV(L2)) e1 e2 e3 e1_azi e1_dip  ... age  to %s\n",
	      first_tracer,last_tracer,message,filename);
      break;
    case ALL_TRACER_DEFORMATION:
      fprintf(stderr,"write_tracer_field: deform from %i to %i at %s in lon lat z F_ij age  to %s\n",
	      first_tracer,last_tracer,message,filename);
      break;
    case ALL_TRACER_LYA:
      fprintf(stderr,"write_tracer_field: Lya from %i to %i at %s in lon lat z l1 l2 l3 t  to %s\n",
	      first_tracer,last_tracer,message,filename);
      break;
    case ALL_TRACER_TI:
      fprintf(stderr,"write_tracer_field: transv. isotr. from texture %i to %i at %s in lon lat z ar at ap perc_ani  to %s\n",
	      first_tracer,last_tracer,message,filename);
      break;
    case ALL_TRACER_RPHI:
      fprintf(stderr,"write_tracer_field: Rayleigh 2/4phi. from texture %i to %i at %s in lon lat z l n period azi amp period azi amp period azi amp  to %s\n",
	      first_tracer,last_tracer,message,filename);
      break;
    case ALL_TRACER_SAV:
      fprintf(stderr,"write_tracer_field: C(6,6) stiffness matrix. from texture %i to %i at %s in lon lat z C_ij  to %s\n",
	      first_tracer,last_tracer,message,filename);
      break;
    default:
      fprintf(stderr,"write_tracer_field: mode %i undefined\n",mode);
      exit(-1);
      break;
    }
  for(i=first_tracer;i <= last_tracer;i++){
    /* 

    start main tracer (i) loop
    
    */
    if(!model->tracer[i].discarded){
      if(nstate == -1)// use last state
	usestate = model->tracer[i].nstate-1;
      else 
	usestate = nstate;
      if(model->tracer[i].nstate < usestate){
	fprintf(stderr,"write_tracer_field: error: tracer %i has no state %i, maxstate %i\n",
		i,nstate,model->tracer[i].nstate);
	exit(-1);
      }
      switch(mode){
      case ALL_TRACER_RTPT:// r theta phi time
	fprintf(out,"%15.7e %15.7e %15.7e %g\n",
		model->tracer[i].state[usestate].x[FSTRACK_R],
		model->tracer[i].state[usestate].x[FSTRACK_THETA],
		model->tracer[i].state[usestate].x[FSTRACK_PHI],
		model->tracer[i].state[usestate].t);
	break;
      case ALL_TRACER_LLZT:// lon lat depth time
	print_llzt(&(model->tracer[i].state[usestate]),out);
	fprintf(out,"\n");
	break;
      case ALL_TRACER_LLZAT:// lon lat depth attribute(s) time
	print_llz(&(model->tracer[i].state[usestate]),out);
	for(j=0;j<NR_T_ATTR;j++)
	  fprintf(out,"%11g ",model->tracer[i].attr[j]);
	fprintf(out,"%11g\n",model->tracer[i].state[usestate].t);
	break;
      case ALL_TRACER_STRAIN_COMP:/*
				    location and left-stretch strain L
				    
				    lon lat depth L upper right half age dz dx
				    
				    spherical system

				  */


	/* 
	   convert the strain from Cartesian to spherical system
	*/
	cart_to_polar_mat_at_r((model->tracer[i].state[usestate].x+3),
			       s,model->tracer[i].state[usestate].x); 

	calc_l2_left_stretch(s,l2);calc_sqrt_sym_matrix(l2,l); 

	lonlatz_from_xp(model->tracer[i].state[usestate].x,xp);
	fprintf(out,"%11.4e %11.4e %11.4e %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e %10g %11g %11g\n",
		xp[FSTRACK_X],xp[FSTRACK_Y],xp[FSTRACK_Z],
		l[RR],l[RT],l[RP],l[TT],l[TP],l[PP],
		(model->tracer[i].state[usestate].t - model->tracer[i].state[0].t),
		ZDEPTH(model->tracer[i].state[usestate].x[FSTRACK_R])-
		ZDEPTH(model->tracer[i].state[0].x[FSTRACK_R]),
		dist_on_sphere(model->tracer[i].state[usestate].x,
			       model->tracer[i].state[0].x));
	break;
      case ALL_TRACER_ISA:
	/*
	 
	location and infinite strain axes which should have been precomputed
       
	lon lat depth isa_R isa_THETA isa_phi success_flag pi_par
       
	spherical system
	
	*/
	/* convert the ISA from Cartesian to spherical */
	cart_vec2polar_vec_at_xp(model->tracer[i].state[usestate].isa,evec,
				 model->tracer[i].state[usestate].x);
	print_llz(&(model->tracer[i].state[usestate]),out);
	fprintf(out,"%13.6e %13.6e %13.6e %1i %11g\n",
		evec[FSTRACK_R],evec[FSTRACK_THETA],evec[FSTRACK_PHI],
		model->tracer[i].state[usestate].isa_exists,
		model->tracer[i].state[usestate].pipar);
	break;
      case ALL_TRACER_STRAIN_EIGEN_ANGLE:
      case ALL_TRACER_STRAIN_EIGEN:/* 
				      location and eigenvectors of
				      left-stretch strain
				      
				      lon lat depth e1 e2 e3 e1r e1t \
				      e1p e2r ... age 
				      
				      (spherical system) or use angles
					 
				      lon lat depth e1 e2 e3  \
				      e1azi e1dip e2azi e2dip ... age 
					 
				   */
	/* left-stretch is computed in Cartesian system  */
	calc_eigensystem_sym(model->tracer[i].state[usestate].left_stretch,
			     eval,ecvec,TRUE);
	for(j=0;j<3;j++)	/* convert vectors to spherical system */
	  cart_vec2polar_vec_at_xp((ecvec+j*3),(evec+j*3),
				   model->tracer[i].state[usestate].x);
	print_llz(&(model->tracer[i].state[usestate]),out);
	if(mode == ALL_TRACER_STRAIN_EIGEN) /* normal vectors */
	  fprintf(out,"%13.6e %13.6e %13.6e %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %11g\n",
		  save_sqrt(eval[2]),save_sqrt(eval[1]),save_sqrt(eval[0]),
		  evec[6+FSTRACK_R],evec[6+FSTRACK_THETA],evec[6+FSTRACK_PHI],
		  evec[3+FSTRACK_R],evec[3+FSTRACK_THETA],evec[3+FSTRACK_PHI],
		  evec[  FSTRACK_R],evec[  FSTRACK_THETA],evec[  FSTRACK_PHI],
		  (model->tracer[i].state[usestate].t-model->tracer[i].state[0].t));
	else			/* angles */
	  fprintf(out,"%13.6e %13.6e %13.6e %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %11g\n",
		  save_sqrt(eval[2]),save_sqrt(eval[1]),save_sqrt(eval[0]),
		  azi_from_etep(evec[6+FSTRACK_THETA],evec[6+FSTRACK_PHI]),
		  dip_from_er_unityvec(evec[6+FSTRACK_R]),
		  azi_from_etep(evec[3+FSTRACK_THETA],evec[3+FSTRACK_PHI]),
		  dip_from_er_unityvec(evec[3+FSTRACK_R]),
		  azi_from_etep(evec[FSTRACK_THETA],evec[FSTRACK_PHI]),
		  dip_from_er_unityvec(evec[FSTRACK_R]),
		  (model->tracer[i].state[usestate].t-model->tracer[i].state[0].t));
	break;
      case ALL_TRACER_TI:/* 
			    
			 best fit symmetry axis of transverse isotropy hexagonal fit to 
			 stiffness matrix 
			 
			 lon lat depth ar at aphi percent_anisotrpy  iso_frac hex_frac tet_frac orth_frac mon_fra tri_frac

			 format (spherical system)

			 */

	print_llz(&(model->tracer[i].state[usestate]),out);
	process_ti_tens_frac(model->tracer[i].state[usestate].x,
			     model->tracer[i].state[usestate].Sav,
			     TRUE,out);	/* always rotate */
	fprintf(out,"\n");
	break;
      case ALL_TRACER_SAV:/* 

			  average stiffness matrix  in upper triangular format
			  
			  lon lat depth C_ij (j fast), i.e.
			  
			  SAV_11 SAV_12 SAV_13 ... SAV_16 \
			  SAV_22 SAV_23 ... SAV_26 \
			  ... 
			  SAV_66

			  rotate the stiffness matrix from the gloval Caretsian to a
			  local Cartesian frame, where z is up, x is south, and y is
			  East

			  */
	rotate_mat_to_reg_cart_system(model->tracer[i].state[usestate].Sav, 
				      model->tracer[i].state[usestate].x,Savr);
	/* print */
	print_llz(&(model->tracer[i].state[usestate]),out);
	print_sym_6by6(Savr,out);
	break;
      case ALL_TRACER_RPHI:/* 

			   vsv/vsh factors and 
			   
			   Rayleigh wave 2/4 phi terms
			   
			   format (spherical system)
			   
			   */
	if(!model->sw_sens_init){
	  fprintf(stderr,"print_tracer_data: error: sensitivity kernels not initialized\n");
	  exit(-1);
	}    
	print_llz(&(model->tracer[i].state[usestate]),out);
	for(j=0;j < N_SW_SENS;j++){ /* loop through periods */
	  period = model->swsens[j].p;
	  /* 
	     
	     rotate tensor into SW system and compute rphi terms
	     
	  */
	  compute_phi_terms_from_Sav(model->tracer[i].state[usestate].Sav,
				     period,model->tracer[i].state[usestate].x,
				     afactor,bfactor,&fac_l,&fac_n,
				     &fac_a,&fac_c,&fac_f,
				     kernel,(kernel+1),(kernel+3),(kernel+4),(kernel+5),
				     ba,hf,gl,ca,cn,
				     model,DREX_SW_CART_CONVENTION);
	  compute_azideg_amp_from_afactors(afactor,azi,amp);
	  if(j==0)		/* 
				   vsv/vsh factors

				   l=\rho v_sv^2
				   n=\rho v_sh^2

				*/
	    fprintf(out,"%.4e %.4e\t",fac_l,fac_n);
	  fprintf(out,"%5.1f\t%7.2f %7.2e %7.2f %7.2e\t",
		  period,azi[0],amp[0],azi[1],amp[1]);
	}
	fprintf(out,"\n");
	break;
      case ALL_TRACER_DEFORMATION:/*
				    location and deformation matrix
				    
				    lon lat depth F11 .. F33 age

				  */
	print_llz(&(model->tracer[i].state[usestate]),out);
	/* 
	   convert the deformation tensor from Cartesian to spherical
	   system
	*/
	cart_to_polar_mat_at_r((model->tracer[i].state[usestate].x+3),s,
			       model->tracer[i].state[usestate].x);
	for(j=0;j < 9;j++)
	  fprintf(out,"%13.6e ",s[j]);
	fprintf(out,"%g\n",(model->tracer[i].state[usestate].t - 
			    model->tracer[i].state[0].t));
	break;
      case ALL_TRACER_LYA:/*
			    location and lyapunov exponents

			    lon lat depth l1 l2 l3 t

			  */
#ifdef DEBUG
	if(fabs(model->tracer[i].state[usestate].t)<EPS_PREC){
	  fprintf(stderr,"error, tracer %i state %i has time 0\n",
		  i,usestate);
	  exit(-1);
	}
#endif
	print_llz(&(model->tracer[i].state[usestate]),out);
	for(j=0;j<3;j++)	/* those were stored in the ISA vector */
	  fprintf(out,"%13.6e ",model->tracer[i].state[usestate].isa[j]/
		  model->tracer[i].state[usestate].t);
	fprintf(out,"%g\n",
		(model->tracer[i].state[usestate].t-model->tracer[i].state[0].t));
	break;
      default:
	fprintf(stderr,"write_tracer_field: mode %i undefined\n",mode);
	exit(-1);
	break;
      }
    }
  }
  fclose(out);
}

/* 
   
print tracer state, general values (mostly for debugging)

last=0:  first    state
last=1:  last     state
last<0, the -last state

*/
void print_tracer_stats(struct trc *tracer,int last,FILE *out)
{
  COMP_PRECISION mstrain,val[3],xp[3];
  int astate;
  char label[10];
  if(last == 0){		/* initial */
    astate = 0;
    sprintf(label,"init:");
  }else if(last == 1){				/* final */
    astate = tracer->nstate - 1;
    sprintf(label,"fin: ");
  }else if(last < 0){
    astate = -last;
    if(astate >= tracer->nstate){
      fprintf(stderr,"print_tracer_stats: out of range: %i vs. %i states\n",
	      astate,tracer->nstate);
      exit(-1);
    }
    sprintf(label,"mid: ");
  }else{
    fprintf(stderr,"print_tracer_stats: last out of range: %i vs. %i states\n",
	    last,tracer->nstate);
    exit(-1);
  }
  /* 
     compute max strain 
  */
  mstrain = max_strain_from_def((tracer->state[astate].x+3),val);
  lonlatz_from_xp(tracer->state[astate].x,xp);
  fprintf(out,"tracer: %s mstrain: %11.3e time: %11g l/l/z: %11g %11g %11g EV(L): %11g %11g %11g PEV: %11g <ODF[%i]>: %g\n",
	  label,mstrain,tracer->state[astate].t,
	  xp[FSTRACK_X],xp[FSTRACK_Y],xp[FSTRACK_Z],sqrt(val[2]),sqrt(val[1]),
	  sqrt(val[0]),sqrt(val[2]*val[1]*val[0]),
	  tracer->state[astate].npole,
	  mean(tracer->state[astate].pdens,
	       tracer->state[astate].npole));
}

/* 
   print 

   lon lat depth time  

   to out WITHOUT NEWLINE
*/
void print_llz(struct stt *state, FILE *out)
{
  COMP_PRECISION xp[3];
  lonlatz_from_xp(state->x,xp);
  fprintf(out,"%13.5e %13.5e %13.5e ",xp[FSTRACK_X],xp[FSTRACK_Y],xp[FSTRACK_Z]);
}
/* 
   print lon lat depth time 

   to out WITHOUT NEWLINE
*/
void print_llzt(struct stt *state, FILE *out)
{
  print_llz(state,out);
  fprintf(out,"%8.3f ",state->t);
}
/* time lon lat depth NO NEWLINE */
void print_tllz(struct stt *state, FILE *out)
{
  fprintf(out,"%8.3f ",state->t);
  print_llz(state,out);

}
/* 


compute the best-fitting, fast TI axis and print to out



*/
void process_ti_tens_frac(COMP_PRECISION *x,
			  COMP_PRECISION *sav,
			  my_boolean rotate_coord,
			  FILE *out)
{
  COMP_PRECISION kmod,gmod,dc_iso[36],dc_hex[36],
    dc_tet[36],dc_ort[36],dc_mon[36],dc_tri[36],
    sav_scc[36],tiamp,symm_frac[6],vel[9],scc_irot[9];
  COMP_PRECISION evec[9],ecvec[9];
  int j,tioff;
  static my_boolean warned = FALSE;
  int scca_old_mode = 0;	/* use the fixed, new mode or, for
				   backward compatibility, use old
				   mode */
  /* 
     compute the best-fit hexagonal axis and percent anisotropy in
     the original, global Cartesian system
  */
  drex_decsym(sav,&kmod,&gmod,vel,symm_frac,
	      ecvec,dc_iso,dc_hex,dc_tet,dc_ort,dc_mon,
	      dc_tri,sav_scc,scc_irot,&scca_old_mode);
  /* 
     vs anisotropy in percent 
  */
  tiamp = 200.0*(vel[7]-vel[8])/(vel[7]+vel[8]);
  if(tiamp > 0)		/* fast symmetry */
    tioff = 0;
  else{			/* slow symmetry */
    tiamp = -tiamp;
    tioff = 3;
  }
  /* 
     convert vector to spherical system 
  */
  if(!rotate_coord){
    if(!warned){		/* exception */
      fprintf(stderr,"process_ti_tens_frac: WARNING: TI axes is not rotated with xp, global Cartesian\n");
      warned = TRUE;
    }
    a_equals_b_vector3d(evec,(ecvec+tioff));
  }else{
    /* 
       normal operation, rotate vector into spherical 

    */
    cart_vec2polar_vec_at_xp((ecvec+tioff),evec,x);
  }
  /* output, vector and amplitude */
  fprintf(out,"%13.6e %13.6e %13.6e\t %13.6e\t",
	  evec[FSTRACK_R],evec[FSTRACK_THETA],evec[FSTRACK_PHI],tiamp);
  /* symmetry components */
  for(j=0;j<6;j++)
    fprintf(out,"%7.5f ",symm_frac[j]);
}

/* 

rotate a matrix at xp[3] into the regular, Cartesian system with x, y, z being S, E, U

where xp[3] are the spherical coordinate r,theta,phi
*/
void rotate_mat_to_reg_cart_system(COMP_PRECISION *sav, COMP_PRECISION *xp,
				   COMP_PRECISION *savr)
{
  /* get the Euler angles to go from Cartesian format to x=South, y=East, z=up */
  static int coord_convention = DREX_REG_CART_CONVENTION;
  COMP_PRECISION alpha,beta,gamma;
  /* get the  Euler angles */
  drex_calc_euler_rad_from_xp(xp,&alpha,&beta,&gamma,&coord_convention);
  /* rotate the matrix */
  drex_rotate_6x6_rad_ftrn(sav,savr,&alpha,&beta,&gamma);
}
