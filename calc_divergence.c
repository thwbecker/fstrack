#include "fstrack_flow.h"

/*

  debugging routine that calculates the divergence
  using central differences given a discrete flow
  field as read in by read_vel_grids
  
  $Id: calc_divergence.c,v 1.6 2016/09/05 04:44:47 becker Exp $

*/


void calc_divergence(struct mod *model)
{
  COMP_PRECISION ddx[3],// d_r v_r, d_phi v_phi, d_theta v_theta
    *cottheta,*invsintheta,theta,invr,dx[3],div;
  int i,j,k,index,id[3][2],nm[3];
  // 1/sin(theta) and 1/tan(theta) factors
  cottheta=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*MDP->n[FSTRACK_THETA]);
  invsintheta=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*MDP->n[FSTRACK_THETA]);
  if(!cottheta || !invsintheta)MEMERROR("calc_divergence");
  for(theta=MDP->dtheta*0.5,i=0;i<MDP->n[FSTRACK_THETA];i++,theta+=MDP->dtheta){
    invsintheta[i] = 1.0/sin(theta);
    cottheta[i]=1.0/tan(theta);
  }
  for(i=0;i<3;i++)// last entries
    nm[i]=MDP->n[i]-1;
  dx[FSTRACK_PHI]= 2.0*MDP->dphi;
  //for(i=0;i<MDP->n[FSTRACK_R];i++){
  for(i=MDP->n[FSTRACK_R]/2;i<=MDP->n[FSTRACK_R]/2;i++){
    invr = 1.0/MDP->rlevels[i];
    if(i == 0){// forward difference in r
      id[FSTRACK_R][1]=1;
      id[FSTRACK_R][0]=0;
    }else if(i == nm[FSTRACK_R]){// backward in r
      id[FSTRACK_R][1]= i;
      id[FSTRACK_R][0]= i - 1;
    }else{// central
      id[FSTRACK_R][1]=i + 1;
      id[FSTRACK_R][0]=i - 1;
    }
    dx[FSTRACK_R] = MDP->rlevels[id[FSTRACK_R][1]] - MDP->rlevels[id[FSTRACK_R][0]];
    for(j=0;j<MDP->n[FSTRACK_THETA];j++){
      if(j == 0){// forward difference in theta
	id[FSTRACK_THETA][1]=1;
	id[FSTRACK_THETA][0]=0;
	dx[FSTRACK_THETA] = MDP->dtheta;
      }else if(j == nm[FSTRACK_THETA]){// backward in theta
	id[FSTRACK_THETA][1]=j;
	id[FSTRACK_THETA][0]=j-1;
	dx[FSTRACK_THETA] = MDP->dtheta;
      }else{// central 
	id[FSTRACK_THETA][1]=j+1;
	id[FSTRACK_THETA][0]=j-1;
	dx[FSTRACK_THETA] = 2.0 * MDP->dtheta;
      }
      for(k=0;k<MDP->n[FSTRACK_PHI];k++){
	if(k==0){
	  id[FSTRACK_PHI][1] = 1;
	  id[FSTRACK_PHI][0] = nm[FSTRACK_PHI];
	}else if(k == nm[FSTRACK_PHI]){
	  id[FSTRACK_PHI][1] = 0;
	  id[FSTRACK_PHI][0] = nm[FSTRACK_PHI]-1;
	}else{
	  id[FSTRACK_PHI][1] = k + 1;
	  id[FSTRACK_PHI][0] = k - 1;
	}
	// d_r v_r
	ddx[FSTRACK_R] = ((COMP_PRECISION)MDP->vr[VOFF(id[FSTRACK_R][1],j,k)] - (COMP_PRECISION)MDP->vr[VOFF(id[FSTRACK_R][0],j,k)])/
	  dx[FSTRACK_R];
	// d_t v_theta
	ddx[FSTRACK_THETA] = ((COMP_PRECISION)MDP->vt[VOFF(i,id[FSTRACK_THETA][1],k)] - (COMP_PRECISION)MDP->vt[VOFF(i,id[FSTRACK_THETA][0],k)])/
	  dx[FSTRACK_THETA];
	// d_p v_phi
	ddx[FSTRACK_PHI] = ((COMP_PRECISION)MDP->vp[VOFF(i,j,id[FSTRACK_PHI][1])] - (COMP_PRECISION)MDP->vp[VOFF(i,j,id[FSTRACK_PHI][0])])/
	  dx[FSTRACK_PHI];
	
	div  = ddx[FSTRACK_R];
	index = VOFF(i,j,k);
	div += invr*(2.0*(COMP_PRECISION)MDP->vr[index] + ddx[FSTRACK_THETA] + 
		     (COMP_PRECISION)MDP->vt[index] * cottheta[j] +
		     ddx[FSTRACK_PHI] * invsintheta[j] );
	fprintf(stderr,"%3i %3i %3i div: %11g\n",i,j,k,div);
      }
    }
  }

  free(cottheta);free(invsintheta);
}

