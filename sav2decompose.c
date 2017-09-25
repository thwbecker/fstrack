/* 

convert Sav stiffness tensor format to decomposition parts


sav2decompose [input_mode, 3] [olivine_fraction, 0.7] ... 
    [VHR_avg_fac, 0] [density, 3.353] [delta_T, 0.0] [scca_old_mode, 1]

input format:

lon lat z upper_triangl_C_ij (j fast), i.e.
			  
           SAV_11 SAV_12 SAV_13 ... SAV_16 \
                  SAV_22 SAV_23 ... SAV_26 \
                            ... 
                                    SAV_66


output for regular operation modes:

ani hex tet ort mon tri c_33^0 c_44^0 epsilon*c_33^0 gamma*c_44^0 delta*c_33^0

if density is < 0, will use PREM

delta_T: temperature anomaly, to be added to regular temperature model

scca_old_mode:

scca_old_mode: 1: only find best-fit, and align
               0: find best and worst fit, and use those as a coordinate system

$Id: sav2decompose.c,v 1.16 2016/09/05 04:44:47 becker Exp $


*/
#include "fstrack.h"
#include "prem.h"
/* opmodes */
/*
			   
			   operational modes (read tensor and do something with it)
			   
			   3:  read matrix from file and decompose into symmetry components
			   
			   4:  read matrix and print ti axes

			   15: hexagonal factors a la Chevrot: x y z azi inc eps gamma delta

			       the eps, gamma, delta factors are in a best-fitting hexagonal (SCC) system!
			   
			   16: print best-fit hexagonal, in original system
			   24: best-fit hex in original system, short output format

			   17: print best-fit orthorhombic, SCC system

			   18: print tensor in SCC reference frame

			   22: compute some of the Love values for the full tensor, original system

			   23: full decomposition, fancy output


			   debug modes

			   0:  Jules' matrices 
			   5:  K & R  matrices
			   6:  f(p,T) old matrices at room p,T
			   10: f(p,T) new matrices at room p,T
			   9:  f(p,T) old matrices at room p = 5 GPa, T=1573K
			   13: f(p,T) new matrices at room p = 5 GPa, T=1573K
			   7:  T dependence of matrices
			   11: T dependence of new matrices
			   8:  p dependence of matrices
			   12: p dependence of new matrices
			   1:  depth dependent fractions on several levels
			   14: new depth dependent fractions on several levels
			   2:  depth dependent on Jules' levels

			   20: old Cij(z) 
			   21: new Cij(z)

			   25: old Cij(z from input) 
			   26: new Cij(z from input)

	
*/
/* operation */
#define SAV2DMODE_DECOMPOSE 3
#define SAV2DMODE_TI 4

#define SAV2DMODE_HEX_CHEVROT 15
#define SAV2DMODE_BESTFIT_HEX 16
#define SAV2DMODE_BESTFIT_HEX_SHORT 24
#define SAV2DMODE_BESTFIT_ORTH 17
#define SAV2DMODE_SCC 18
#define SAV2DMODE_COMPUTE_LOVE 22
#define SAV2DMODE_FULL_DECOMPOSE 23
 
/* debug */
#define SAV2DMODE_JULES 0
#define SAV2DMODE_KR 5
#define SAV2DMODE_OLEN_ROOM 6
#define SAV2DMODE_OLEN_NEW_ROOM 10
#define SAV2DMODE_OLEN_UM 9
#define SAV2DMODE_OLEN_NEW_UM 13
#define SAV2DMODE_DT_OLEN 7
#define SAV2DMODE_DT_OLEN_NEW 11
#define SAV2DMODE_DP_OLEN 8
#define SAV2DMODE_DP_OLEN_NEW 12
#define SAV2DMODE_DZ_OLEN 1
#define SAV2DMODE_DZ_OLEN_NEW 14
#define SAV2DMODE_DZJ_OLEN 2
/*  */
#define SAV2DMODE_OLD_CIJ_Z 20
#define SAV2DMODE_NEW_CIJ_Z 21

#define SAV2DMODE_OLD_CIJ_ZIN 25
#define SAV2DMODE_NEW_CIJ_ZIN 26

/* 

usage:  
 

          sav2decompose mode


*/
void print_tens_decomp(COMP_PRECISION *,my_boolean, int);

int main(int argc, char **argv)
{
  COMP_PRECISION sav[36],savi[36],savir[36],savt[36],x[3],remainder,perc,sav1[36],sav2[36],
    tiaxis[6],kmod,gmod,vs1,vs2,vsmean,vel[9],*zlevel=NULL,sav_scc[36],vsh,vsv,
    symm_frac[6],ddummy,ciso[36],t,p,z,dz,evec[3],tiamp,xp[3],scc_irot[9],
    chex[36],ctet[36],cort[36],cmon[36],ctri[36],xloc[3],vs[2],vp[2],
    inc,azi,cbfhex[36],cbfort[36],gl[2],cn[2],ca[2],ba[2],hf[2],lmod,nmod,
    amod,cmod,fmod;
  int i,j,n,nr,nz,tioff;
  FILE *in,*out;
  char types[2];
  /* 

  default values

  */
  /* olivine fraction */
  COMP_PRECISION ol_frac = 0.7;	/* for mix */
  /* averaging scheme */
  COMP_PRECISION vhrf = 0.0;			/* 0: voigt 0.5: VHR 1: Reuss  */
  my_boolean prem_dens = FALSE;
  COMP_PRECISION density = 3.353; /* for velocities */
  COMP_PRECISION delta_T = 0.0;	/* temperature offset  */
  my_boolean pfull = FALSE;
  int mode = SAV2DMODE_DECOMPOSE;		
  int scca_old_mode = 1;	/* default is old mode */
  /* 
     to rotate enstatite before additon 
  */
  /* 
     Euler angles in degrees for en rotation before 
     Voigt averaging for debugging modes  1 and 2
  */
  COMP_PRECISION euler[3] = ENSTATITE_EULER;
  struct prem_model prem[1];
  my_boolean verbose = FALSE;

  in = stdin;

  if(argc>1)
    sscanf(argv[1],"%i",&mode);	/* mode */
  if(argc>2)
    sscanf(argv[2],FLT_FORMAT,&ol_frac); /* olivine fraction */
  if(argc>3)			/* voigt reuss hill averaging */
    sscanf(argv[3],FLT_FORMAT,&vhrf);
  if(argc>4)			/* density */
    sscanf(argv[4],FLT_FORMAT,&density);
  if(argc>5)			/* temperature anomaly */
    sscanf(argv[5],FLT_FORMAT,&delta_T);
  if(argc>6)			
    sscanf(argv[6],"%i",&scca_old_mode);
  
  if(density < 0){
    prem_dens = TRUE;
    fprintf(stderr,"%s: WARNING: varying density with z depth and PREM\n",argv[0]);
    prem_read_model(PREM_MODEL_FILE,prem,FALSE);
  }
  if(scca_old_mode)
    fprintf(stderr,"%s: using old mode of SCCA alignment\n",argv[0]);
  else
    fprintf(stderr,"%s: using new mode of SCCA alignment\n",argv[0]);
  if(delta_T != 0.0)
    fprintf(stderr,"%s: WARNING: temperature offset of %g\n",
	    argv[0],delta_T);
  
  if((verbose) || (mode != SAV2DMODE_DECOMPOSE)||
     (ol_frac != 0.7)|| (vhrf != 0.) ){
    fprintf(stderr,"%s: running in mode %i, olivine fraction %g%%, VHR: %g\n",
	    argv[0],mode,ol_frac*100,vhrf);
  }else{
    if(vhrf != 0.)
      fprintf(stderr,"%s: WARNIGN: VHR: %g\n",argv[0],vhrf);
  }
  if(ol_frac > 1){
    fprintf(stderr,"%s: olivine fraction should be between 0 and 1 (%g)\n",
	    argv[0],ol_frac);
    exit(-1);
  }

  if(mode == SAV2DMODE_DZJ_OLEN){
    /* init Jules' zlevels */
    nz=8;    
    my_vecrealloc(&zlevel,nz,"sav2decompose");
    zlevel[0]=60;zlevel[1]=115;zlevel[2]=185;zlevel[3]=220;zlevel[4]=265;
    zlevel[5]=310;zlevel[6]=355;zlevel[7]=400;
  }
  /* 
     switch input modes
  */
  switch(mode){
    /* 

    operational modes 

    */
  case SAV2DMODE_FULL_DECOMPOSE:
  case SAV2DMODE_COMPUTE_LOVE:
  case SAV2DMODE_DECOMPOSE:	
  case SAV2DMODE_TI:
  case SAV2DMODE_HEX_CHEVROT:
  case SAV2DMODE_BESTFIT_HEX:
  case SAV2DMODE_BESTFIT_HEX_SHORT:
  case SAV2DMODE_BESTFIT_ORTH:
  case SAV2DMODE_SCC:		/* all of these need to read in the tensor */
    /*  */
    if(verbose)
      fprintf(stderr,"%s: reading tensor from %s\n",argv[0],"stdin");
    n=0;
    /* read in location and depth */
    while(fscanf(in,THREE_FLT_FORMAT,xloc,(xloc+1),(xloc+2))==3){
      if(!read_sym_6by6(sav,in)){
	fprintf(stderr,"%s: error\n",argv[0]);
	exit(-1);
      }
      if(prem_dens){
	/* get PREM density for non-water layers */
	prem_get_rho(&density,1.0-xloc[2]/6371.0,prem);
	density /= 1000.0;
	if(density < 2.6)
	  density = 2.6;
      }
      /* 

      decompose the tensor
      
      */
      drex_decsym(sav,&kmod,&gmod,vel,symm_frac,tiaxis,
		  ciso,chex,ctet,cort,cmon,ctri,sav_scc,scc_irot,
		  &scca_old_mode);
      /* 
	 get anisotropy defined as (v_s1 - v_s2)/<v_s> 
	 for fast symmetry axis, vs1 > vs2
      */
      vs1 = vel[7]/sqrt(density);
      vs2 = vel[8]/sqrt(density);
      vsmean = (vs1+vs2)/2.0;
      
      switch(mode){
      case SAV2DMODE_DECOMPOSE:
	/* 
	   general output mode
	*/
	if(verbose){
	  fprintf(stderr,"%s: hex: vs1: %11g vs2: %11g vs anisotropy: %g %%\n",
		  argv[0],vs1,vs2,(vs1-vs2)/vsmean*100);
	  fprintf(stderr,"%s: tensor fraction: iso hex tet ort mono\n",argv[0]);
	}
	//print_6by6_nice(sav,stdout);fprintf(stderr,"\n");
	for(i=0;i < 6;i++)	
	  /*   
	       TENSOR COMPONENT DECOMPOSITION IN SCC SYSTEM

	       1   2   3   4   5   6
	       ani hex tet ort mon tri 
	       
	  */
	  fprintf(stdout,"%8.5f ",symm_frac[i]);
	/* 
	   7      8        9             10           11
	   c_{11/33}^0 c_{44/77}^0 epsilon*c_33^0 gamma*c_44^0 delta*c_33^0 
	   ciso(4,4) = ciso[21]


	*/
	fprintf(stdout,"%11g %11g %11g %11g %11g ",
		ciso[0],ciso[21],vel[4],vel[5],vel[6]);
	/*  */
	fprintf(stdout,"\n");
	break;
      case SAV2DMODE_TI:
	/* 
	   TI axes 
	*/
	/* anisotropy in percent */
	tiamp = (vs1-vs2)/vsmean*100.0;
	/* 
	   select the fast axis, best-fitting might be slow
	*/
	/*
	  if(tiamp > 0)
	  tioff = 0;
	  else{
	  tiamp = -tiamp;
	  tioff = 3;
	  }
	*/
	tioff = 0;		/* always pick the best-axis, even if
				   slow */
	
	fprintf(stdout,"%13.6e %13.6e %13.6e\t",
		xloc[0],xloc[1],xloc[2]);
	/* 

	convert local x,y,z system to r,theta,phi for consistency with
	fstrack output of TI axes

	*/
	fprintf(stdout,"%13.6e %13.6e %13.6e\t %13.6e\t",
		tiaxis[tioff+2],tiaxis[tioff+0],tiaxis[tioff+1],
		tiamp);
	for(j=0;j<6;j++)fprintf(stdout,"%7.5f ",symm_frac[j]);
	fprintf(stdout,"\n");
	break;
      case SAV2DMODE_HEX_CHEVROT:
	/* 

	orientation of TI axes and hexagonal factors, 
	the latter in SCC system
	
	*/
	fprintf(stdout,"%g %g %g\t",xloc[0],xloc[1],xloc[2]);
	/* azimuth and incidence */
	/* 
	   azimuth clockwise from North 
	*/
	azi = atan2(tiaxis[1],-tiaxis[0]);
	if(azi<0)
	  azi += TWOPI;
	inc = asin(-tiaxis[2]);	/* positive is downward */
	fprintf(stdout,"%g %g\t",azi, inc);
	/* gamma epsilon delta */
	fprintf(stdout,"%g %g %g\n",
		vel[5]/ciso[21],vel[4]/ciso[0],vel[6]/ciso[0]);
	break;
      case SAV2DMODE_BESTFIT_HEX: 
	/* 
	   best-fitting hexagonal tensor in orig system
	*/
	compute_best_fit_hex(ciso,chex,scc_irot,1,cbfhex);
	fprintf(stderr,"orig:\n");
	print_6by6_nice(sav,stderr);fprintf(stderr,"\n");
	fprintf(stderr,"best-fit hex (orig)\n");
	print_6by6_nice(cbfhex,stdout);fprintf(stderr,"\n");
	break;
      case SAV2DMODE_BESTFIT_HEX_SHORT: /* short version of above */
	compute_best_fit_hex(ciso,chex,scc_irot,1,cbfhex);
	fprintf(stdout,"%g %g %g\t",xloc[0],xloc[1],xloc[2]);
	print_sym_6by6(cbfhex,stdout);fprintf(stdout,"\n");
	break;
      case SAV2DMODE_BESTFIT_ORTH: /* need to sum everything below
				      orthogonal, in SCC system */
	for(i=0;i<36;i++)
	  cbfort[i] = ciso[i]+chex[i]+ctet[i]+cort[i];
	fprintf(stderr,"orig (SCC):\n"); /* rotated */
	print_6by6_nice(sav_scc,stderr);fprintf(stderr,"\n");
	fprintf(stderr,"iso (SCC):\n"); /* rotated */
	print_6by6_nice(ciso,stderr);fprintf(stderr,"\n");
	fprintf(stderr,"best-fit orthorhombic (SCC):\n");
	print_6by6_nice(cbfort,stdout);fprintf(stderr,"\n");
	break;
      case SAV2DMODE_SCC: /* print tensor in SCC frame */
	for(i=0;i<36;i++)
	  cbfort[i] = ciso[i]+chex[i]+ctet[i]+cort[i]+ctri[i]+cmon[i];
	fprintf(stderr,"reassembled (SCC):\n"); /* rotated */
	print_6by6_nice(cbfort,stderr);fprintf(stderr,"\n");

	fprintf(stderr,"scc (SCC):\n"); /* rotated */
	print_6by6_nice(sav_scc,stdout);fprintf(stderr,"\n");
	break;
	
      case SAV2DMODE_COMPUTE_LOVE: /* compute the Love parameters for an equivalent 
				      VTI medium
				   */
	/* 

	for this, we should/can use the full tensor

	*/
	//fprintf(stderr,"WARNING: not rotated\n");
	drex_compute_swpar_ftrn(sav,gl,cn,ca,ba,hf,
				&lmod,&nmod,&amod,&cmod,&fmod);
	vsh = sqrt(nmod / density);
	vsv = sqrt(lmod / density);
	
	fprintf(stdout,"%11g %11g %11g\t%11g %11g\t%11g %11g %11g\t%11g %11g %11g\t%11g %11g %11g\n",
		xloc[0],xloc[1],xloc[2],
		kmod,gmod,	/* K and G moduli, GPa */
		vsh,vsv,50.0*log10(vsh/vsv),
		nmod/lmod,	/* chi = (vsh/vsv)^2 */
		cmod/amod,	/* phi = (vpv/vph)^2 */
		fmod/(amod-2.0*lmod), /* eta = F/(A-2L) */
		hypot(ba[0],ba[1]), /* sqrt(B_c^2+B_s^2)/A */
		hypot(hf[0],hf[1]), /* sqrt(H_c^2+H_s^2)/F */
		hypot(gl[0],gl[1])  /* sqrt(G_c^2+G_s^2)/L */
		);
	break;
      case SAV2DMODE_FULL_DECOMPOSE: /* nice deompositiin */
	print_tens_decomp(sav,TRUE,scca_old_mode);
	break;
      default:
	fprintf(stderr,"mode error\n");
	exit(-1);
	break;
      }
      n++;
    }
    if(verbose)
      fprintf(stderr,"%s: decomposed %i tensors \n",
	      argv[0],n);
    break;

    /* 
       debugging modes
    */
  case SAV2DMODE_JULES:
  case SAV2DMODE_KR:
  case SAV2DMODE_OLEN_ROOM:
  case SAV2DMODE_OLEN_NEW_ROOM:
  case SAV2DMODE_OLEN_UM:
  case SAV2DMODE_OLEN_NEW_UM:
  case SAV2DMODE_DT_OLEN:
  case SAV2DMODE_DT_OLEN_NEW:
  case SAV2DMODE_DP_OLEN:
  case SAV2DMODE_DP_OLEN_NEW:
    /* 
       test decomposition using standard values 
    */
    /* 
       
    olivine
    
    */
    switch(mode){
    case SAV2DMODE_JULES:
      drex_get_sav_constants(sav,DREX_C_OL_JULES,0,0); /* Jules'  */
      break;
    case SAV2DMODE_KR:
      drex_get_sav_constants(sav,DREX_C_OL_KR,0,0); /* KR'  */
      break;
    case SAV2DMODE_OLEN_ROOM:
      drex_get_sav_constants(sav,DREX_CTP_OL_ESTEY,298.0,0); /* T, p (room) */
      break;
    case SAV2DMODE_OLEN_NEW_ROOM:
      drex_get_sav_constants(sav,DREX_CTP_OL_NEW,298.0,0); /* new T, p (room) */
      break;
    case SAV2DMODE_OLEN_UM:
      drex_get_sav_constants(sav, DREX_CTP_OL_ESTEY,1573.0,5.0); /* T, p (200km) */
      break;
    case SAV2DMODE_OLEN_NEW_UM:
      drex_get_sav_constants(sav,DREX_CTP_OL_NEW,1573.0,5.0); /* new T, p (200km) */
      break;
    case SAV2DMODE_DT_OLEN:
      drex_get_sav_constants(sav1,DREX_CTP_OL_ESTEY,298.0,0); /* dC/dT */
      drex_get_sav_constants(sav,DREX_CTP_OL_ESTEY,299.0,0); 
      sub_a_from_b_vector(sav1,sav,36);
      break;
    case SAV2DMODE_DT_OLEN_NEW:
      drex_get_sav_constants(sav1,DREX_CTP_OL_NEW,298.0,0); /* new dC/dT */
      drex_get_sav_constants(sav,DREX_CTP_OL_NEW,299.0,0); 
      sub_a_from_b_vector(sav1,sav,36);
      break;
    case SAV2DMODE_DP_OLEN:
      drex_get_sav_constants(sav1,DREX_CTP_OL_ESTEY,298.0,0); /* dC/dp */
      drex_get_sav_constants(sav,DREX_CTP_OL_ESTEY,298.0,1); 
      sub_a_from_b_vector(sav1,sav,36);
      break;
    case  SAV2DMODE_DP_OLEN_NEW:
      drex_get_sav_constants(sav1,DREX_CTP_OL_NEW,298.0,0); /* new dC/dp */
      drex_get_sav_constants(sav,DREX_CTP_OL_NEW,298.0,1); 
      sub_a_from_b_vector(sav1,sav,36);
      break;
    default:
      fprintf(stderr,"mode error %i\n",mode);
      exit(-1);
    }
    print_tens_decomp(sav,pfull,scca_old_mode);
    /* 
       enstatite at p, T
    */
    switch(mode){
    case SAV2DMODE_JULES:
      drex_get_sav_constants(savi,DREX_C_EN_JULES,0,0); /* Jules'  */
      break;
    case SAV2DMODE_KR:
      drex_get_sav_constants(savi,DREX_C_EN_KR,0,0); /* KR'  */
      break;
    case SAV2DMODE_OLEN_ROOM:
      drex_get_sav_constants(savi,DREX_CTP_EN_ESTEY,298.0,0); /* T, p @ room */
      break;
    case SAV2DMODE_OLEN_NEW_ROOM:
      drex_get_sav_constants(savi,DREX_CTP_EN_NEW,298.0,0); /* new T, p @ room */
      break;
    case SAV2DMODE_OLEN_UM:
      drex_get_sav_constants(savi,DREX_CTP_EN_ESTEY,1573.0,5.0); /* T, p @ 200 km*/
      break;
    case SAV2DMODE_OLEN_NEW_UM:
      drex_get_sav_constants(savi,DREX_CTP_EN_NEW,1573.0,5.0); /* new T, p @ 200 km*/
      break;
    case SAV2DMODE_DT_OLEN:
      drex_get_sav_constants(sav1,DREX_CTP_EN_ESTEY,298.0,0); /* dC/dT */
      drex_get_sav_constants(savi,DREX_CTP_EN_ESTEY,299.0,0); 
      sub_a_from_b_vector(sav1,savi,36);
      break;
    case SAV2DMODE_DT_OLEN_NEW:
      drex_get_sav_constants(sav1,DREX_CTP_EN_NEW,298.0,0); /* new dC/dT */
      drex_get_sav_constants(savi,DREX_CTP_EN_NEW,299.0,0); 
      sub_a_from_b_vector(sav1,savi,36);
      break;
    case SAV2DMODE_DP_OLEN:
      drex_get_sav_constants(sav1,DREX_CTP_EN_ESTEY,298.0,0); /* dC/dp */
      drex_get_sav_constants(savi,DREX_CTP_EN_ESTEY,298.0,1); 
      sub_a_from_b_vector(sav1,savi,36);
      break;
    case SAV2DMODE_DP_OLEN_NEW:
      drex_get_sav_constants(sav1, DREX_CTP_EN_NEW,298.0,0); /* new dC/dp */
      drex_get_sav_constants(savi, DREX_CTP_EN_NEW,298.0,1); 
      sub_a_from_b_vector(sav1,savi,36);
      break;
    default:
      fprintf(stderr,"mode error %i\n",mode);
      exit(-1);
    }
    print_tens_decomp(savi,pfull,scca_old_mode);
    /* 
       assemblage
    */
    drex_rotate_6x6_deg_ftrn(savi,savir,(euler),(euler+1),(euler+2)); /* rotate ens like below */
    zero_small_entries(savir,36);
    drex_st_vhr_avg_ftrn(sav,savir,&ol_frac,&vhrf,savt); /* get VHR average */
    print_tens_decomp(savt,pfull,scca_old_mode);
    /* file out */
    out = myopen("tmp.sav.out","w",argv[0]);fprintf(out,"0 0 0 ");
    print_sym_6by6(savt,out);
    fclose(out);
    break;
  case SAV2DMODE_DZ_OLEN:			/* closely spaced tables */
  case SAV2DMODE_DZ_OLEN_NEW:
    nz = 103;dz=410/(nz-1);
    my_vecrealloc(&zlevel,nz,"sav2decompose");
    for(i=0;i<nz;i++)
      zlevel[i] = 0 + dz*i;
  case SAV2DMODE_DZJ_OLEN:
    if(!prem_dens)		/* init PREM */
      prem_read_model(PREM_MODEL_FILE,prem,FALSE);
    /*  */
    fprintf(stderr,"%s: output: z T p   K G   ani orth hex mono tet tri\teps gamma delta\t vp0 vs0 vp1 vs1\tchi phi eta\n",argv[0]);
    for(i=0;i < nz;i++){
      z = zlevel[i];
      t = temperature_model(z,1) + delta_T;
      p = pressure_model(z);

      if((mode == SAV2DMODE_DZ_OLEN)||
	 (mode == SAV2DMODE_DZJ_OLEN)){		/* old */
	/* 
	   get single crystal constants, old derivatives
	*/
	drex_get_sav_constants(sav,DREX_CTP_OL_ESTEY,t,p); /* olivine */
	drex_get_sav_constants(savi,DREX_CTP_EN_ESTEY,t,p); /* enstatite */
      }else{			
	/* 
	   new derivatives
	*/
	drex_get_sav_constants(sav,DREX_CTP_OL_NEW,t,p); /* olivine */
	drex_get_sav_constants(savi,DREX_CTP_EN_NEW,t,p); /* enstatite */
      }
      /* 

      rotate enstatite such that en [001] with ol [100]
                                    [100] with ol [100]
      
      */
      drex_rotate_6x6_deg_ftrn(savi,savir,(euler),(euler+1),(euler+2));
      zero_small_entries(savir,36);
      drex_st_vhr_avg_ftrn(sav,savir,&ol_frac,&vhrf,savt); /* get VHR average */
      /* 
	 decompose
      */
      drex_decsym(savt,&kmod,&gmod,vel,symm_frac,tiaxis,ciso,
		  chex,ctet,cort,cmon,ctri,sav_scc,scc_irot,
		  &scca_old_mode);

      /* compute equivalent Love parameters based on the original,
	 full tensor */
      drex_compute_swpar_ftrn(sav,gl,cn,ca,ba,hf,&lmod,
			      &nmod,&amod,&cmod,&fmod);
      /* output of 

      z T p   K G   ani orth hex mono tet tri   eps gamma delta (all SCC) vp[0] vs[0] vp[1] vs[1]   xi phi eta (orig frame)

      */
      /* get density */
      prem_get_rho(&density,1.0-z/6371.0,prem);
      /* get velocities */
      /* 
	 based on isotropic projection
      */
      vp[0] = sqrt(ciso[0]*1e9/density); /* C_11 */
      vs[0] = sqrt(ciso[21]*1e9/density); /* C_44 */
      /* 
	 based on means from hexagonal approximation
      */
      vp[1] = sqrt(vel[2]*1e9/density);
      vs[1] = sqrt(vel[3]*1e9/density);
      fprintf(stdout,"%11g %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f  %9.3f  %9.3f\t%9.5f  %9.5f  %9.5f\t%g %g\t%g %g\t%g %g %g\n",
	      z,t,p,kmod,gmod,symm_frac[0]*100,
	      symm_frac[3]*100,symm_frac[1]*100,symm_frac[2]*100,symm_frac[4]*100,
	      symm_frac[5]*100,
	      vel[4]/ciso[0],	/* eps SCC */
	      vel[5]/ciso[21],	/* gamma SCC */
	      vel[6]/ciso[0],	/* delta SCC */
	      vp[0],vs[0],vp[1],vs[1], /* velocities */
	      nmod/lmod,
	      cmod/amod,
	      fmod/(amod-2.0*lmod) /* chi, phi, eta (original frame) */
	      );

    }
    break;
  case SAV2DMODE_OLD_CIJ_Z:
  case SAV2DMODE_NEW_CIJ_Z:
    nz = 103;dz=410/(nz-1);
    my_vecrealloc(&zlevel,nz,"sav2decompose");
    for(i=0;i < nz;i++)
      zlevel[i] = 0 + dz*i;
    fprintf(stderr,"%s: output: z T p Cij\n",argv[0]);
    for(i=0;i < nz;i++){
      z = zlevel[i];
      t = temperature_model(z,1) + delta_T;
      p = pressure_model(z);
      if(mode == SAV2DMODE_OLD_CIJ_Z){		/* old */
	/* 
	   get single crystal constants, old derivatives
	*/
	drex_get_sav_constants(sav,DREX_CTP_OL_ESTEY,t,p); /* olivine */
	drex_get_sav_constants(savi,DREX_CTP_EN_ESTEY,t,p); /* enstatite */
      }else{			
	/* 
	   new derivatives
	*/
	drex_get_sav_constants(sav,DREX_CTP_OL_NEW,t,p); /* olivine */
	drex_get_sav_constants(savi,DREX_CTP_EN_NEW,t,p); /* enstatite */
      }
      /* 

      rotate enstatite such that en [001] with ol [100]
                                    [100] with ol [100]
      
      */
      drex_rotate_6x6_deg_ftrn(savi,savir,(euler),(euler+1),(euler+2));
      zero_small_entries(savir,36);
      drex_st_vhr_avg_ftrn(sav,savir,&ol_frac,&vhrf,savt); /* get VHR average */
      /* 
	 output of 
	 
	 T p z Cij

      */
      fprintf(stdout,"%11g %9.3f %9.3f ",t,p,z);
      print_sym_6by6(savt,stdout);
	
    }
    break; 
  case SAV2DMODE_OLD_CIJ_ZIN:
  case SAV2DMODE_NEW_CIJ_ZIN:
    fprintf(stderr,"%s: which depth to evaluate tensor? [km]\n",argv[0]);
    fscanf(stdin,FLT_FORMAT,&z);
    fprintf(stderr,"%s: output: z T p Cij\n",argv[0]);
    t = temperature_model(z,1) + delta_T;
    p = pressure_model(z);
    if(mode == SAV2DMODE_OLD_CIJ_Z){		/* old */
      /* 
	 get single crystal constants, old derivatives
      */
      drex_get_sav_constants(sav,DREX_CTP_OL_ESTEY,t,p); /* olivine */
      drex_get_sav_constants(savi,DREX_CTP_EN_ESTEY,t,p); /* enstatite */
    }else{			
      /* 
	 new derivatives
      */
      drex_get_sav_constants(sav,DREX_CTP_OL_NEW,t,p); /* olivine */
      drex_get_sav_constants(savi,DREX_CTP_EN_NEW,t,p); /* enstatite */
    }
    /* 
       
       rotate enstatite such that en [001] with ol [100]
       [100] with ol [100]
       
    */
    drex_rotate_6x6_deg_ftrn(savi,savir,(euler),(euler+1),(euler+2));
    zero_small_entries(savir,36);
    drex_st_vhr_avg_ftrn(sav,savir,&ol_frac,&vhrf,savt); /* get VHR average */
    /* 
       output of 
       
       T p z Cij
       
    */
    fprintf(stdout,"%11g %9.3f %9.3f ",t,p,z);
    print_sym_6by6(savt,stdout);
    break;
  default:
    fprintf(stderr,"%s: mode %i undefined\n",argv[0],mode);
    exit(-1);
    break;
  }
  
}
  
void print_tens_decomp(COMP_PRECISION *sav, my_boolean pfull, int scca_old_mode)
{
  COMP_PRECISION kmod,gmod,symm_frac[6],tiaxis[6],vel[9],scc_irot[9],
    ciso[36],chex[36],ctet[36],cort[36],cmon[36],ctri[36],sav_scc[36];
  
  drex_decsym(sav,&kmod,&gmod,vel,symm_frac,tiaxis,ciso,chex,ctet,cort,cmon,ctri,
	      sav_scc,scc_irot,&scca_old_mode);
  fprintf(stderr,"%g %g\n",kmod,gmod);
  fprintf(stderr,"ani hex tet ort mon tri\n");
  print_vector(symm_frac,6,stderr);
  fprintf(stderr,"ori:\n");
  print_6by6_nice(sav,stdout);fprintf(stderr,"\n");
  if(pfull){
    fprintf(stderr,"ori_scc:\n"); /* original rotated into SCC system */
    print_6by6_nice(sav_scc,stdout);fprintf(stderr,"\n");
    fprintf(stderr,"iso:\n");	/* isotropic projection of original in SCC */
    print_6by6_nice(ciso,stdout);fprintf(stderr,"\n");
    fprintf(stderr,"hex:\n");	/* hexgonal projection  */
    print_6by6_nice(chex,stdout);fprintf(stderr,"\n");
    fprintf(stderr,"tet:\n");	/* tetragonal projection */
    print_6by6_nice(ctet,stdout);fprintf(stderr,"\n");
    fprintf(stderr,"ort:\n");	/* orthorhombic projection */
    print_6by6_nice(cort,stdout);fprintf(stderr,"\n");
    fprintf(stderr,"tri:\n");	/* triclinic projection */
    print_6by6_nice(ctri,stdout);fprintf(stderr,"\n");
  }
}
/* 

given an isotropic (ciso) and hex component (chex) tensor in the SCC
system (rotation matrix: rot) as produced by drex_decsym, return the
best fitting hexagonal tensor (cbfhex) in the original coordinate
frame (orig_sys = 1), or in the rotated SCC frame

 */
void compute_best_fit_hex(COMP_PRECISION *ciso,COMP_PRECISION *chex, 
			  COMP_PRECISION *irot, int orig_sys,
			  COMP_PRECISION *cbfhex)
{
  int i;
  COMP_PRECISION hr[36];
  for(i=0;i<36;i++)
    cbfhex[i] = ciso[i]+chex[i];
  if(orig_sys){
    /* rotate back */
    drex_rotate_6x6_rot_ftrn(cbfhex,irot,hr);
    a_equals_b_vector(cbfhex,hr,36);
    zero_small_entries(cbfhex,36);
  }

}
