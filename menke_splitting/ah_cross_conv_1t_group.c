#include <stdio.h>
#include <math.h>
#include <rpc/rpc.h>
#include "SHahhead.h"

#define DTOR (3.141592654/180.0)
#define TRUE (1)
#define FALSE (0)
#define MAXTRACES (100)

double dot();
double myrand();
char *progname;

main(int argc, char *argv[])
{
	FILE *f1, *f2, *f3, *f4, *f5, *f6;
	XDR in1, in2, out;
	ahhed head1[MAXTRACES], head2[MAXTRACES];
	float *data1[MAXTRACES], *data2[MAXTRACES], *data3, *data4;
	int i, j, n, ii, jj, kk, first, ntraces, nullfile;
	int tau1, maxn;
	double theta, phi, ct, st, cp, sp, a[4], b[4];
	double  E, Emin, taumax, thetamin, phimin;
	double Ea, Eb, R, delta, Rall, Rallmin;
	double bazi[MAXTRACES], vh[MAXTRACES];
	int clag[4], rtflag;
	double min_h, max_h, d_h, min_B, max_B, d_B, min_a1, max_a1, d_a1, min_a2, max_a2, d_a2;
	double c[3][3][3][3], c2[3][3][3][3], rot[3][3], Vs;
	char north[256], east[256];
	double ncoeff[2], ecoeff[2], zcoeff[2], radialcoeff[2], tangentialcoeff[2], time[2];
	double hmin, Bmin, a1min, a2min, vbackground, a1, a2, a3, B, h, rho;
	double v1, v2, v3;
	double ta[2], tb[2];
	int flag;

	progname=argv[0];

	if( argc != 18 ) {
		fprintf(stderr,"usage: %s min_h max_h d_h min_B max_B d_B min_a1 max_a1 d_a1 min_a2 max_a2 d_a2 filelist.txt output.ah error_surface.txt tensor.txt flag > results.txt\n",argv[0]);
		fprintf(stderr,"    min_h, minimum layer thickness, in km\n");
		fprintf(stderr,"    max_h, maximum layer thickness, in km\n");
		fprintf(stderr,"    d_h, thickness increment, in km\n");
		fprintf(stderr,"    min_B, maximum oblateness, >= -1.0\n");
		fprintf(stderr,"    max_B, maximum oblateness, <=  1.0\n");
		fprintf(stderr,"    d_B, oblateness increment\n");
		fprintf(stderr,"    min_a1, minimum first euler angle, in deg\n");
		fprintf(stderr,"    max_a1, maximum first euler angle, in deg\n");
		fprintf(stderr,"       note: maximum range [0, 180] deg]\n");
		fprintf(stderr,"       note: a1 is azimuth of tensor symmetry axis in deg E of N\n");
		fprintf(stderr,"    d_a1, first euler angle increment, in deg\n");
		fprintf(stderr,"    min_a2, minimum second euler angle (dip of axis), in deg\n");
		fprintf(stderr,"    max_a2, maximum second euler angle, in deg\n");
		fprintf(stderr,"       note: maximum range [0, 180] deg]\n");
		fprintf(stderr,"       note: a2 is dip of tensor symmetry axis in deg\n");
		fprintf(stderr,"    d_a2, second euler angle increment, in deg\n");
		fprintf(stderr,"    filelist.txt, a file containing one line for each observation\n");
		fprintf(stderr,"    	with four blank-separated fioelds on each line:\n");
		fprintf(stderr,"    	field 1: backazimuth, in degrees\n");
		fprintf(stderr,"    	field 2: horizontal phase velocity, in km/s\n");
		fprintf(stderr,"    	field 3: path of AH format north component file\n");
		fprintf(stderr,"    	field 4: path of AH format east  component file\n");
		fprintf(stderr,"flag, y to include a free-surface correction, n otherwise.\n");

		fprintf(stderr,"\n");
		fprintf(stderr,"purpose: estimate simple one-layer anisotropy parameters\n");
		fprintf(stderr,"    from a group of radially-polarized (e.g. SKS) shear phases\n");
		fprintf(stderr,"    using the Menke & Levin (2002) cross-convolution Method\n");
		fprintf(stderr,"\n");
		fprintf(stderr,"Synopsis of method:\n");
		fprintf(stderr,"    The observed north A(t)=s(t)*a(t) and east B(t)=s(t)*b(t) data\n");
		fprintf(stderr,"    are assumed related to an unknown source wavelet, s(t) and\n");
		fprintf(stderr,"    impulse respose fuctions, a(t) and b(t).  The response functions\n");
		fprintf(stderr,"    are calculated for a simple anisotropic layer of thickness, h,\n");
		fprintf(stderr,"    and standard anisotropic parameters Vs=4.5 km/s, Vp=1.8*Vs,\n");
		fprintf(stderr,"    oblateness parameter B, C=0.0, E=B and rho=3300 kg/m3, and euler\n");
		fprintf(stderr,"    angles, a1, a2, a3=0.0. The oblateness parameter of the hexagonally\n");
		fprintf(stderr,"    symmetric tensor (symmetry axis is horizontal and north/south\n");
		fprintf(stderr,"    before the rotation) can vary between [-1, 1], although \n");
		fprintf(stderr,"    values between [-0.1, 0.1] are typical of the earth\n");
		fprintf(stderr,"\n");
		fprintf(stderr,"    This program determines h, B, a1 and a2 by a grid search by\n");
		fprintf(stderr,"    minimizing the mean (over the group of N timeseries) of the\n");
		fprintf(stderr,"    amplitude-normalized squared L2 norm, R, of the cross-convolution,\n");
		fprintf(stderr,"            e(t) = A(t)*b(t) - B(t)*a(t).\n");
		fprintf(stderr,"    Here R = (1/N) (sum) ||e||**2 / [ ||A(t)*b(t)||**2 + ||B(t)*a(t)||**2 ]\n");
		fprintf(stderr,"\n");
		fprintf(stderr,"output.ah, ah file with best-fitting\n");
		fprintf(stderr,"    A(t)*b(t)\n");
		fprintf(stderr,"    B(t)*a(t)\n");
		fprintf(stderr,"\n");
		fprintf(stderr,"error_surface.txt, a line with a zero on it followed by a table of:\n");
		fprintf(stderr,"    h, layer thickness in km\n");
		fprintf(stderr,"    B, oblateness (dimensionless)\n");
		fprintf(stderr,"    a1, euler angle 1, in deg\n");
		fprintf(stderr,"    a2, euler angle 2, in deg\n");
		fprintf(stderr,"    R, value of amplitude-normalized squared L2 norm\n");
		fprintf(stderr,"    Note: no output written if filenmame is /dev/null\n");
		fprintf(stderr,"          output written to stderr filename is -\n");
		fprintf(stderr,"\n");
		fprintf(stderr,"results.txt, a table of:\n");
		fprintf(stderr,"    h, best-fitting layer thickness in km\n");
		fprintf(stderr,"    B, best-fitting oblateness parameter (dimensionless)\n");
		fprintf(stderr,"    a1, a2, best-fitting euler angles in degrees\n");
		fprintf(stderr,"    E, minimum value of amplitude-normalized squared L2 norm\n");
		fprintf(stderr,"\n");
		fprintf(stderr,"tensor.txt, resulting tensor, in standard format\n");

		fprintf(stderr,"\n");

		euler_blerb();
		exit(-1);
		}

	if( sscanf(argv[1],"%le",&min_h) != 1 ) {
		fprintf(stderr,"error: %s: can read min_h from command line\n", argv[0] );
		exit(-1);
		}

	if( sscanf(argv[2],"%le",&max_h) != 1 ) {
		fprintf(stderr,"error: %s: can read max_h from command line\n", argv[0] );
		exit(-1);
		}

	if( sscanf(argv[3],"%le",&d_h) != 1 ) {
		fprintf(stderr,"error: %s: can read d_h from command line\n", argv[0] );
		exit(-1);
		}

	if( sscanf(argv[4],"%le",&min_B) != 1 ) {
		fprintf(stderr,"error: %s: can read min_B from command line\n", argv[0] );
		exit(-1);
		}

	if( sscanf(argv[5],"%le",&max_B) != 1 ) {
		fprintf(stderr,"error: %s: can read max_B from command line\n", argv[0] );
		exit(-1);
		}

	if( sscanf(argv[6],"%le",&d_B) != 1 ) {
		fprintf(stderr,"error: %s: can read d_B from command line\n", argv[0] );
		exit(-1);
		}

	if( sscanf(argv[7],"%le",&min_a1) != 1 ) {
		fprintf(stderr,"error: %s: can read min_first_euler_angle from command line\n", argv[0] );
		exit(-1);
		}

	if( sscanf(argv[8],"%le",&max_a1) != 1 ) {
		fprintf(stderr,"error: %s: can read max_first_euler_angle from command line\n", argv[0] );
		exit(-1);
		}

	if( sscanf(argv[9],"%le",&d_a1) != 1 ) {
		fprintf(stderr,"error: %s: can read first_euler_angle_increment from command line\n", argv[0] );
		exit(-1);
		}

	if( sscanf(argv[10],"%le",&min_a2) != 1 ) {
		fprintf(stderr,"error: %s: can read min_second_euler_angle from command line\n", argv[0] );
		exit(-1);
		}

	if( sscanf(argv[11],"%le",&max_a2) != 1 ) {
		fprintf(stderr,"error: %s: can read max_second_euler_angle from command line\n", argv[0] );
		exit(-1);
		}

	if( sscanf(argv[12],"%le",&d_a2) != 1 ) {
		fprintf(stderr,"error: %s: can read first_second_angle_increment from command line\n", argv[0] );
		exit(-1);
		}

	if( (f1=fopen(argv[13],"r"))==NULL ) {
		fprintf(stderr,"error: %s: cant open filelist.txt <%s> for read\n", argv[0], argv[13]);
		exit(-1);
		}

	if( (f2=fopen(argv[14],"w"))==NULL ) {
		fprintf(stderr,"error: %s: cant open output_ah_file <%s> for write\n", argv[0], argv[14]);
		exit(-1);
		}
	xdrstdio_create(&out,  f2,  XDR_ENCODE);

	nullfile = !strcmp(argv[15],"/dev/null");
	if( !nullfile ) {
		if( !strcmp(argv[15],"-") ) f3=stderr;
		else if( (f3=fopen(argv[15],"w"))==NULL ) {
			fprintf(stderr,"error: %s: cant open error_surface_text_file <%s> for write\n", argv[0], argv[15]);
			exit(-1);
			}
		}

	if( (f6=fopen(argv[16],"w"))==NULL ) {
		fprintf(stderr,"error: %s: cant open tensor.txt <%s> for write\n", argv[0], argv[16]);
		exit(-1);
		}

	if( argv[17][0]=='y' || argv[17][0]=='Y' ) flag=1;
	else if( argv[17][0]=='n' || argv[17][0]=='N' ) flag=0;
	else {
		fprintf(stderr,"error: %s: cbad free surface flag <%s> (must be y or n)\n", argv[0], argv[17]);
		exit(-1);
		}

	j=0;
	while(  fscanf(f1,"%le%le%s%s",&(bazi[j]),&(vh[j]),north,east) == 4 ) {

		if( (f4=fopen(north,"r"))==NULL ) {
			fprintf(stderr,"error: %s: cant open north <%s> for read\n", argv[0], north);
			exit(-1);
			}
		if( (f5=fopen(east,"r"))==NULL ) {
			fprintf(stderr,"error: %s: cant open east <%s> for read\n", argv[0], east);
			exit(-1);
			}

		xdrstdio_create(&in1,  f4,  XDR_DECODE);
		xdrstdio_create(&in2,  f5,  XDR_DECODE);

		if( xdr_getrecord2(&(head1[j]),&(data1[j]), &in1) != 0 ) {
			fprintf(stderr,"error: %s: cant read north ah_file\n", argv[0], north );
			exit(-1);
			}

		if( j==0 ) {
			delta=head1[j].record.delta;
			maxn=head1[j].record.ndata;
			}
		else if( head1[j].record.ndata > maxn ) {
			maxn=head1[j].record.ndata;
			}

		if( xdr_getrecord2(&(head2[j]),&(data2[j]), &in2) != 0 ) {
			fprintf(stderr,"error: %s: cant read ah_file_2\n", argv[0], east );
			exit(-1);
			}

		xdr_destroy(&in1); xdr_destroy(&in2); fclose(f4); fclose(f5);

		if( (head1[j].record.delta != delta) || (head2[j].record.delta != delta) ) {
			fprintf(stderr,"error: %s: deltas dont match (timeseries %d)\n", argv[0], j);
			exit(-1);
			}

		if( head1[j].record.ndata != head2[j].record.ndata ) {
			fprintf(stderr,"error: %s: ndatas dont match (timeseries %d)\n", argv[0], j);
			exit(-1);
			}
		j++;
		}
		
	ntraces=j;

	if( (data3=(float*)malloc(maxn*sizeof(float))) == NULL ) {
		fprintf(stderr,"error: %s: out of memory!\n", argv[0]);
		exit(-1);
		}
	if( (data4=(float*)malloc(maxn*sizeof(float))) == NULL ) {
		fprintf(stderr,"error: %s: out of memory!\n", argv[0]);
		exit(-1);
		}

	first=TRUE;

	for( B=min_B; B<=max_B; B+=d_B ) {

		/* setup tensor Vp=1.8*Vs, Vs=4500, B, C=0, E=B, rho=3300 */
		/* scale factor of 1e-11 gets you units of Mbar */
		Vs=4500.0;
		rho=3300.0;
		one_axis_tensor(1.8*Vs, Vs, B, 0.00, B, rho*1.0e-11, c2 );

		/* one_axis_tensor produces a hexagonal tensor with vertical */
		/* symmetry axis. Rotate this to North/South                 */
		euler( 0.0, 90.0*DTOR, 0.0, rot );
		rotate3x3x3x3(c2,rot,c);

	for( h=min_h; h<=max_h; h+=d_h ) {
	for( a1=min_a1; a1<=max_a1; a1+=d_a1 ) {
	for( a2=min_a2; a2<=max_a2; a2+=d_a2 ) {

	Rall = 0.0;
	for( j=0; j<ntraces; j++ ) {

		vbackground=4500.0;
		a3=0.0;
		rho=3300.0; /* note scaling so result in Mbar */
		if ( !single_layer_anisotropy( bazi[j], 1000.0*vh[j], vbackground, h*1000.0, c, rho, a1, a2, a3,
        		ncoeff, ecoeff, zcoeff, radialcoeff, tangentialcoeff, time, flag ) ) {
			fprintf(stderr,"%s: error in call to single_layer_anisotropy\n", argv[0]);
			exit(-1);
			}

		a[0] = ncoeff[1];
		a[1] = ncoeff[0];

		b[0] = ecoeff[1];
		b[1] = ecoeff[0];

		clag[0]=0; clag[1]=(int) (0.5+((time[0]-time[1])/delta));

		n=head1[j].record.ndata;

		for( ii=0; ii<n; ii++ ) {
			data3[ii]=0.0; data4[ii]=0.0;
			for( jj=0; jj<2; jj++ ) {
				kk = ii-clag[jj];
				if( kk<0 ) kk+=n; /*circular*/
				data3[ii]+=b[jj]*data1[j][kk];
				data4[ii]+=a[jj]*data2[j][kk];
				}
			}

		E = 0.0; Ea = 0.0; Eb = 0.0;
		for( ii=0; ii<n; ii++ ) {
			Ea += data3[ii]*data3[ii];
			Eb += data4[ii]*data4[ii];
			E  += (data3[ii]-data4[ii]) * (data3[ii]-data4[ii]);
			}

		R = E/(Ea+Eb);

		Rall += R;

		} /* end for j */

		if( first ) {
			Rallmin=Rall;
			hmin=h; Bmin=B; a1min=a1; a2min=a2;
			first=FALSE;
			}
		else if( Rall<Rallmin ) {
			Rallmin=Rall;
			hmin=h; Bmin=B; a1min=a1; a2min=a2;
			first=FALSE;
			}
		if( !nullfile ) {
			fprintf(f3,"%f\t%f\t%f\t%f\t%f\n",h,B,a1,a2,Rall/(double)ntraces );
			}

		}}}} /*end for B, h, a1, a2 */

	fprintf(stdout,"h\tB\ta1\ta2\tR\n");
	fprintf(stdout,"%f\t%f\t%f\t%f\t%f\n",hmin,Bmin,a1min,a2min,Rallmin/(double)ntraces );

	h=hmin; B=Bmin; a1=a1min; a2=a2min; a3=0.0;

	/* setup tensor Vp=1.8*Vs, Vs=4500, B, C=0, E=B, rho=3300 */
	/* scale factor of 1e-11 gets you units of Mbar */
	Vs=4500.0;
	rho=3300.0;
	one_axis_tensor(1.8*Vs, Vs, B, 0.00, B, rho*1.0e-11, c );

	/* one_axis_tensor produces a hexagonal tensor with vertical */
	/* symmetry axis. Rotate this to North/South                 */
	euler( 0.0, 90.0*DTOR, 0.0, rot );
	rotate3x3x3x3(c2,rot,c);

	euler( a1*DTOR, a2*DTOR, a3*DTOR, rot );
	rotate3x3x3x3(c,rot,c2);

	/* write out tensor */
        fprintf(f6,"%.12f\n", rho );
        get_velocity( c2, rho, 90.0, 0.0, &v1, &v2, &v3);
        fprintf(f6,"%.12f %.12f %.12f\n", v1, v2, v3 );
        get_velocity( c2, rho, 90.0, 90.0, &v1, &v2, &v3);
        fprintf(f6,"%.12f %.12f %.12f\n", v1, v2, v3 );
        get_velocity( c2, rho, 0.0, 0.0, &v1, &v2, &v3);
        fprintf(f6,"%.12f %.12f %.12f\n", v1, v2, v3 );
        write3x3x3x3(f6,c2);

	Rall = 0.0;
	for( j=0; j<ntraces; j++ ) {

		vbackground=4500.0;
		a3=0.0;
		rho=3300.0;
		if ( !single_layer_anisotropy( bazi[j], 1000.0*vh[j], vbackground, h*1000.0, c, rho, a1, a2, a3,
        		ncoeff, ecoeff, zcoeff, radialcoeff, tangentialcoeff, time, flag ) ) {
			fprintf(stderr,"%s: error in call to single_layer_anisotropy\n", argv[0]);
			exit(-1);
			}

		a[0] = ncoeff[1];
		a[1] = ncoeff[0];;

		b[0] = ecoeff[1];
		b[1] = ecoeff[0];

		clag[0]=0; clag[1]=(int) (0.5+((time[0]-time[1])/delta));

		n=head1[j].record.ndata;

		for( ii=0; ii<n; ii++ ) {
			data3[ii]=0.0; data4[ii]=0.0;
			for( jj=0; jj<2; jj++ ) {
				kk = ii-clag[jj];
				if( kk<0 ) kk+=n; /*circular*/
				data3[ii]+=b[jj]*data1[j][kk];
				data4[ii]+=a[jj]*data2[j][kk];
				}
			}

		maxamp(&head1,data3);
        	if( xdr_putrecord(&head1,data3,&out) != 0) {
        		fprintf(stderr,"Error writing record %d %s\n",j,progname);
               		exit(-1);
               		}
		maxamp(&head1,data4);
       		if( xdr_putrecord(&head1,data4,&out) != 0) {
               		fprintf(stderr,"Error writing record %d %s\n",j,progname);
               		exit(-1);
               		}
		} /*end for j */

	exit(0);
	}




/* single_layer_anisotropy, calculate the arrival time and polarizations of
   a radially-polarized wave traversing an anisotropic layer

	input:
	bazi, backazimuth in degrees east of north
	vh, horizontal velocity in meters per second
	vbackground, background isotropic shear velocity in meters per second
	h, thickness of anisotropic layer in meters
	c2, standard 3x3x3x3 christoffel tensor
	rho, density in kilograms per cubic meter
	a1, a2, a3, standard euler angles for tensor in degrees
	flag, free surface correction included if TRUE

	output:
	n[2], amplitudes of north component, in meters (assuming unit initial polarization)
	e[2], amplitudes of east component, in meters
	v[2], amplitudes of vertical (positive up) component, in meters
	radial[2], amplitudes of radial-horizontal component, in meters
	tangential[2], amplitudes of tangential-horizontal component, in meters
	t[2], arrival times, in seconds

	returns: True (1) on no error, False (0) otherwise
*/
	
int single_layer_anisotropy( double bazi, double vh, double vbackground, double h,
	double c2[3][3][3][3], double rho, double a1, double a2, double a3,
        double *n, double *e, double *z, double *radial, double *tangential, double *t, int flag ) {

	int i, mode, status;
	double temp, ph, px, py, azi, angle_inc;
	double c[3][3][3][3];
	double rot[3][3];
	double szup[3][2], szdown[3][2], polup[3][3][2], poldown[3][3][2];
	double ctemp[2], slow[2][3][2], tan_in[3];
	double time[2], A[2][2], xcomp[2], ycomp[2], zcomp[2];
	double x[3], xi[3], s;
	double P[4], normal[3], initpol[3];
	double sazi, cazi;
	double slow_out[6][3][2], pol_out[6][3][2], A_out[6][2];


	/* top surface of anisotropic layer */
	x[0]=0.0; x[1]=0.0; x[2]=h;
	normal[0]=0.0; normal[1]=0.0; normal[2]=1.0;
	eqn_of_plane2( normal, x, P );

	/* below the anisotropic layer the plane S wave has polarization     */
	/* direction (CW from N) of bazi and horizontal phase velocity of vh */
	ph = 1.0/vh;
	azi=bazi+180.0; /* azimuth from backazimuth */
	azi=90.0-azi;   /* azimuth CCW from east from azimuth CW from N */
	if( azi > 360.0 ) azi-=360.0; else if (azi < 0.0) azi+=360.0;

	/* horizontal slowness determined by azimuth */
	sazi=sin(DTOR*azi); cazi=cos(DTOR*azi);
	px=ph*cazi;
	py=ph*sazi;

	/* initial polarization in background material*/
	angle_inc = asin( ph * vbackground );
	initpol[0]=cos(angle_inc)*cazi;
	initpol[1]=cos(angle_inc)*sazi;
	initpol[2]=sin(angle_inc);

	euler( a1*DTOR, a2*DTOR, a3*DTOR, rot );
	rotate3x3x3x3(c2,rot,c);

	/* calculate vertical slowness of plane wave of known horizontal slowness */ 
	status= slowness( c, rho, px, py, szup, szdown, polup, poldown );
	if( status != 1 ) {
		fprintf(stderr,"error: vertical slowness calculation failed\n");
		return(0);
		}

	for( mode=0; mode<2; mode++ ) {
		/* complete slowness vector of mode */
		slow[mode][0][0]=px;            slow[mode][0][1]=0.0;
		slow[mode][1][0]=py;            slow[mode][1][1]=0.0;
		slow[mode][2][0]=szup[mode][0]; slow[mode][2][1]=szup[mode][1];
		/* no evanescent waves allowed in this calculation */
		if( isevanescent(slow[mode]) ) {
			fprintf(stderr,"error: evanescent wave encountered\n");
			return(0);
			}

		/* project initial polarization onto mode to get its amplitude */
		A[mode][0]=0.0; A[mode][1]=0.0;
		for( i=0; i<3; i++ ) {
			cmul( polup[mode][i][0],polup[mode][i][1], initpol[i],0.0,
				&ctemp[0],&ctemp[1] );
			A[mode][0] += ctemp[0];
			A[mode][1] += ctemp[1];
			}

		/* ray starts at [x, y, z] = [0, 0, 0] */
		x[0]=0.0; x[1]=0.0; x[2]=0.0;

        	/* tangent to ray is parallel to slowness vector */
        	temp=sqrt( slow[mode][0][0]*slow[mode][0][0]
          		+  slow[mode][1][0]*slow[mode][1][0]
			+  slow[mode][2][0]*slow[mode][2][0] );
        	for(i=0; i<3; i++ ) tan_in[i]=slow[mode][i][0]/temp;

		/* take ray to top of anisotropic layer */
		if( !ray_exits_layer( P, x, tan_in, normal, xi, &s ) ) {
			fprintf(stderr,"error: ray doesnt intersect top of layer\n");
			return(0);
	       		}

		/* traveltime calculation */
		time[mode] = slow[mode][0][0]*(xi[0]-x[0])
			   + slow[mode][1][0]*(xi[1]-x[1])
			   + slow[mode][2][0]*(xi[2]-x[2]);

		xcomp[mode] = A[mode][0] * polup[mode][0][0];
		ycomp[mode] = A[mode][0] * polup[mode][1][0];
		zcomp[mode] = A[mode][0] * polup[mode][2][0];

		if (flag) {
			if( !scatter_fs( c, rho, c, 0.0, 1, slow[mode], polup[mode], A[mode],
	  			slow_out, pol_out, A_out) ) {
				fprintf(stderr,"error: free surface correction failed\n");
				exit(-1);
				}
			for(i=0; i<3; i++ ) {
				xcomp[mode] += A_out[i][0] * pol_out[i][0][0];
				ycomp[mode] += A_out[i][0] * pol_out[i][1][0];
				zcomp[mode] += A_out[i][0] * pol_out[i][2][0];
				}
			}

		/* (x, y, z) corresponds to (e, n, z) */
		n[mode]=ycomp[mode];
		e[mode]=xcomp[mode];
		z[mode]=zcomp[mode];
		t[mode]=time[mode];

                radial[mode] =     cazi*e[mode]    + sazi*n[mode];
                tangential[mode] = (-sazi*e[mode]) + cazi*n[mode];

		} /* next mode */


	return(1);
	}


/* euler rotation matrix, from Corbin & Stehle (1960) */
/* angles in radians */
/* rotation thru phi about 3-axis        */
/*          thru theta about new 1-axis  */
/*	    thru psi about new 3-axis    */

euler( phi, theta, psi, s )
double phi, theta, psi;
double s[3][3];
{
	s[0][0] = cos(psi)*cos(phi) - sin(psi)*cos(theta)*sin(phi);
	s[0][1] = cos(psi)*sin(phi) + sin(psi)*cos(theta)*cos(phi);
	s[0][2] = sin(psi)*sin(theta);

	s[1][0] = (-sin(psi)*cos(phi)-cos(psi)*cos(theta)*sin(phi));
	s[1][1] = (-sin(psi)*sin(phi)+cos(psi)*cos(theta)*cos(phi));
	s[1][2] = cos(psi)*sin(theta);

	s[2][0] = sin(theta)*sin(phi);
	s[2][1] = (-sin(theta)*cos(phi));
	s[2][2] = cos(theta);

	return(1);
	}

int eqn_of_plane( x1, x2, x3, P )
double x1[3], x2[3], x3[3];
double P[4];
{
	double x2mx1[3], x3mx1[3], v[3];
	int i;

	/* eqn is (x-x1) (dot) ( (x2-x1) (cross) (x3-x1) )     */
	/* see Rektorys Survey of Applicable Mathematics p 240 */

	for( i=0; i<3; i++ ) {
		x2mx1[i] = x2[i]-x1[i];
		x3mx1[i] = x3[i]-x1[i];
		}

	cross( x2mx1, x3mx1, v );
	P[0] = v[0];
	P[1] = v[1];
	P[2] = v[2];
	P[3] = (-dot(x1,v));
	return(1);
	}

int eqn_of_plane2( n, x1, P )
double n[3], x1[3];
double P[4];
{
	/* n = grad(Ax + By + Cz + D) = [A, B, C]  */
	/* then D = -(Ax1 + By1 + Cz1)
	/* see Rektorys Survey of Applicable Mathematics p 240 */

	P[0] = n[0];
	P[1] = n[1];
	P[2] = n[2];
	P[3] = (-(P[0]*x1[0]+P[1]*x1[1]+P[2]*x1[2]));
	return(1);
	}


/* dot product between vectors a and b */
double dot( a, b )
double a[3], b[3];
{
	return( a[0]*b[0] + a[1]*b[1] + a[2]*b[2] );
	}

/* cross product c = a (cross) b */
cross( a, b, c )
double a[3], b[3], c[3];
{
	c[0] = a[1]*b[2] - b[1]*a[2];
	c[1] = a[2]*b[0] - b[2]*a[0];
	c[2] = a[0]*b[1] - b[0]*a[1];
	return(1);
	}

/* finds intersection, xi, of ray with the plane, P. The ray is
   specified by a point, x, and a tangent unit vector t.
   The point must be below the plane, in the sense that the computed
   arclength, s, along the tangent from xi to the plane is positive.
   The normal of the plane is also calculated. Returns 1 on
   success, 0 on error. */

int ray_exits_layer( P, x, t, normal, xi, s )
double P[4], x[3], t[3];
double normal[3], xi[3], *s;
{
	int i;
	double a, xc[3];

	for( i=0; i<3; i++ ) xc[i]=x[i];

	/* calculate intersection of ray with plane, P.  The
	   intersection is found by simultaneously solving equation of plane,
	   Ax + By + Cz + D = 0 and position on ray as function of arclength, a,
	   x = xc[0]+a*t[0], y = xc[1]+a*t[1], z = xc[2]+a*t[2] */

	normal_to_plane( P, normal );
	a = (-(P[0]*xc[0]+P[1]*xc[1]+P[2]*xc[2]+P[3])/(P[0]*t[0]+P[1]*t[1]+P[2]*t[2]));
	for( i=0; i<3; i++ ) xi[i] = xc[i] + a*t[i];
	if( finite(a) && a > 0.0 ) { /* arclength has right sign */
		*s = a;
		return(1);
		}
	else {
		return(0);
		}

	}

/* normal to plane */
normal_to_plane( P, n )
double P[4], n[3];
{
	double a;

	/* normal is in direction grad(Ax+By+Cz+D)=(A,B,C) */
	n[0]=P[0], n[1]=P[1]; n[2]=P[2];
	a = sqrt(dot(n,n));
	n[0]/=a; n[1]/=a; n[2]/=a;
	return(1);
	}

int rotate3x3x3x3(c,s,c2)
double c[3][3][3][3];
double s[3][3];
double c2[3][3][3][3];
{
	int i, j, k, l;
	int p, q, r, t;

	for(i=0;i<3;i++) for(j=0;j<3;j++) for(k=0;k<3;k++) for(l=0;l<3;l++) {
		c2[i][j][k][l]=0.0;
		for(p=0;p<3;p++) for(q=0;q<3;q++) for(r=0;r<3;r++) for(t=0;t<3;t++) {
			c2[i][j][k][l] += s[i][p]*s[j][q]*s[k][r]*s[l][t]*c[p][q][r][t];
			}
		}
	return(1);
	}


/* complex multiplication: c = a*b */
int cmul( ar, ai, br, bi, cr, ci )
double ar, ai, br, bi;
double *cr, *ci;
{
	*cr = ar*br - ai*bi;
	*ci = ar*bi + ai*br;
	return(1);
	}

/* reads c, rho, v0, v1, v2 from file with given filename */
/* c: anistropic tensor, in Mbar           */
/* rho: density in km/m3                   */
/* v1: (vP, VS-fast, Vs-slow) along 1-axis */
/* v2: (vP, VS-fast, Vs-slow) along 1-axis */
/* v3: (vP, VS-fast, Vs-slow) along 1-axis */
/* returns 1 on success, 0 on error        */

int load3x3x3x3(filename, c, rho, v1, v2, v3)
char *filename;
double c[3][3][3][3];
double v1[3], v2[3], v3[3];
double *rho;
{
	FILE *f;
	int i, j, k, l, n;
	double t;
	
	if( (f=fopen(filename,"r")) == NULL ) return(0);

	if( fscanf(f,"%le", rho ) != 1 ) {fclose(f); return(0); }

	if( fscanf(f,"%le%le%le", &v1[0], &v1[1], &v1[2]) != 3 ) {fclose(f); return(0);}
	if( fscanf(f,"%le%le%le", &v2[0], &v2[1], &v2[2]) != 3 ) {fclose(f); return(0);}
	if( fscanf(f,"%le%le%le", &v3[0], &v3[1], &v3[2]) != 3 ) {fclose(f); return(0);}

	for( n=0; n<81; n++ ) {
		if( fscanf(f,"%d%d%d%d%le", &i, &j, &k, &l, &t) != 5 ) {fclose(f); return(0);}
		if( i<0 || i>2 ) {fclose(f); return(0);}
		if( j<0 || j>2 ) {fclose(f); return(0);}
		if( k<0 || k>2 ) {fclose(f); return(0);}
		if( l<0 || l>2 ) {fclose(f); return(0);}
		c[i][j][k][l] = t;
		}

	fclose(f);
	return(1);
	}

int isevanescent( slow )
double slow[3][2];
{
	int i;
	double s, cdot();

	/* a wave is evanescent if its slowness has a significant imaginary part */
	s = 1.0e-4 * sqrt(cdot(slow,slow));
	for( i=0; i<3; i++ ) if( fabs(slow[i][1]) > s ) return(1);
	return(0);
	}


/* finds the vertical slowness, sz (in s/m) and polarization, pol
for up and downgoing waves of a given mode (0, 1, 2) = (S-slow, S-fast, P)
of propagation into an anisotropic medium with tensor c (in Mbar), and
density rho (in kg/m3), given that the incident wave has horizontal
slowness (sx, sy) (in s/m). An upgoing propagating wave has real(sz)>0
and an upgoint (=decaying upward) evanscent wave has imag(sz)>0.
Returns 1 on success, 0 on failure (e.g no convergence to solution) */

int slowness( c, rho, sx, sy, szup, szdown, polup, poldown )
double c[3][3][3][3], rho, sx, sy;
double szup[3][2], szdown[3][2], polup[3][3][2], poldown[3][3][2];
{
	int i, j, p, q, k, n;
	int status[6];
	int na, na6, nt;

	double slow[3];
	double M[3][3][3];
	double f=1.0e11; /*conversion from Mbar to N/m2 */
	double d[7];
	double a[10][2], a6[10][2], b[2][2], t[10][2], x[2];
	double norm, scale, final_f[6];
	double root[6][2];
	double cslow[3][2];
	double size, dotprod;
	double disk;
	srandom(12475209);

	slow[0] = sx;
	slow[1] = sy;
	slow[2] = 0.0; /*never used*/

	/* matrix M(ip) = c(ijpq) s(j) s(q) - rho delta(ip) is quadradic in sz*/
	for( i=0; i<3; i++) for(p=0; p<3; p++ ) {
		for( k=0; k<=2; k++ ) M[i][p][k] = 0.0;
		if( i == p ) {
			M[i][p][0] -= rho;
			}
		for( j=0; j<3; j++ ) for( q=0; q<3; q++ ) {
			if( (j!=2) && (q!=2) ) { /*constant term*/
				M[i][p][0] += f*c[i][j][p][q] * slow[j] * slow[q];
				}
			else if( (j!=2) && (q==2) ) { /*linear term */
				M[i][p][1] += f*c[i][j][p][q] * slow[j];
				}
			else if( (j==2) && (q!=2) ) { /*linear term */
				M[i][p][1] += f*c[i][j][p][q] * slow[q];
				}
			else if( (j==2) && (q==2) ) { /*quadratic term */
				M[i][p][2] += f*c[i][j][p][q];
				}
			else {
				fprintf(stderr, "Slowness Internal Error!\n" );
				exit(-1);
				}
			}
		}


	detpoly(M,d);

	/*normalize det[M(sz)]=0 as follows: 
		switch from sz to scale*sz
		multiply by norm
		*/
	norm=1.0e-11;
	scale= 5000.0;
	for(k=0; k<=6; k++ ) {
		d[k] *= norm;
		norm /= scale;
		}
	
	/* find roots using the following algorithm:
		repeat:
			Newton'n method to find a root, x0
			reduce order by synthetic division by (x-x0)
			refine root against orgininal polynomial
		 */
	na=6;
	na6=6;
	for( i=0; i<=na; i++ ) {
		a[i][0]=d[i];
		a[i][1]=0.0;
		a6[i][0]=d[i];
		a6[i][1]=0.0;
		}
	for( k=0; k<6; k++ ) {
		disk=0.1;
		for( n=0; n<100; n++ ) {
			x[0]=disk*(myrand()-0.5);
			x[1]=disk*(myrand()-0.5);
			/* first search for a root in the reduced polynomial */
			if( (status[k]=findroot( a, na, x, 1.0e-10, &(final_f[0]))) == 1 ) {
				if( n>5 ) fprintf(stderr,"warning: findroot very slow\n");
				break;
				}
			}
		if(status[k] != 1) {
			fprintf(stderr,"scatter: coarse findroot failed (%e)\n",final_f[0]);
			fprintf(stderr,"the polynomial\n");
			for( n=0; n<=na; n++ ) {
				fprintf(stderr,"%d (%e, %e)\n", n, a6[n][0], a6[n][1] );
				}
			return(0);
			}
		/* but then refine it against the original polynomial */
		if( (status[k]=findroot( a6, na6, x, 1.0e-6, &(final_f[0]))) != 1 ) {
			fprintf(stderr,"slowness: refined findroot failed\n");

			return(0);
			}
		root[k][0]=x[0]/scale;
		root[k][1]=x[1]/scale;
		final_f[k]=final_f[0];
		if( k!=5 ) {
			b[0][0]=(-x[0]);
			b[0][1]=(-x[1]);
			b[1][0]=1.0;
			b[1][1]=0.0;
			monodiv( a, na, b, t, &nt );
			for( i=0; i<=nt; i++ ) {
				a[i][0]=t[i][0];
				a[i][1]=t[i][1];
				}
			na=nt;
			}
		}

	sortroots( root, szup, szdown );

	/* find corresponding polarizations using the equation
	         (C-ijpq s-j s-q - rho delta-ip) p-j = 0
								*/
	for(k=0; k<3; k++ ) {
		cslow[0][0] = slow[0];
		cslow[0][1] = 0.0;
		cslow[1][0] = slow[1];
		cslow[1][1] = 0.0;
		cslow[2][0] = szup[k][0];
		cslow[2][1] = szup[k][1];
		findp( c, rho, cslow, polup[k] );

		cslow[2][0] = szdown[k][0];
		cslow[2][1] = szdown[k][1];
		findp( c, rho, cslow, poldown[k] );
		}

	/* handle special degenerate case where both shear waves have
	   the same vertical slowness.  Then although the polarizations
	   of the two shear waves ought to be perpendicular, findp()
	   will not necessarily return perpendicular (or even different)
	   vectors */

	size = sqrt( (szup[0][0]-szup[1][0])*(szup[0][0]-szup[1][0]) +
		     (szup[0][1]-szup[1][1])*(szup[0][1]-szup[1][1]) ) /
		sqrt( szup[0][0]*szup[0][0] + szup[0][1]*szup[0][1] );
	if( size < 0.001 ) { /*assume degenerate*/

		i=0;
		/* find two polarizations that are not parallel */
		while( issamecvec(polup[0],polup[1]) ) {
			cslow[2][0] = szup[1][0];
			cslow[2][1] = szup[1][1];
			findp( c, rho, cslow, polup[1] );
			++i;
			if( i>10 ) {
				fprintf(stderr,"scatter: findp failed\n");
				return(0);
				}
			}
		orthogonalize( polup[0], polup[1] );
		}

	size = sqrt( (szdown[0][0]-szdown[1][0])*(szdown[0][0]-szdown[1][0]) +
		     (szdown[0][1]-szdown[1][1])*(szdown[0][1]-szdown[1][1]) ) /
		sqrt( szdown[0][0]*szdown[0][0] + szdown[0][1]*szdown[0][1] );
	if( size < 0.001 ) { /*assume degenerate*/
		i=0;
		/* find two polarizations that are not parallel */
		while( issamecvec(poldown[0],poldown[1]) ) {
			cslow[2][0] = szdown[1][0];
			cslow[2][1] = szdown[1][1];
			findp( c, rho, cslow, poldown[1] );
			++i;
			if( i>10 ) {
				return(0);
				}
			}
		orthogonalize( poldown[0], poldown[1] );
		}

	/* make sure polarizations are unit length */
	for( i=0; i<3; i++ ) {
		cnormalize( poldown[i] );
		cnormalize( polup[i] );
		}

	return(1);
	}

/*orthogonalize two complex vectors that define a plane so */
/* that they have Sv and Sh directionality                 */
int orthogonalize( a , b )
double a[3][2], b[3][2];
{
	int i, j;
	double t, max, xr, xi, yr, yi, c[3][2], d[3][2];

	/* c is normal to plane of a and b   */
	ccross( a, b, c );

	/* divide thru by largest element. a constant phase complex */
	/* vector is thus made real                                 */
	max=(-1.0); j=0;
	for( i=0; i<3; i++ ) {
		t = c[i][0]*c[i][0]+c[i][1]*c[i][1];
		if( max<t ) {
			max=t; j=i; yr=c[i][0]; yi=c[i][1];
			}
		}
	for( i=0; i<3; i++ ) {
		cdiv( c[i][0],c[i][1], yr,yi, &xr,&xi );
		c[i][0]=xr; c[i][1]=xi;
		}

	cnormalize(c);

	/* choose a and b such that a has Sv polarization and b has Sh polarization */
	if( (fabs(c[2][0])>0.9999) && (fabs(c[2][1])<0.0001) ) {
		/* a, b in horizontal plane */
		for(i=0;i<3;i++) for(j=0;j<2;j++) {
			a[i][j]=0.0;
			b[i][j]=0.0;
			}
		a[0][0]=1.0; b[1][0]=1.0;
		}
	else {
		/* d is vertical; SH: b=dxc and Sv: a=cxb */
		for(i=0;i<3;i++) for(j=0;j<2;j++) d[i][j]=0.0;
		d[2][0]=1.0;
		ccross( c, d, b );
		ccross( c, b, a );
		}
	return(1);
	}


sortroots( root, up, down )
double root[6][2], up[3][2], down[3][2];
{
	/* notes:
		A. upgoing is paired with exponentially decaying upward
			up:  + real or + imag
			dn:  - real or - imag
		B. modes ordered (0, 1, 2) = (S-slow, S-fast, P)
		*/

	int i, nup, ndown, cmproot();

	/* sort roots into supgoing modes and downgoing modes */
	nup=0; ndown=0;
	for(i=0; i<6; i++ ) {
		if( fabs(root[i][0]) >= fabs(root[i][1]) ) { /*propagating*/
			if( root[i][0] >= 0.0 ) {
				up[nup][0] = root[i][0];
				up[nup][1] = root[i][1];
				++nup;
				}
			else {
				down[ndown][0] = root[i][0];
				down[ndown][1] = root[i][1];
				++ndown;
				}
			}
		else { /*evanescent*/
			if( root[i][1] >= 0.0 ) {
				up[nup][0] = root[i][0];
				up[nup][1] = root[i][1];
				++nup;
				}
			else {
				down[ndown][0] = root[i][0];
				down[ndown][1] = root[i][1];
				++ndown;
				}
			}
		}

	/* sort slownesses within each directional set */
	qsort( up, 3, 2*sizeof(double), cmproot );
	qsort( down, 3, 2*sizeof(double), cmproot );
	}

int cmproot( a, b )
double a[2], b[2];
{
	/* compares phase velocities, ranking the slower ones
	   below the faster ones */
	int isap, isbp;

	if( fabs(a[0]) >= fabs(a[1]) ) isap=1; else isap=0; /*propagating or evanescent?*/
	if( fabs(b[0]) >= fabs(b[1]) ) isbp=1; else isbp=0;
	if( (isap==1) && (isbp==1) ) { /*propagating sorted by real part*/
		if( fabs(a[0]) > fabs(b[0]) ) return(-1);
		else if( fabs(a[0]) < fabs(b[0]) ) return(1);
		else return(0);
		}
	else if( (isap==0) && (isbp==0) ) { /*evanescent sorted by imag part*/
		if( fabs(a[1]) > fabs(b[1]) ) return(1);
		else if( fabs(a[1]) < fabs(b[1]) ) return(-1);
		else return(0);
		}
	else if( (isap==0) && (isbp==1) ) return(1); /*propagating sorted ahead of evanescent*/
	else return(-1);
	}


/*finds a root, x, of a polynomial equation with complex coefficients
	f(x) = a0 + a1 x + a2 x*x ... = 0
  by Newton's method. x must contain a starting guess
  iterations stop when |f|<=epsilon */
findroot( a, na, x, epsilon, final_f)
double a[10][2], x[2], epsilon, *final_f;
int na;
{
	int i;
	double t, dx[2], f[2], dfdx[2];

	for( i=0; i<50; i++ ) {
		evalpoly( a, na, x, f, dfdx );
		cdiv( f[0],f[1], dfdx[0],dfdx[1], &(dx[0]),&(dx[1]) );
		x[0] -= dx[0];
		x[1] -= dx[1];
		t = fabs(f[0]) + fabs(f[1]);
		if( (i>5) && (t<=epsilon) ) { *final_f=t; return(1); }
		}
	*final_f=t;
	return(0);
	}

/* evaluate polynomial with complex coefficients and calculate its derivative */
evalpoly( a, na, x, f, dfdx )
double a[10][2], x[2], f[2], dfdx[2];
int na;
{
	int i;
	double t[2], v[2];
	double h1[7][2], h2[7][2], h3[7][2], h4[7][2], h5[7][2];


	/* Horner's Scheme, see Rektorys p 60 */
	for( i=0; i<=na; i++ ) {
		h1[na-i][0]=a[i][0];
		h1[na-i][1]=a[i][1];
		}
	h2[0][0]=0.0; h2[0][1]=0.0;
	h3[0][0]=h1[0][0]; h3[0][1]=h1[0][1];
	for( i=1; i<=na; i++ ) {
		/* h2[i] = x*h3[i-1] */
		cmul( x[0],x[1], h3[i-1][0],h3[i-1][1], &(h2[i][0]),&(h2[i][1]) );
		/* h3[i] = h1[i]+h2[i] */
		h3[i][0] = h1[i][0]+h2[i][0];
		h3[i][1] = h1[i][1]+h2[i][1];
		}
	h4[0][0]=0.0; h4[0][1]=0.0;
	h5[0][0]=h3[0][0]; h5[0][1]=h3[0][1];
	for( i=1; i<na; i++ ) {
		cmul( x[0],x[1], h5[i-1][0],h5[i-1][1], &(h4[i][0]),&(h4[i][1]) );
		h5[i][0] = h3[i][0]+h4[i][0];
		h5[i][1] = h3[i][1]+h4[i][1];
		}

	f[0]=h3[na][0]; f[1]=h3[na][1];
	dfdx[0]=h5[na-1][0]; dfdx[1]=h5[na-1][1];

	return(1);
	}



/* monomial division
   input: complex polynomial a (of order na)
          complex monomial b
   returned: complex polynomial c=a/b (of order nc) */
monodiv( a, na, b, c, nc )
double a[10][2], b[2][2], c[10][2];
int na, *nc;
{
	int i;
	double t[2];

	*nc = na-1;
	cdiv( a[0][0],a[0][1], b[0][0],b[0][1], &(c[0][0]), &(c[0][1]) );

	for( i=1; i<=(*nc); i++ ) {
		cmul( b[1][0],b[1][1], c[i-1][0],c[i-1][1], &(t[0]), &(t[1]) );
		t[0] = a[i][0]-t[0];
		t[1] = a[i][1]-t[1];
		cdiv( t[0],t[1], b[0][0],b[0][1], &(c[i][0]), &(c[i][1]) );
		}

	return(1);
	}

/* complex division: c = a/b */
int cdiv( ar, ai, br, bi, cr, ci )
double ar, ai, br, bi;
double *cr, *ci;
{
	double t;
	t = br*br + bi*bi;
	*cr = (ar*br + ai*bi)/t;
	*ci = (ai*br - ar*bi)/t;
	return(1);
	}

int findp( c, rho, s, pol )
double c[3][3][3][3], rho, s[3][2], pol[3][2];
{
	int i, j, p, q, itriag, ierror, tries;
	int k;
	double a, br, bi, sum;
	double M[3][3][2], Mc[3][3][2], T[3][3][2];
	double t[2], newpol[3][2];

	/* use the 'inverse power method' to find the polarizations. It
	   is based on the iterating the equation

		M-ip * pol(new)-p =  pol(old)-i

	Since renormaization occurs, we only need Minverse up to a constant
	multiplicative factor and can use just the matrix of minors */

	/* compute M-ip = C-ijpq s-j s-q - rho delta-ip for s perturbed by ep*/
	for( i=0; i<3; i++ ) for( p=0; p<3; p++ ) {
		if( i==p ) {
			M[i][p][0] = (-rho);
			M[i][p][1] = 0.0;
			}
		else {
			M[i][p][0] = 0.0;
			M[i][p][1] = 0.0;
			}
		for( j=0; j<3; j++ ) for( q=0; q<3; q++ ) {
			cmul( s[j][0],s[j][1], s[q][0],s[q][1], &(t[0]), &(t[1]) );
			M[i][p][0] += 1.0e11 * c[i][j][p][q] * t[0];
			M[i][p][1] += 1.0e11 * c[i][j][p][q] * t[1];
			}
		}

	for( i=0; i<3; i++ ) for( p=0; p<3; p++ ) {
		Mc[i][p][0]=M[i][p][0];
		Mc[i][p][1]=M[i][p][1];
		}
	
	tries=0;
	do { /*this loop handles failure of proceedure due to singular matrix*/
		/* some random starting polarizations to seed the iterative calculation */
		for(i=0;i<3;i++) for(j=0;j<2;j++) {
			pol[i][j]=2.0*(myrand()-0.5);
			}

		/* iterate a few times */
		itriag=1;
		for( i=0; i<3; i++ ) for( p=0; p<3; p++ ) { /* perturb M slightly   */
							    /* to avoid singularity */
			M[i][p][0]=(1.0+1.0e-10*myrand())*Mc[i][p][0];
			M[i][p][1]=(1.0+1.0e-10*myrand())*Mc[i][p][1];
			}
		for( q=0; q<10; q++ ) {

			/* M-ip pol(new)-p = Minverse-ip pol(old)-i */
			for( i=0; i<3; i++ ) {
				newpol[i][0]=pol[i][0];
				newpol[i][1]=pol[i][1];
				}
			cgauss(M, newpol, 3, 3, 0.0, &ierror, itriag );
			itriag=0;
			ierror=0;

			/*renormalize to a unit vector*/
			sum=0.0;
			for( i=0; i<3; i++ ) {
				sum += newpol[i][0]*newpol[i][0] + newpol[i][1]*newpol[i][1];
				}
			sum = sqrt(sum);
			for( i=0; i<3; i++ ) {
				pol[i][0] = newpol[i][0]/sum;
				pol[i][1] = newpol[i][1]/sum;
				}

			for( i=0; i<3; i++ ) {
				if( !finite(pol[i][0]) || !finite(pol[i][1]) ) {
					ierror=1;
					}
				}
			if( ierror!=0 ) {
				break;
				}
			}
		tries++;
		} while( (ierror!=0) && (tries<10) );

	if( ierror!=0 ) {
		fprintf(stderr,"iteration fails in findp\n");
		return(0);
		}

	/* normalize so that largest element has zero phase */
	/* this beautifies the result, removing an overall complex phase */
	q=0;
	sum=(-1.0);
	for( i=0; i<3; i++ ) {
		a=sqrt( pol[i][0]*pol[i][0] + pol[i][1]*pol[i][1] );
		if( a>sum ) {
			q=i;
			sum=a;
			}
		}
	cdiv( sum, 0.0, pol[q][0], pol[q][1], &(t[0]), &(t[1]) );
	for( i=0; i<3; i++ ) {
		cmul( pol[i][0], pol[i][1], t[0], t[1], &(pol[i][0]), &(pol[i][1]) );
		}

	/*checkpolar(Mc,pol);*/

	return(1);
	}

/* complex dot product */
double cdot( a, b )
double a[3][2], b[3][2];
{
	int k;
	double dot;

	dot = 0.0;
	for( k=0; k<3; k++ ) {
		dot += a[k][0]*b[k][0] + a[k][1]*b[k][1];
		}

	return(dot);
	}

/* comples cross product c = a (cross) b*/
int ccross( a, b, c )
double a[3][2], b[3][2], c[3][2];
{
	double xr, xi, yr, yi;
	cmul( a[1][0],a[1][1], b[2][0],b[2][1], &xr,&xi);
	cmul( a[2][0],a[2][1], b[1][0],b[1][1], &yr,&yi);
	c[0][0]=xr-yr;
	c[0][1]=xi-yi;
	cmul( a[2][0],a[2][1], b[0][0],b[0][1], &xr,&xi);
	cmul( a[0][0],a[0][1], b[2][0],b[2][1], &yr,&yi);
	c[1][0]=xr-yr;
	c[1][1]=xi-yi;
	cmul( a[0][0],a[0][1], b[1][0],b[1][1], &xr,&xi);
	cmul( a[1][0],a[1][1], b[0][0],b[0][1], &yr,&yi);
	c[2][0]=xr-yr;
	c[2][1]=xi-yi;
	return(1);
	}

double myrand()
{
	double x;
	x=( ((double)(random())) / ((double)2147483647) );
	return(x);
	}

/* normalize complex vector to unit length */
int cnormalize( a )
double a[3][2];
{
	int i, j;
	double t;
	t = sqrt(cdot(a,a));
	for(i=0;i<3;i++) for(j=0;j<2;j++) a[i][j] /= t;
	return(1);
	}

#define MAX_TABLE_COLS (100)

cgauss(a,vec,n,nstore,test,ierror,itriag)
double *a, *vec, test;
int n, nstore, *ierror, itriag;
{
 
/* subroutine cgauss, by William Menke, 9-23-97 */
/* based on gauss, written by Menke many years ago */
 
/* a subroutine to solve a system of n linear equations with complex */
/* coefficients in n unknowns, where n doesn't exceed MAX_TABLE_COLS */
/* gaussian reduction with partial pivoting is used                  */
/* 	a		(sent, destroyed)	n by n matrix		        */
/*	vec		(sent, overwritten)	n vector, replaced w/ solution  */
/* 	nstore		(sent)			dimension of a	                */
/*	test		(sent)			div by zero check number        */
/*	ierror		(returned)		zero on no error                */
/*	itriag		(sent)			matrix triangularized only      */
/*						 on TRUE useful when solving    */
/*						 multiple systems with same a   */
	static int isub[MAX_TABLE_COLS], l1;
	int line[MAX_TABLE_COLS], iet, ieb, i, j, k, l, j2;
	double big, testa, br, bi, sumr, sumi, *offset, rp, ip, len;
	

	iet=0;	/* initial error flags, one for triagularization*/
	ieb=0;	/* one for backsolving */

/* triangularize the matrix a*/
/* replacing the zero elements of the triangularized matrix */
/* with the coefficients needed to transform the vector vec */

	if (itriag) {	/* triangularize matrix */
 
 		for( j=0; j<n; j++ ) {      /*line is an array of flags*/
			line[j]=0; 
			/* elements of a are not moved during pivoting*/
			/* line=0 flags unused lines */
			}    /*end for j*/
			
		for( j=0; j<n-1; j++ ) {
			/*  triangularize matrix by partial pivoting */
                       big = 0.0; /* find biggest element in j-th column*/
				  /* of unused portion of matrix*/
		       for( l1=0; l1<n; l1++ ) {
                               if( line[l1]==0 ) {
				       offset = a + 2*(l1*nstore+j);
				       rp = *offset; ip = *(offset+1);
				       testa = sqrt( rp*rp + ip*ip );
                                       if (testa>big) {
						i=l1;
						big=testa;
						} /*end if*/
					} /*end if*/
				} /*end for l1*/
                       if( big<=test) {   /* test for div by 0 */
                               iet=1;
                               } /*end if*/
 
                       line[i]=1;  /* selected unused line becomes used line */
                       isub[j]=i;  /* isub points to j-th row of tri. matrix */
 
		       offset = a+2*(i*nstore+j);
		       rp = *offset; ip = *(offset+1);
		       len = rp*rp + ip*ip;
		       sumr =  rp/len;
		       sumi = -ip/len;
			/*reduce matrix towards triangle */
		       for( k=0; k<n; k++ ) {
				if( line[k]==0 ) {
					offset = a+2*(k*nstore+j);
		       			rp = *offset; ip = *(offset+1);
					br = rp*sumr - ip*sumi;
					bi = rp*sumi + ip*sumr;
				       	for( l=j+1; l<n; l++ ) {
					       offset = a+2*(i*nstore+l);
					       rp = *offset; ip = *(offset+1);
					       offset = a+2*(k*nstore+l);
                                               *offset -= (br*rp-bi*ip);
					       *(offset+1) -= (br*ip+bi*rp);
					       } /*end for l*/
				       offset = a+2*(k*nstore+j);
				       *offset = br;
				       *(offset+1) = bi;
					} /*end if*/
				} /*end for k*/
			} /*end for j*/
 
               for( j=0; j<n; j++ ) {
			/*find last unused row and set its pointer*/
			/*  this row contians the apex of the triangle*/
			if( line[j]==0) {
				l1=j;   /*apex of triangle*/
				isub[n-1]=j;
				break;
				} /*end if*/
			} /*end for j*/
 
		} /*end if itriag true*/
		
	/*start backsolving*/
	
	for( i=0; i<n; i++ ) {	/* invert pointers. line(i) now gives*/
				/* row no in triang matrix of i-th row*/
				/* of actual matrix */
		line[isub[i]] = i;
		} /*end for i*/
 
	for( j=0; j<n-1; j++) { /*transform the vector to match triang. matrix*/
	       offset = vec+2*isub[j];
               br = *offset;
	       bi = *(offset+1);
               for( k=0; k<n; k++ ) {
                      if (line[k]>j) {	/* skip elements outside of triangle*/
				offset = a+2*(k*nstore+j);
				rp = *offset; ip = *(offset+1);
				offset = vec + 2*k;
				*offset -= (br*rp - bi*ip );
				*(offset+1) -= (br*ip + bi*rp );
				} /*end if*/
			} /*end for k*/
		} /*end for j*/
 
      offset = a+2*(l1*nstore+(n-1));
      br= *offset; bi = *(offset+1); /*apex of triangle*/
      if( sqrt(br*br+bi*bi)<=test) {
		/*check for div by zero in backsolving*/
		ieb=2;
		} /*end if*/
      offset = vec+2*isub[n-1];
      rp = *offset; ip = *(offset+1);
      len = br*br + bi*bi;
      *offset = (rp*br + ip*bi)/len;
      *(offset+1) = (ip*br-rp*bi)/len;
 
      for( j=n-2; j>=0; j-- ) {	/* backsolve rest of triangle*/
		offset = vec+2*isub[j];
		sumr = *offset; sumi = *(offset+1);
		for( j2=j+1; j2<n; j2++ ) {
			offset = a+2*(isub[j]*nstore+j2);
			rp = *offset; ip = *(offset+1);
			offset = vec+2*isub[j2];
			sumr -= ( (*(offset))*rp - (*(offset+1))*ip );
			sumi -= ( (*(offset))*ip + (*(offset+1))*rp );
			} /*end for j2*/
		offset = a+2*(isub[j]*nstore+j);
		br = *offset; bi=*(offset+1);
               if( sqrt(br*br+bi*bi) <= test) {
			/* test for div by 0 in backsolving */
			ieb=2;
			} /*end if*/
		offset = vec+2*isub[j];
		len = br*br + bi*bi;
		*offset = (sumr*br+sumi*bi)/len;
		*(offset+1) = (sumi*br-sumr*bi)/len; /*solution returned in vec*/
		} /*end for j*/

/*put the solution vector into the proper order*/

      for( i=0; i<n; i++ ) {    /* reorder solution */
		for( k=i; k<n; k++ ) {  /* search for i-th solution element */
			if( line[k]==i ) {
				j=k;
				break;
				} /*end if*/
			} /*end for k*/
       	       /* swap solution and pointer elements*/
               br = *(vec+2*j); bi = *(vec+2*j+1);
               *(vec+2*j) = *(vec+2*i); *(vec+2*j+1) = *(vec+2*i+1);
               *(vec+2*i) = br; *(vec+2*i+1) = bi;
               line[j] = line[i];
		} /*end for i*/
 
      *ierror = iet + ieb;   /* set final error flag*/
}


/* are complex vector a and b the same ? */
int issamecvec( a, b )
double a[3][2], b[3][2];
{
	int i, j;
	double t1, t2;

	t1 = sqrt(cdot(a,a));
	t2 = sqrt(cdot(a,a));
	if( t2>t1 ) t1=t2;
	t1 *= 1.0e-4;
	for( i=0;i<3;i++) for( j=0;j<2;j++ ) if( fabs(a[i][j]-b[i][j])>t1 ) return(0);
	
	return(1);
	}


/* determinant of 3x3 matrix, each element of which is a quadratic
polynomial with real coefficients*/
int detpoly( M, d )
double M[3][3][3], d[7];
{
	int k;

	for( k=0; k<=6; k++ ) d[k]=0.0;

	sumdetpoly( M, d, 0,0, 1,1, 2,2, 1.0 );
	sumdetpoly( M, d, 0,1, 1,2, 2,0, 1.0 );
	sumdetpoly( M, d, 0,2, 2,1, 1,0, 1.0 );
	sumdetpoly( M, d, 0,2, 1,1, 2,0, -1.0 );
	sumdetpoly( M, d, 0,1, 1,0, 2,2, -1.0 );
	sumdetpoly( M, d, 0,0, 2,1, 1,2, -1.0 );
	return(1);
	}

/* adds one term, a 6-th order polynomial, to determinant */
sumdetpoly( M, d, i,j, k,l, m,n, s )
double M[3][3][3], d[7];
int i,j, k,l, m,n;
double s;
{
	int p, q, r;
	for(p=0; p<=2; p++) for(q=0; q<=2; q++ ) for( r=0; r<=2; r++ ) {
		d[p+q+r] += s * M[i][j][p] * M[k][l][q] * M[m][n][r];
		}
	return(1);
	}

int one_axis_tensor(double Vp, double Vs, double B, double C, double E, double rho, double c[3][3][3][3] ) {
	double A, D, lc[6][6];

        A = Vp*Vp*rho;
        B = A*B;
        C = A*C;
        D = Vs*Vs*rho;
        E = D*E;

        ABCDEtolc( A, B, C, D, E, lc );
        lctobc( lc, c );

	return(1);
	}

/* 6x6 to 3x3x3x3 tensor transformation */
int lctobc( lc, bc )
double lc[6][6];
double bc[3][3][3][3];
{
        int i, j, m, n, p, q;

        /* from Fuchs, K, Phys. Earth Planet. Int. 31, 93-118, 1983 */

        for(m=0;m<3;m++) for(n=0;n<3;n++) for(p=0;p<3;p++) for(q=0;q<3;q++) bc[m][n][p][q] = 0.0;

        for(m=0;m<3;m++) for(n=0;n<3;n++) for(p=0;p<3;p++) for(q=0;q<3;q++) {
                if( m==n ) i=(m+n+2)/2; else i=9-m-n-2;
                if( p==q ) j=(p+q+2)/2; else j=9-p-q-2;
                bc[m][n][p][q] = lc[i-1][j-1];
                }

        return(1);
        }

/* Anisotropy with one axis of symmetry, from Vadim's write-up */
int ABCDEtolc( A, B, C, D, E, lc )
double A, B, C, D, E, lc[6][6];
{
        int i, j;

        for( i=0; i<6; i++ ) for( j=0; j<6; j++ ) lc[i][j]=0.0;

        lc[0][0] = A - B + C;
        lc[0][1] = lc[1][0] = A - B + C - 2.0*(D-E);
        lc[0][2] = lc[2][0] = lc[2][1] = lc[1][2] = A - 3.0*C - 2.0*(D+E);
        lc[1][1] = A - B + C;
        lc[2][2] = A + B + C;
        lc[3][3] = lc[4][4] = D + E;
        lc[5][5] = D - E;

        return(1);
        }

int euler_blerb() {
fprintf(stderr,"euler angles, a1, a2, a3 in degrees (see Corbin & Stehle, 1960):\n");
fprintf(stderr,"	CCW rotation thru a1 about 3-axis\n");
fprintf(stderr,"	             thru a2 about new 1-axis\n");
fprintf(stderr,"	             thru a3 about new 3-axis\n");
fprintf(stderr,"        note that these angles rotate the coordinate system, not\n");
fprintf(stderr,"                 the object! So a a1 of 30 deg rotates an object in that\n");
fprintf(stderr,"		 coordinate system by -30 deg.\n");
return(1);
}

int write3x3x3x3(f,bc)
FILE *f;
double bc[3][3][3][3];
{
        int m, n, p, q;

        for(m=0;m<3;m++) for(n=0;n<3;n++) for(p=0;p<3;p++) for(q=0;q<3;q++) {
                fprintf(f, "%d %d %d %d %.12f\n", m, n, p, q, bc[m][n][p][q] );
                }
        return(1);
        }


/* input: media parameters c (Mbar), rho (kg/m2)
          propagation directions (theta, phi) in degrees
   output: quasi-P, 2 quasi-S phase velocities            */

get_velocity( c, rho, theta, phi, vp, vs1, vs2 )
double c[3][3][3][3], rho;
double theta, phi;
double *vp, *vs1, *vs2;
{
	double st[3];
	double MM[3][3];
	double b[6];
	double lambda[3], vec[3][3];
	double pp;

	/* make a unit vector pointing in (theta, phi) direction */
	tangent( theta, phi, st );

	/* make the 3x3 matrix MM*/
	do_MM( c, st, MM );

	/* compute its eigenvalues and vectors */
        b[0]=MM[0][0]; b[1]=MM[0][1]; b[2]=MM[0][2];
	b[3]=MM[1][1]; b[4]=MM[1][2]; b[5]=MM[2][2];
        eigen3( b, lambda, vec );

	*vp =  sqrt(lambda[2]/rho);
	*vs1 = sqrt(lambda[1]/rho);
	*vs2 = sqrt(lambda[0]/rho);

	}

do_MM( c, t, MM )
double c[3][3][3][3], t[3], MM[3][3];
{
	int i, j, p, q;

	/* note factor of 1e11, conversion from Mbar to N/m2 */

	for( i=0; i<3; i++) for(p=0; p<3; p++ ) {
		MM[i][p] = 0.0;
		for( j=0; j<3; j++ ) for( q=0; q<3; q++ ) {
			MM[i][p] += 1.e11*c[i][j][p][q] * t[j] * t[q];
			}
		}
	}

int eigen3( double b[6], double eval[3], double evec[3][3] ) {

/*
c
c subroutine eigen3, by william menke, march, 1980.
c Fortran with Semicolons C conversion, January 2003.
c
c this subroutine computes the eigenvalues and eigenvectors
c of a real, symmetric, three by three matrix. note that if
c the matrix has several equal eigenvalues, the corresponding
c eigenvectors may not be orthogonal. parameters :
c
c b ...... a vector containing the independent elements of
c          the matrix, say m, in the order :
c
c          b =   m  , m  , m  , m  , m  , m
c                 00   01   02   11   12   22
c
c eval ... the three eigenvalues, smallest first.
c
c evec ... a matrix whose columns are the three eigenvectors.
c
*/
	int i, j, itt, ivec;
        double y[3][3], temp[3];
	double p, q, r, ax, bx, xm, temp2;
	double th1, th2, th3, t1, t2, xl, sum;
	double test, maxval, diff1, diff2;

        p = - ( b[0] + b[3] + b[5] );
        q = b[0]*b[3] + b[0]*b[5] + b[3]*b[5]
          - b[2]*b[2] - b[1]*b[1] - b[4]*b[4];
        r = b[5]*b[1]*b[1] + b[0]*b[4]*b[4]    + b[3]*b[2]*b[2]
     *                     - 2.*b[1]*b[2]*b[4] - b[0]*b[3]*b[5];

        ax = (3.0*q-p*p)/3.0;
        bx = (2.0*p*p*p - 9.0*p*q + 27.0*r)/27.0;

	temp2 = -ax/3.0;
	if( temp2<0.0 ) temp2=0.0;
        xm = 2. * sqrt( temp2 );
	temp2=3.*bx/(ax*xm);
	if( temp2>1.0 ) {
		temp2=1.0;
		}
	else if( temp2<(-1.0) ) {
		temp2=(-1.0);
		}	
        th1 = acos( temp2 ) / 3.0;
        th2 = th1 + 0.6666666667*3.141592654;
        th3 = th1 + 1.3333333333*3.141592654;

        eval[0] = xm * cos(th1) - p/3.0;
        eval[1] = xm * cos(th2) - p/3.0;
        eval[2] = xm * cos(th3) - p/3.0;

        for( j=0; j<3; j++ ) {
       	 	for( i=0; i<2; i++ ) {
        		t1 = eval[i];
        		t2 = eval[i+1];
       			if( t1 > t2 ) {
        			eval[i]=t2;
        			eval[i+1]=t1;
				}
			}
		}

        evec[0][0] = 0.64638;
        evec[1][0] = 0.16349;
        evec[2][0] = 0.46921;
        evec[0][1] = 0.27336;
        evec[1][1] =-0.88332;
        evec[2][1] = 0.36313;
        evec[0][2] =-0.99843;
        evec[1][2] =-0.26391;
        evec[2][2] = 0.72759;

        for( ivec=0; ivec<3; ivec++ ) {

		xl = 1.00000001*eval[ivec];
		y[0][0] = (b[3]-xl)*(b[5]-xl) - b[4]*b[4];
		y[1][1] = (b[0]-xl)*(b[5]-xl) - b[2]*b[2];
		y[2][2] = (b[0]-xl)*(b[3]-xl) - b[1]*b[1];
		y[0][1] = b[2]*b[4] - (b[5]-xl)*b[1];
		y[0][2] = b[1]*b[4] - (b[3]-xl)*b[2];
		y[1][2] = b[1]*b[2] - (b[0]-xl)*b[4];
		y[1][0] = y[0][1];
		y[2][0] = y[0][2];
		y[2][1] = y[1][2];

		for( itt=0; itt<3; itt++ ) {

			for( i=0; i<3; i++ ) {
				sum = 0.0;
				for( j=0; j<3; j++ ) {
					sum += y[i][j]*evec[j][ivec];
					}
        			temp[i] = sum;
				}

        		for( i=0; i<3; i++ ) {
				evec[i][ivec]=temp[i];
				}

        		sum = 0.0;
        		for( i=0; i<3; i++ ) {
       				test = fabs( evec[i][ivec] );
        			if(test>sum) sum=test;
				}
       			if( sum==0.0 ) sum=1.0;

        		for( i=0; i<3; i++ ) {
        			evec[i][ivec] /= sum;
				}

			}

       		sum=0.0;
       		for( i=0; i<3; i++ ) {
        		sum += evec[i][ivec]*evec[i][ivec];
			}
        	sum = sqrt(sum);
        	if(sum==0.0) sum=1.0;
        	for( i=0; i<3; i++ ) {
        		evec[i][ivec] /= sum;
			}

		}

	maxval=fabs(eval[0]);
	test=fabs(eval[1]);
	if(test>maxval) maxval=test;
	test=fabs(eval[2]);
	if(test>maxval) maxval=test;
	if( maxval<=0.0 ) maxval=1.0;
	diff1 = fabs(eval[0]-eval[1])/maxval;
	diff2 = fabs(eval[1]-eval[2])/maxval;

	if( (diff1<=1.0e-8) && (diff2<=1.0e-8) ) {
		evec[0][0] = 1.0;
		evec[1][0] = 0.0;
		evec[2][0] = 0.0;
		evec[0][1] = 0.0;
		evec[1][1] = 1.0;
		evec[2][1] = 0.0;
		evec[0][2] = 0.0;
		evec[1][2] = 0.0;
		evec[2][2] = 1.0;
		}
	else if( (diff1<=1.0e-8) && (diff2>1.0e-8) ) {
		evec[0][0] = evec[1][1]*evec[2][2]-evec[2][1]*evec[1][2];
		evec[1][0] = evec[2][1]*evec[0][2]-evec[0][1]*evec[2][2];
		evec[2][0] = evec[0][1]*evec[1][2]-evec[1][1]*evec[0][2];
		}
	else if( (diff1>1.0e-8) && (diff2<=1.0e-8) ) {
		evec[0][2] = evec[1][0]*evec[2][1]-evec[2][0]*evec[1][1];
		evec[1][2] = evec[2][0]*evec[0][1]-evec[0][0]*evec[2][1];
		evec[2][2] = evec[0][0]*evec[1][1]-evec[1][0]*evec[0][1];
		}
	
        return(1);
        }

int tangent( theta, phi, t )
double theta, phi, t[3];
{
        t[0] = sin( DTOR * theta ) * cos( DTOR * phi );
        t[1] = sin( DTOR * theta ) * sin( DTOR * phi );
        t[2] = cos( DTOR * theta );
        }

int tangentr( theta, phi, t )
double theta, phi, t[3];
{
        t[0] = sin( theta ) * cos( phi );
        t[1] = sin( theta ) * sin( phi );
        t[2] = cos( theta );
        }

int viresp( double baz, double theta, double a[2], double b[2]) {

	double ct, st, cbaz, sbaz, ap[2], bp[2];
	int jj;

	ct = cos(DTOR*(baz-theta));
	st = sin(DTOR*(baz-theta));

	cbaz = cos(DTOR*baz);
	sbaz = sin(DTOR*baz);

        ap[0] = ct*ct;
        ap[1] = st*st;

        bp[0] = ct*st;
        bp[1] = (-bp[0]);

	for( jj=0; jj<2; jj++ ) {
        	a[jj] = cbaz*ap[jj] + sbaz*bp[jj];
      		b[jj] = sbaz*ap[jj] - cbaz*bp[jj];
                }

	return(1);
	}

/* reflection and transmission coefficients for a plane wave incident
from below (dir_in = 1) or above (dir_in = -1) upon a horizontal interface
between two anisotropic halfspaces.  The halfspaces have real anisotropy
tensors c_below and c_above (in Mbar) and real density rho_below and
rho_above (in kg/m3).  The incident wave has complex slowness slow_in
= [(sx, 0.0), (sy,0.0), (sz-r, sz-i)]-transpose (in s/m), complex
polarization unit vector pol_in, and complex amplitude A_in. Corresponding
values for the scattered waves are arranged as:

	0: slow_out[0][i], pol_out[0][i], A_out[0]: S-slow-down
	1: slow_out[1][i], pol_out[1][i], A_out[1]: S-fast-down
	2: slow_out[2][i], pol_out[2][i], A_out[2]: P-down

	3: slow_out[3][i], pol_out[3][i], A_out[3]: S-slow-up
	4: slow_out[4][i], pol_out[4][i], A_out[4]: S-fast-up
	5: slow_out[5][i], pol_out[5][i], A_out[5]: P-up

returns 1 on success, 0 on failure */
/* like scatterss, but treats the medium with rho=0.0 as vacuum */
scatter_fs( c_below, rho_below, c_above, rho_above, dir_in, slow_in, pol_in, A_in,
	  slow_out, pol_out, A_out)
double c_below[3][3][3][3], c_above[3][3][3][3];
double rho_below, rho_above;
int dir_in;
double slow_in[3][2], pol_in[3][2], A_in[2];
double slow_out[6][3][2], pol_out[6][3][2], A_out[6][2];
{
	int i, j, k, wave, p, q, ierror;
	double normal[3];
	double szup[3][2], szdown[3][2], polup[3][3][2], poldown[3][3][2];
	double M[3][3][2], x[3][2], t[2];

	if( (dir_in!=(1)) && (dir_in!=(-1)) ) {
		fprintf(stderr,"scatter_fs: bad value of dir_in\n");
		return(0);
		}
	else if( (rho_above==0.0) && (rho_below==0.0) ) {
		fprintf(stderr,"scatter_fs: two vacuums in contact\n");
		return(0);
		}
	else if( (rho_below==0.0) && (dir_in==(1)) ) {
		fprintf(stderr,"scatter_fs: upgoing wave in vacuum\n");
		return(0);
		}
	else if( (rho_above==0.0) && (dir_in==(-1)) ) {
		fprintf(stderr,"scatter_fs: downgoing wave in vacuum\n");
		return(0);
		}
	else if( (rho_above!=0.0) && (rho_below!=0.0) ) {
		fprintf(stderr,"scatter_fs: no free surface\n");
		return(0);
		}

	/* slowness and polarizations of scattered waves in bottom media*/
	if( dir_in < 0 ) {
		for( i=0; i<3; i++ ) for( j=0; j<2; j++ ) {
			slow_out[i][0][j] = 0.0;
			slow_out[i][1][j] = 0.0;
			slow_out[i][2][j] = 0.0;
			}
		for( i=0; i<3; i++ ) for( j=0; j<3; j++ ) for( k=0; k<2; k++ ) {
			pol_out[i][j][k] = 0.0;
			}
		}
	else {
		if( slowness( c_below, rho_below, slow_in[0][0], slow_in[1][0],
			szup, szdown, polup, poldown ) != 1 ) {
			fprintf(stderr,"slowness (below) fails in scatter\n");
			return(0);
			}
		for( i=0; i<3; i++ ) for( j=0; j<2; j++ ) {
			slow_out[i][0][j] = slow_in[0][j];
			slow_out[i][1][j] = slow_in[1][j];
			slow_out[i][2][j] = szdown[i][j];
			}
		for( i=0; i<3; i++ ) for( j=0; j<3; j++ ) for( k=0; k<2; k++ ) {
			pol_out[i][j][k] = poldown[i][j][k];
			}
		}

	/* slowness and polarizations of scattered waves in top media*/
	if( dir_in > 0 ) {
		for( i=0; i<3; i++ ) for( j=0; j<2; j++ ) {
			slow_out[i+3][0][j] = 0.0;
			slow_out[i+3][1][j] = 0.0;
			slow_out[i+3][2][j] = 0.0;
			}
		for( i=0; i<3; i++ ) for( j=0; j<3; j++ ) for( k=0; k<2; k++ ) {
			pol_out[i+3][j][k] = 0.0;
			}
		}
	else {
		if( slowness( c_above, rho_above, slow_in[0][0], slow_in[1][0],
			szup, szdown, polup, poldown ) != 1 ) {
			fprintf(stderr,"slowness (above) fails in scatter\n");
			return(0);
			}
		for( i=0; i<3; i++ ) for( j=0; j<2; j++ ) {
			slow_out[i+3][0][j] = slow_in[0][j];
			slow_out[i+3][1][j] = slow_in[1][j];
			slow_out[i+3][2][j] = szup[i][j];
			}
		for( i=0; i<3; i++ ) for( j=0; j<3; j++ ) for( k=0; k<2; k++ ) {
			pol_out[i+3][j][k] = polup[i][j][k];
			}
		}

	/* interface normal */
	normal[0]=0.0; normal[1]=0.0; normal[2]=1.0;

	/* now solve for amplitudes.  Each wave has formula (see Aki & Richards [1980], p 185)
		displacement[i] = A pol[i] exp{ -i omega*(time - slowness tangent[j]x[j]) }
		traction[i] = stress[i][j] n[j]
			    = c[i][j][p][q] strain[p][q] n[j]
			    = 0.5 * c[i][j][p][q] (u[p],[q]+u[q],[p]) n[j]
			    = c[i][j][p][q] u[p],[q] n[j]                   (since cijpq=cijqp)
			    = c[i][j][p][q] i omega slowness A pol[p] tangent[q] n[j]
	so its sufficient to require continuity of
		traction/(i omega): c[i][j][p][q] slowness A pol[p] tangent[q] n[j]
	which is to say, in the case of incidence from the bottom, the 3x3 system is:
		bottom = 0
		-(incident + S-slow-down + S-fast-down + P-down) = 0
		-S-slow-down -S-fast-down -P-down = incident
	and in the case of incidence from the top, the 6x6 system is:
		top = 0
		S-slow-up + S-fast-up + P-up + incident = 0
		S-slow-up + S-fast-up + P-up = -incident
									*/

	if( dir_in>0 ) { /*incidence from below */
		/* set up the 3x3 matrix M */
		for( wave=0; wave<3; wave++ ) for( i=0; i<3; i++ ) {
			M[i][wave][0] = 0.0;
			M[i][wave][1] = 0.0;
			for( j=0; j<3; j++ ) for( p=0; p<3; p++ ) for( q=0; q<3; q++ ) {
				cmul( pol_out[wave][p][0], pol_out[wave][p][1],
					slow_out[wave][q][0], slow_out[wave][q][1],
					&(t[0]), &(t[1]) );
				M[i][wave][0] -= c_below[i][j][p][q] * normal[j] * t[0];
				M[i][wave][1] -= c_below[i][j][p][q] * normal[j] * t[1];
				}
			}
		for( i=0; i<3; i++ ) {
			x[i][0] = 0.0;
			x[i][1] = 0.0;
			for( j=0; j<3; j++ ) for( p=0; p<3; p++ ) for( q=0; q<3; q++ ) {
				cmul(A_in[0],A_in[1], pol_in[p][0],pol_in[p][1], &(t[0]),&(t[1]));
				cmul( t[0],t[1], slow_in[q][0],slow_in[q][1], &(t[0]),&(t[1]) );
				x[i][0] += c_below[i][j][p][q] * normal[j] * t[0];
				x[i][1] += c_below[i][j][p][q] * normal[j] * t[1];
				}
			}
		}
	else if( dir_in<0 ) { /* incidence from above */
		for( wave=0; wave<3; wave++ ) for( i=0; i<3; i++ ) {
			M[i][wave][0] = 0.0;
			M[i][wave][1] = 0.0;
			for( j=0; j<3; j++ ) for( p=0; p<3; p++ ) for( q=0; q<3; q++ ) {
				cmul( pol_out[wave+3][p][0], pol_out[wave+3][p][1],
					slow_out[wave+3][q][0], slow_out[wave+3][q][1],
					&(t[0]), &(t[1]) );
				M[i][wave][0] += c_above[i][j][p][q] * normal[j] * t[0];
				M[i][wave][1] += c_above[i][j][p][q] * normal[j] * t[1];
				}
			}
		for( i=0; i<3; i++ ) {
			x[i][0] = 0.0;
			x[i][1] = 0.0;
			for( j=0; j<3; j++ ) for( p=0; p<3; p++ ) for( q=0; q<3; q++ ) {
				cmul( A_in[0],A_in[1], pol_in[p][0],pol_in[p][1], &(t[0]),&(t[1]) );
				cmul( t[0],t[1], slow_in[q][0],slow_in[q][1], &(t[0]),&(t[1]) );
				x[i][0] -= c_above[i][j][p][q] * normal[j] * t[0];
				x[i][1] -= c_above[i][j][p][q] * normal[j] * t[1];
				}
			}

		}


	/* solve 3x3 system by Gauss-Jordan reduction*/
	scaleeqn3x3(M,x);
	cgauss(M, x, 3, 3, 1.0e-6, &ierror, 1 );
	if( ierror != 0 ) {
		fprintf(stderr,"cgauss fails in fs_scatter\n");
		return(0);
		}

	if( dir_in>0 ) {
		for( wave=0; wave<3; wave++ ) for(j=0; j<2; j++) {
			A_out[wave][j] = x[wave][j];
			A_out[wave+3][j] = 0.0;
			}
		}
	else if( dir_in<0 ) {
		for( wave=0; wave<3; wave++ ) for(j=0; j<2; j++) {
			A_out[wave][j] = 0.0;
			A_out[wave+3][j] = x[wave][j];
			}
		}

	return(1);
	}


/*scales 3x3 linear system of eqns so that all rows have equal weight*/
scaleeqn3x3(M,x)
double M[3][3][2], x[3][2];
{
	int i, j;
	double f, s;
	for( i=0; i<3; i++ ) {
		f=0.0;
		for( j=0; j<3; j++ ) {
			s=sqrt( M[i][j][0]*M[i][j][0]+M[i][j][1]*M[i][j][1]);
			if( s>f ) f=s;
			}
		if( f==0.0 ) continue;
		for( j=0; j<3; j++ ) {
			M[i][j][0] /= f;
			M[i][j][1] /= f;
			}
		x[i][0] /= f;
		x[i][1] /= f;
		}
	return(1);
	}

