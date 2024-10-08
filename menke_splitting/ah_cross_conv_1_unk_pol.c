#define TRUE (1)
#define FALSE (0)
#define DTOR (3.141592654/180.0)

#define DTHETA (1)
#define DPOL (1)

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <rpc/rpc.h>
#include "SHahhead.h"

char *progname;

main(argc,argv)
int argc;
char *argv[];
{
	FILE *f1, *f2, *f3, *f4;
	XDR in1, in2, out;
	ahhed head1, head2;
	float *data1, *data2, *data3, *data4, *data5, *data6, *data7;
	double  E, e, Emin, tlag;
	int i, j, n, ii, jj, ierror, lag1, lag1p, lagmin, maxlag, itheta, first, rtflag, ipol, firstingroup;
	double ea, eb, Ea, Eb, Eamin, Ebmin, R, Rmin, theta, thetamin;
	double a0, a1, b0, b1, a0p, a1p, b0p, b1p, a0pp, a1pp, b0pp, b1pp;
	float distance, az, baz;
	double phi, cbaz, sbaz, ctheta, stheta, pol, spol, cpol, polmin, polmingroup;

	progname=argv[0];

	if( argc != 7 ) {
		fprintf(stderr,"usage: %s rt max_lag radial_ah_file tangential_ah_file output_ah_file error_surface_text_file > results_text_file\n",argv[0]);
		fprintf(stderr,"or\n");
		fprintf(stderr,"usage: %s ne max_lag north_ah_file east_ah_file output_ah_file error_surface_text_file > results_text_file\n",argv[0]);
		fprintf(stderr,"    where max_lag is in seconds, and where the flag rt or ne\n");
		fprintf(stderr,"    specifies whether radial/tangential or north/east data are supplied.\n");
		fprintf(stderr,"    Note: AH headers must contain valid receiver/event locations.\n");
		fprintf(stderr,"\n");
		fprintf(stderr,"purpose: estimate simple one-layer anisotropy parameters\n");
		fprintf(stderr,"    from a single near-vertically incident S phase of unknown initial\n");
		fprintf(stderr,"    polarization, pol, using the Menke & Levin (2002) cross-convolution Method\n");
		fprintf(stderr,"\n");
		fprintf(stderr,"Synopysis of method:\n");
		fprintf(stderr,"    When the rt flag is given and the radial horizontal component\n");
		fprintf(stderr,"    (ie. pointing back to the earthquake, A) and tangential horizontal\n");
		fprintf(stderr,"    (ie. 90 deg CCW from radial, B) component data are supplied, then we assume\n");
		fprintf(stderr,"    radial A(t) = s(t)*a(t) and tangential B(t) = s(t)*b(t) where\n");
		fprintf(stderr,"    where   s(t) = unknown source wavelet, and a(t), b(t) are the\n");
		fprintf(stderr,"    respose fuctions for a single anisotropic layer at normal incidence.\n");
		fprintf(stderr,"    The layer is parameterized by its splitting delay, t1, and the\n");
		fprintf(stderr,"    azimuth, phi, (in degrees CW from north) of its fast axis.\n");

		fprintf(stderr,"    This program determines pol, t1 and theta by a grid search by\n");
		fprintf(stderr,"    minimizing the amplitude-normalized squared L2 norm, R, of the cross-convolution,\n");
		fprintf(stderr,"            e(t) = A(t)*b(t) - B(t)*a(t).\n");
		fprintf(stderr,"    Here R = ||e||**2 / [ ||A(t)*b(t)||**2 + ||B(t)*a(t)||**2 ]\n");
		fprintf(stderr,"\n");
		fprintf(stderr,"    When the ne flag is given and the north and east horizontal\n");
		fprintf(stderr,"    components are supplied, then all calculateion are performed\n");
		fprintf(stderr,"    in the north/east coordinate system and e(t) is redefined in terms\n");
		fprintf(stderr,"    of north and east components\n");
		fprintf(stderr,"\n");
		fprintf(stderr,"output_ah_file, ah file with best-fitting\n");
		fprintf(stderr,"    A(t)*b(t)\n");
		fprintf(stderr,"    B(t)*a(t)\n");
		fprintf(stderr,"\n");
		fprintf(stderr,"error_surface_text_file, the backazimuth followed by a table of\n");
		fprintf(stderr,"    phi, azimuth of fast axis in deg east of north\n");
		fprintf(stderr,"    t1, anisotropic delay time in seconds\n");
		fprintf(stderr,"    pol, best-fitting polarization angle, in degrees east of north\n");
		fprintf(stderr,"        for given values of phi and t1\n");
		fprintf(stderr,"    R, value of amplitude-normalized squared L2 norm\n");
		fprintf(stderr,"\n");
		fprintf(stderr,"results_text_file, a table of:\n");
		fprintf(stderr,"    backazimuth of earthquake in deg east of north\n");
		fprintf(stderr,"    phi, best-fitting azimuth of fast axis in deg east of north\n");
		fprintf(stderr,"    t1, best-fitting anisotropic delay time in seconds\n");
		fprintf(stderr,"    pol, best-fitting polarization angle, in degrees east of north\n");
		fprintf(stderr,"    R, minimum value of amplitude-normalized squared L2 norm\n");
		fprintf(stderr,"\n");
		exit(-1);
		}

	if( !strcmp(argv[1],"rt") ) {
		rtflag=TRUE;
		}
	else if( !strcmp(argv[1],"ne") ) {
		rtflag=FALSE;
		}
	else {
		fprintf(stderr,"error: %s: flag must be rt or ne\n", argv[0] );
		exit(-1);
		}

	if( sscanf(argv[2],"%le",&tlag) != 1 ) {
		fprintf(stderr,"error: %s: can read max_lag from command line\n", argv[0] );
		exit(-1);
		}

	if( (f1=fopen(argv[3],"r"))==NULL ) {
		fprintf(stderr,"error: %s: cant open radial (or north) file <%s> for read\n", argv[0], argv[3]);
		exit(-1);
		}

	if( (f2=fopen(argv[4],"r"))==NULL ) {
		fprintf(stderr,"error: %s: cant open tangential (or east) file <%s> for read\n", argv[0], argv[4]);
		exit(-1);
		}

	if( (f3=fopen(argv[5],"w"))==NULL ) {
		fprintf(stderr,"error: %s: cant open output_ah_file <%s> for write\n", argv[0], argv[5]);
		exit(-1);
		}

	if( (f4=fopen(argv[6],"w"))==NULL ) {
		fprintf(stderr,"error: %s: cant open error_surface_text_file <%s> for write\n", argv[0], argv[6]);
		exit(-1);
		}

	xdrstdio_create(&in1,  f1,  XDR_DECODE);
	xdrstdio_create(&in2,  f2,  XDR_DECODE);
	xdrstdio_create(&out,  f3,  XDR_ENCODE);

	fprintf(stdout,"baz\tphi\tt1\tpol\tR\n");

	j=0;
	while( xdr_getrecord2(&head1,&data1,&in1) == 0 ) {

		if( xdr_getrecord2(&head2,&data2,&in2) != 0 ) {
			fprintf(stderr,"error: %s: cant read record in ah_file_2\n", argv[0] );
			exit(-1);
			}

		if( head1.record.delta != head2.record.delta ) {
			fprintf(stderr,"error: %s: deltas dont match\n", argv[0]);
			exit(-1);
			}

		if( head1.record.ndata != head2.record.ndata ) {
			fprintf(stderr,"error: %s: ndatas dont match\n", argv[0]);
			exit(-1);
			}

     		menke_distaz( head1.event.lat, head1.event.lon,
			      head1.station.slat, head1.station.slon,
			      &distance, &az, &baz );
		cbaz = cos( DTOR*baz );
		sbaz = sin( DTOR*baz );
		fprintf(f4,"%f\n", baz);

		n=head1.record.ndata;
		maxlag = (int) (0.5+tlag/head1.record.delta);

		if( (data3=(float*)malloc(n*sizeof(float))) == NULL ) {
			fprintf(stderr,"error: %s: out of memory!\n", argv[0]);
			exit(-1);
			}
		if( (data4=(float*)malloc(n*sizeof(float))) == NULL ) {
			fprintf(stderr,"error: %s: out of memory!\n", argv[0]);
			exit(-1);
			}
		if( (data5=(float*)malloc(n*sizeof(float))) == NULL ) {
			fprintf(stderr,"error: %s: out of memory!\n", argv[0]);
			exit(-1);
			}
		if( (data6=(float*)malloc(n*sizeof(float))) == NULL ) {
			fprintf(stderr,"error: %s: out of memory!\n", argv[0]);
			exit(-1);
			}
		if( (data7=(float*)malloc(n*sizeof(float))) == NULL ) {
			fprintf(stderr,"error: %s: out of memory!\n", argv[0]);
			exit(-1);
			}

		first=TRUE;

		/*theta: CCW from pol to fast */
		for( itheta=0; itheta<180; itheta+=DTHETA ) { 

			theta = (double) itheta;
			ctheta=cos(DTOR*theta); stheta=sin(DTOR*theta);
			a0pp = ctheta*ctheta;
			a1pp = stheta*stheta;
			b0pp = stheta*ctheta;
			b1pp = (-b0pp);

		for( lag1p=1; lag1p<=maxlag; lag1p++ ) {

			lag1=(lag1p);
			firstingroup=TRUE;

		/* pol: CW from N to polarization direction */
		for( ipol=0; ipol<180; ipol+=DPOL ) {

			pol = (double) ipol;
			a0 = a0pp; a1=a1pp; b0=b0pp; b1=b1pp;

			/* response in coordinate system aligned with pol */

			/* rotate to coordinate system aligned with radial */
			/* rotation angle is baz-pol                       */
			cpol=cos(DTOR*(baz-pol)); spol=sin(DTOR*(baz-pol));
			a0p = cpol*a0 - spol*b0;
			b0p = spol*a0 + cpol*b0;
			a1p = cpol*a1 - spol*b1;
			b1p = spol*a1 + cpol*b1;

			if( rtflag ) {
				a0=a0p; a1=a1p; b0=b0p; b1=b1p;
				}
			else {
				a0 = cbaz*a0p + sbaz*b0p;
				a1 = cbaz*a1p + sbaz*b1p;
				b0 = sbaz*a0p - cbaz*b0p;
				b1 = sbaz*a1p - cbaz*b1p;
				}

			E = 0.0; Ea=0.0; Eb=0.0;
			for( ii=0; ii<n; ii++ ) {
				jj = ii-lag1;
				if( jj<0 ) jj+=n; /* circular */
				e =  b0 * data1[ii];
				e += a0 * (-data2[ii]);
				ea = b0 * data1[ii];
				eb = a0 * data2[ii];
				e  += b1 * data1[jj] + a1 * (-data2[jj]);
				ea += b1 * data1[jj];
				eb += a1 * data2[jj];
				data3[ii]=(float)ea;
				data4[ii]=(float)eb;
				data5[ii]=(float)e;
				E += e*e; Ea += ea*ea; Eb += eb*eb;
				}

			R = E/(Ea+Eb);

			if( first ) {
				Emin=E; Eamin=Ea; Ebmin=Eb; Rmin=R;
				lagmin=lag1; thetamin=theta; polmin=pol; first=FALSE;
				for( ii=0; ii<n; ii++ ) {
					data6[ii]=data3[ii];
					data7[ii]=data4[ii];
					}
				}
			else if( R<Rmin ) {
				Emin=E; Eamin=Ea; Ebmin=Eb; Rmin=R;
				lagmin=lag1; thetamin=theta; polmin=pol;
				for( ii=0; ii<n; ii++ ) {
					data6[ii]=data3[ii];
					data7[ii]=data4[ii];
					}
				}

			if( firstingroup ) {
				polmingroup=pol; firstingroup=FALSE;
				}
			else if( R<Rmin ) {
				polmingroup=pol;
				}

			}

			phi = polmingroup-theta; if( phi<0.0 ) phi+=360.0; else if( phi>=360.0 ) phi-=360.0;
			if( phi>180.0 ) phi-=180.0;

			fprintf(f4,"%f\t%f\t%f\t%f\n", phi, head1.record.delta*(double)lag1,
					polmingroup, R );

			}}


		phi = polmin-thetamin; if( phi<0.0 ) phi+=360.0; else if( phi>=360.0 ) phi-=360.0;
		if( phi>180.0 ) phi-=180.0;
		
		fprintf(stdout,"%f\t%f\t%f\t%f\t%f\n", baz, phi, head1.record.delta*(double)lagmin, polmin, Rmin);
		fflush(stdout);

		maxamp(&head1,data6);
                if( xdr_putrecord(&head1,data6,&out) != 0) {
                        fprintf(stderr,"Error writing record %d %s\n",j,progname);
                        exit(-1);
                        }
		maxamp(&head1,data7);
                if( xdr_putrecord(&head1,data7,&out) != 0) {
                        fprintf(stderr,"Error writing record %d %s\n",j,progname);
                        exit(-1);
                        }
 
		free(data1); free(data2); free(data3); free(data4);
		free(data5); free(data6); free(data7);
		j++;
		}
	}


