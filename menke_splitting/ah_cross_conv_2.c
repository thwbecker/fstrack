#define TRUE (1)
#define FALSE (0)
#define DTOR (3.141592654/180.0)

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <rpc/rpc.h>
#include "SHahhead.h"

double lagdot( float *, int, float *, int, int );

char *progname;

main(argc,argv)
int argc;
char *argv[];
{
	FILE *f1, *f2, *f3, *f4;
	XDR in1, in2, out;
	ahhed head1, head2;
	float *data1, *data2, *data3, *data4, *data5, *data6;
	int i, j, n, ii, jj, kk, maxtau, itheta, iphi, itau1, itau2, first, firstingroup, rtflag;
	int tau1, tau2, tau1min, tau2min, tau1mingroup, tau2mingroup;
	double theta, phi, ct, st, cp, sp, a[4], b[4], ap[4], bp[4];
	double  E, Emin, Emingroup, taumax, thetamin, phimin;
	double Ea, Eb, Eamin, Ebmin, Eamingroup, Ebmingroup, R, Rmin, Rmingroup;
	int clag[4];
	float distance, az, baz;
	double cbaz, sbaz, faz1, faz2;

	progname=argv[0];

	if( argc != 7 ) {
		fprintf(stderr,"usage: %s rt max_lag radial_ah_file tangential_ah_file output_ah_file error_surface_text_file > results_text_file\n",argv[0]);
		fprintf(stderr,"or\n");
		fprintf(stderr,"usage: %s ne max_lag north_ah_file east_ah_file output_ah_file error_surface_text_file > results_text_file\n",argv[0]);
		fprintf(stderr,"    where max_lag is in seconds, and where the flag rt or ne\n");
		fprintf(stderr,"    specifies whether radial/tangential or north/east data are supplied.\n");
		fprintf(stderr,"    Note: AH headers must contain valid receiver/event locations.\n");
		fprintf(stderr,"\n");
		fprintf(stderr,"purpose: estimate simple two-layer anisotropy parameters\n");
		fprintf(stderr,"    from a single radially-polarized (e.g. SKS) phase\n");
		fprintf(stderr,"    using the Menke & Levin (2002) cross-convolution Method\n");
		fprintf(stderr,"\n");
		fprintf(stderr,"Synopysis of method:\n");
		fprintf(stderr,"    When the rt flag is given and the radial horizontal component\n");
		fprintf(stderr,"    (ie. pointing back to the earthquake, A) and tangential horizontal\n");
		fprintf(stderr,"    (ie. 90 deg CCW from radial, B) component data are supplied, then we assume\n");
		fprintf(stderr,"    radial A(t) = s(t)*a(t) and tangential B(t) = s(t)*b(t) where\n");
		fprintf(stderr,"    where   s(t) = unknown source wavelet, and a(t), b(t) are the\n");
		fprintf(stderr,"    respose fuctions for a two anisotropic layers at normal incidence.\n");
		fprintf(stderr,"    Each layer is parameterized by its delay time, t, and the azimuth,\n");
		fprintf(stderr,"    phi, of its fast direction.\n");
		fprintf(stderr,"\n");
		fprintf(stderr,"    This program determines t1 and theta by a grid search by\n");
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
		fprintf(stderr,"    phi1, azimuth of fast axis of the bottom layer in deg east of north\n");
		fprintf(stderr,"    phi2, azimuth of fast axis of the top layer in deg east of north\n");
		fprintf(stderr,"    t1, for the given combination of (phi1, phi2), the best-fitting\n");
		fprintf(stderr,"        anisotropic delay time of the bottom layer in seconds\n");
		fprintf(stderr,"    t2, for the given combination of (phi1, phi2), the best-fitting\n");
		fprintf(stderr,"        anisotropic delay time of the top layer in seconds\n");
		fprintf(stderr,"    R, minimum value of amplitude-normalized squared L2 norm for this (phi1, phi2)\n");
		fprintf(stderr,"\n");
		fprintf(stderr,"results_text_file, a table of:\n");
		fprintf(stderr,"    backazimuth of earthquake in deg east of north\n");
		fprintf(stderr,"    phi1, best-fitting azimuth of fast axis of the bottom layer in deg east of north\n");
		fprintf(stderr,"    t1, best-fitting anisotropic delay time of the bottom layer in seconds\n");
		fprintf(stderr,"    phi2, best-fitting azimuth of fast axis of the top layer in deg east of north\n");
		fprintf(stderr,"    t2, best-fitting anisotropic delay time of the top layer in seconds\n");
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

	if( sscanf(argv[2],"%le",&taumax) != 1 ) {
		fprintf(stderr,"error: %s: can read max_tau from command line\n", argv[0] );
		exit(-1);
		}

	if( (f1=fopen(argv[3],"r"))==NULL ) {
		fprintf(stderr,"error: %s: cant open cant open radial (or north) file <%s> for read\n", argv[0], argv[3]);
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

	fprintf(stdout,"baz\tphi1\ttau1\tphi2\ttau2\tR\n");

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
		maxtau = (int) (0.5+taumax/head1.record.delta);

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

		first=TRUE;

		for( itheta=0; itheta<180; itheta++ ) {
			theta = DTOR * (double) itheta;
			ct = cos(theta); st = sin(theta);
		for( iphi=0; iphi<180; iphi++ ) {
			phi =  DTOR * (double) iphi;
			cp = cos(phi);   sp = sin(phi);
			ap[0] = ct*ct*cp*cp + ct*st*cp*sp;
			ap[1] = ct*ct*sp*sp - ct*st*cp*sp;
			ap[2] = st*st*cp*cp - ct*st*cp*sp;
			ap[3] = st*st*sp*sp + ct*st*cp*sp;
			bp[0] = ct*ct*cp*sp + ct*st*sp*sp;
			bp[1] = (-ct*ct*cp*sp) + ct*st*cp*cp;
			bp[2] = st*st*cp*sp - ct*st*sp*sp;
			bp[3] = (-st*st*cp*sp) - ct*st*cp*cp;
			if( rtflag ) {
				for( jj=0; jj<4; jj++ ) {
					a[jj]=ap[jj];
					b[jj]=bp[jj];
					}
				}
			else {
				for( jj=0; jj<4; jj++ ) {
					a[jj] = cbaz*ap[jj] + sbaz*bp[jj];
					b[jj] = sbaz*ap[jj] - cbaz*bp[jj];
					}
				}
		firstingroup=TRUE;
		for( itau1=0; itau1<=maxtau; itau1++ ) { 
			tau1=itau1; 
		for( itau2=0; itau2<=maxtau; itau2++ ) { 
			tau2=itau2;
			clag[0]=0; clag[1]=tau2; clag[2]=tau1; clag[3]=tau1+tau2;

			for( ii=0; ii<n; ii++ ) {
				data3[ii]=0.0; data4[ii]=0.0;
				for( jj=0; jj<4; jj++ ) {
					kk = ii-clag[jj];
					if(kk<0) kk+=n; /*circular*/
					data3[ii]+=b[jj]*data1[kk];
					data4[ii]+=a[jj]*data2[kk];
					}
				}

			E = 0.0; Ea = 0.0; Eb = 0.0;
			for( ii=0; ii<n; ii++ ) {
				Ea += data3[ii]*data3[ii];
				Eb += data4[ii]*data4[ii];
				E  += (data3[ii]-data4[ii]) * (data3[ii]-data4[ii]);
				}

			R = E/(Ea+Eb);

			/* this is the global minimum */
			if( first ) {
				for( ii=0; ii<n; ii++ ) { data5[ii]=data3[ii]; data6[ii]=data4[ii]; }
				Emin=E; Eamin=Ea; Ebmin=Eb; Rmin=R;
				tau1min=tau1; tau2min=tau2;
				thetamin=theta; phimin=phi; first=FALSE;
				}
			else if( R<Rmin ) {
				for( ii=0; ii<n; ii++ ) { data5[ii]=data3[ii]; data6[ii]=data4[ii]; }
				Emin=E; Eamin=Ea; Ebmin=Eb; Rmin=R;
				tau1min=tau1; tau2min=tau2;
				thetamin=theta; phimin=phi; first=FALSE;
				}

			/* this is the minimum with respect to the current (theta, phi)*/
			if( firstingroup ) {
				Emingroup=E; Eamingroup=Ea; Ebmingroup=Eb; Rmingroup=R;
				tau1mingroup=tau1; tau2mingroup=tau2;
				firstingroup=FALSE;
				}
			else if( R<Rmingroup ) {
				Emingroup=E; Eamingroup=Ea; Ebmingroup=Eb; Rmingroup=R;
				tau1mingroup=tau1; tau2mingroup=tau2;
				firstingroup=FALSE;
				}
			}}
			faz1 = baz-(theta/DTOR); if( faz1<0.0 ) faz1+=360.0; else if( faz1>=360.0 ) faz1-=360.0;
			if( faz1>180.0 ) faz1-=180.0;
			faz2 = baz-(phi/DTOR); if( faz2<0.0 ) faz2+=360.0; else if( faz2>=360.0 ) faz2-=360.0;
			if( faz2>180.0 ) faz2-=180.0;
			fprintf(f4,"%f\t%f\t%f\t%f\t%f\t%f\n", baz, faz1, faz2,
				head1.record.delta*tau1mingroup, head1.record.delta*tau2mingroup, Rmingroup );
			}}

		faz1 = baz-(thetamin/DTOR); if( faz1<0.0 ) faz1+=360.0; else if( faz1>=360.0 ) faz1-=360.0;
		if( faz1>180.0 ) faz1-=180.0;
		faz2 = baz-(phimin/DTOR); if( faz2<0.0 ) faz2+=360.0; else if( faz2>=360.0 ) faz2-=360.0;
		if( faz2>180.0 ) faz2-=180.0;
		fprintf(stdout,"%f\t%f\t%f\t%f\t%f\n", faz1, head1.record.delta*tau1min,
                                                   faz2, head1.record.delta*tau2min, Rmin );
		fflush(stdout);
		maxamp(&head1,data5);
                if( xdr_putrecord(&head1,data5,&out) != 0) {
                        fprintf(stderr,"Error writing record %d %s\n",j,progname);
                        exit(-1);
                        }
		maxamp(&head1,data6);
                if( xdr_putrecord(&head1,data6,&out) != 0) {
                        fprintf(stderr,"Error writing record %d %s\n",j,progname);
                        exit(-1);
                        }
 
		free(data1);free(data2);free(data3);free(data4);free(data5);free(data6);
		j++;
		}
	}



