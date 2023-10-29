#define TRUE (1)
#define FALSE (0)
#define DTOR (3.141592654/180.0)
#define MAXTRACES (100)

#define DTHETA (1)
#define DPHI   (1)

#include <stdio.h>
#include <math.h>
//#include <sunmath.h>
#include <malloc.h>
#include <rpc/rpc.h>
#include "SHahhead.h"

double lagdot( float *, int, float *, int, int );

char *progname;

main(int argc, char *argv[])
{
	FILE *f1, *f2, *f3, *f4;
	XDR in1, in2, out;
	ahhed head1[MAXTRACES], head2[MAXTRACES];
	float *data1[MAXTRACES], *data2[MAXTRACES], *data3, *data4;
	int i, j, n, ii, jj, kk, maxtau, itheta, iphi, itau1, itau2, first, firstingroup, ntraces;
	int tau1, tau2, tau1min, tau2min, maxn, tau1min0, tau2min0;
	double theta, phi, ct, st, cp, sp, ap[4], bp[4], a[4], b[4];
	double  E, Emin, taumax, thetamin, phimin;
	double Ea, Eb, R, delta, Rall, Rallmin, Rallmin0;
	float distance, az, baz[MAXTRACES];
	int clag[4], rtflag;
	double cbaz[MAXTRACES], sbaz[MAXTRACES];
	double tau1mingroup, tau2mingroup, Rallmingroup, faz1, faz2;

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
		fprintf(stderr,"    from a group of radially-polarized (e.g. SKS) phases\n");
		fprintf(stderr,"    using the Menke & Levin (2002) cross-convolution Method\n");
		fprintf(stderr,"\n");
		fprintf(stderr,"Synopysis of method:\n");
		fprintf(stderr,"    When the rt flag is given and the radial horizontal component\n");
		fprintf(stderr,"    (ie. pointing back to the earthquake, A) and tangential horizontal\n");
		fprintf(stderr,"    (ie. 90 deg CCW from radial, B) component data are supplied, then we assume\n");
		fprintf(stderr,"    radial A(t) = s(t)*a(t) and tangential B(t) = s(t)*b(t) where\n");
		fprintf(stderr,"    where s(t) = unknown source wavelet, and a(t), b(t) are the\n");
		fprintf(stderr,"    respose fuctions for two anisotropic layers at normal incidence.\n");
		fprintf(stderr,"    This program determines the azimuths, phi1 and phi2, and delay\n");
		fprintf(stderr,"    times, t1 and t2, of the bottom and top layers, respecitively,\n");
		fprintf(stderr,"    by a grid search that minimizes the mean (over the group of N\n");
		fprintf(stderr,"    timeseries) of the amplitude-normalized squared L2 norm, R,\n");
		fprintf(stderr,"    of the cross-convolution,\n");
		fprintf(stderr,"            e(t) = A(t)*b(t) - B(t)*a(t).\n");
		fprintf(stderr,"    Here R = (1/N) (sum) ||e||**2 / [ ||A(t)*b(t)||**2 + ||B(t)*a(t)||**2 ]\n");
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
		fprintf(stderr,"error_surface_text_file, a line with a zero on it followed by a table of:\n");
		fprintf(stderr,"    phi1, azimuth of fast axis of the bottom layer in deg east of north\n");
		fprintf(stderr,"    phi2, azimuth of fast axis of the top layer in deg east of north\n");
		fprintf(stderr,"    t1, best-fitting anisotropic delay time of the bottom layer in seconds\n");
		fprintf(stderr,"        with respect to this particular (phi1, phi2).\n");
		fprintf(stderr,"    t1, best-fitting anisotropic delay time of the top layer in seconds\n");
		fprintf(stderr,"        with respect to this particular (phi1, phi2).\n");
		fprintf(stderr,"    R, corresponding minimum value of amplitude-normalized squared L2 norm\n");
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
		fprintf(stderr,"error: %s: cant open ah_file_1 <%s> for read\n", argv[0], argv[3]);
		exit(-1);
		}

	if( (f2=fopen(argv[4],"r"))==NULL ) {
		fprintf(stderr,"error: %s: cant open ah_file_2 <%s> for read\n", argv[0], argv[4]);
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

	j=0;
	while( xdr_getrecord2(&(head1[j]),&(data1[j]),&in1) == 0 ) {

		if( j==0 ) {
			delta=head1[j].record.delta;
			maxtau = (int) (0.5+taumax/delta);
			maxn=head1[j].record.ndata;
			}
		else if( head1[j].record.ndata > maxn ) {
			maxn=head1[j].record.ndata;
			}

		if( xdr_getrecord2(&(head2[j]),&(data2[j]), &in2) != 0 ) {
			fprintf(stderr,"error: %s: cant read record %d in ah_file_2\n", argv[0], j );
			exit(-1);
			}

		if( (head1[j].record.delta != delta) || (head2[j].record.delta != delta) ) {
			fprintf(stderr,"error: %s: deltas dont match (timeseries %d)\n", argv[0], j);
			exit(-1);
			}

		if( head1[j].record.ndata != head2[j].record.ndata ) {
			fprintf(stderr,"error: %s: ndatas dont match (timeseries %d)\n", argv[0], j);
			exit(-1);
			}

		menke_distaz( head1[j].event.lat, head1[j].event.lon,
			      head1[j].station.slat, head1[j].station.slon,
			      &distance, &az, &(baz[j]) );

		baz[j] = DTOR * baz[j];
		cbaz[j]=cos(baz[j]); sbaz[j]=sin(baz[j]);
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

	fprintf(f4,"%f\n", 0.0 );

	for( itheta=0; itheta<180; itheta+=DTHETA ) {
	theta = DTOR * (double) itheta;
	for( iphi=0; iphi<180; iphi+=DPHI ) {
	phi =  DTOR * (double) iphi; 
	firstingroup=TRUE;
		for( itau1=0; itau1<=maxtau; itau1++ ) {
		tau1=itau1; 
		for( itau2=0; itau2<=maxtau; itau2++ ) {
		tau2=itau2;
		clag[0]=0; clag[1]=tau2; clag[2]=tau1; clag[3]=tau1+tau2;

		Rall = 0.0;
		for( j=0; j<ntraces; j++ ) {
		
			n=head1[j].record.ndata;

			ct = cos(baz[j]-theta); st = sin(baz[j]-theta);
			cp = cos(baz[j]-phi);   sp = sin(baz[j]-phi);
	
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
					a[jj] = cbaz[j]*ap[jj] + sbaz[j]*bp[jj];
					b[jj] = sbaz[j]*ap[jj] - cbaz[j]*bp[jj];
					}
				}

			for( ii=0; ii<n; ii++ ) {
				data3[ii]=0.0; data4[ii]=0.0;
				for( jj=0; jj<4; jj++ ) {
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

			}

		/* this is the global minimum */
		if( first ) {
			Rallmin=Rall;
			tau1min=tau1; tau2min=tau2;
			thetamin=theta; phimin=phi; first=FALSE;
			}
		else if( Rall<Rallmin ) {
			Rallmin=Rall;
			tau1min=tau1; tau2min=tau2;
			thetamin=theta; phimin=phi; first=FALSE;
			}

		/* this is the minimum over the current (theta, phi) */
		if( firstingroup ) {
			Rallmingroup=Rall;
			tau1mingroup=tau1; tau2mingroup=tau2;
			firstingroup=FALSE;
			}
		else if( Rall<Rallmingroup ) {
			Rallmingroup=Rall;
			tau1mingroup=tau1; tau2mingroup=tau2;
			firstingroup=FALSE;
			}

		}}
		fprintf(f4,"%f\t%f\t%f\t%f\t%f\n", theta/DTOR, delta*tau1mingroup, 
						   phi/DTOR,   delta*tau2mingroup,
						   Rallmingroup/(double)ntraces );
		fflush(f4);
	}}

	fprintf(stdout,"t0\tphi1\tt2\tphi2\tR\n");
	fprintf(stdout,"%f\t%f\t%f\t%f\t%f\n", delta*(double)tau1min, thetamin/DTOR,
			                       delta*(double)tau2min, phimin/DTOR,
					       Rallmin/(double)ntraces );
	fflush(stdout);

	theta=thetamin; phi=phimin;
	tau1=tau1min; tau2=tau2min;
	clag[0]=0; clag[1]=tau2; clag[2]=tau1; clag[3]=tau1+tau2;

	for( j=0; j<ntraces; j++ ) {
		
		n=head1[j].record.ndata;

		ct = cos(baz[j]-theta); st = sin(baz[j]-theta);
		cp = cos(baz[j]-phi);   sp = sin(baz[j]-phi);
	
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
				a[jj] = cbaz[j]*ap[jj] + sbaz[j]*bp[jj];
				b[jj] = sbaz[j]*ap[jj] - cbaz[j]*bp[jj];
				}
			}

		for( ii=0; ii<n; ii++ ) {
			data3[ii]=0.0; data4[ii]=0.0;
			for( jj=0; jj<4; jj++ ) {
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

		}

	exit(0);
	}



