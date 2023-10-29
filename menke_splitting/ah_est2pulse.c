#define TRUE (1)
#define FALSE (0)

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
	FILE *f1, *f2, *f3;
	XDR in1, in2, out;
	ahhed head1, head2;
	float *data1, *data2, *data3, *data4;
	double GTG[3][3], GTd[3], E, e, Emin, soln[3], tlag;
	int i, j, n, ii, jj, ierror, lag1, lag1p, lagmin, maxlag;
	double ea, eb, Ea, Eb, Eamin, Ebmin, R, Rmin;

	progname=argv[0];

	if( argc != 5 ) {
		fprintf(stderr,"usage: %s max_lag A(t)_ah_file B(t)_ah_file e(t)_ah_file > results.txt\n", argv[0]);
		fprintf(stderr,"       with max_lag in seconds.\n");
		fprintf(stderr,"\n");
		fprintf(stderr,"purpose: estimate the 2-pulse operator that relates two\n");
		fprintf(stderr,"    timeseries A(t) and B(t)\n");
		fprintf(stderr,"\n");
		fprintf(stderr,"    Suppose A(t) = s(t)*a(t) and B(t) = s(t)*b(t)\n");
		fprintf(stderr,"    where   a(t) = (sum i=0,1) ai delta(t-ti) and\n");
		fprintf(stderr,"                b(t) = (sum i=0,1) bi delta(t-ti)\n");
		fprintf(stderr,"    and where t0=0 and s(t) is unknown.  This program determines\n");
		fprintf(stderr,"    the amplitudes ai and bi (up to a multiplicative constant,\n");
		fprintf(stderr,"    e.g. b0=1) by the Menke & Levin cross-conolution method.  The\n");
		fprintf(stderr,"    idea is to minimize the norm, R, of the prediction\n");
		fprintf(stderr,"    error, e(t) = A(t)*b(t) - B(t)*a(t), where\n");
		fprintf(stderr,"    R = ||e||**2 / [ ||A(t)*b(t)||**2 + ||B(t)*a(t)||**2 ]\n");
		fprintf(stderr,"    The best-fitting ti, ai, bi, R are written to stdout\n");
		fprintf(stderr,"    and A(t)*b(t) and B(t)*a(t) written to the e(t)_ah_file.\n");
		fprintf(stderr,"\n");
		fprintf(stderr,"    When using this routine to to solve for a one-layer anisotropy operator\n");
		fprintf(stderr,"    I recommend that A(t) represent the tangential-horizontal component\n");
		fprintf(stderr,"    and B(t) represent the radial-horizontal component.\n");
 		fprintf(stderr,"\n");
		fprintf(stderr,"    When using this routine to to solve for receiver function operator\n");
		fprintf(stderr,"    I recommend that A(t) represent the radial-horizontal component\n");
		fprintf(stderr,"    and B(t) represent the vertical component.\n");

		exit(-1);
		}

	fprintf(stdout,"t0\tt1\ta0\ta1\tb0\tb1\tR\n");

	if( sscanf(argv[1],"%le",&tlag) != 1 ) {
		fprintf(stderr,"error: %s: can read max_lag from command line\n", argv[0] );
		exit(-1);
		}

	if( (f1=fopen(argv[2],"r"))==NULL ) {
		fprintf(stderr,"error: %s: cant open A(t)_ah_file <%s> for read\n", argv[0],argv[2]);
		exit(-1);
		}

	if( (f2=fopen(argv[3],"r"))==NULL ) {
		fprintf(stderr,"error: %s: cant open B(t)_ah_file <%s> for read\n", argv[0], argv[3]);
		exit(-1);
		}

	if( (f3=fopen(argv[4],"w"))==NULL ) {
		fprintf(stderr,"error: %s: cant open e(t)_file <%s> for write\n", argv[4]);
		exit(-1);
		}

	xdrstdio_create(&in1,  f1,  XDR_DECODE);
	xdrstdio_create(&in2,  f2,  XDR_DECODE);
	xdrstdio_create(&out, f3, XDR_ENCODE);

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

		for( lag1p=1; lag1p<=maxlag; lag1p++ ) {

			lag1=(lag1p);

			GTG[0][0] =   lagdot( data1, lag1, data1, lag1, n );
			GTG[0][1] = (-lagdot( data1, lag1, data2, 0,    n ));
			GTG[0][2] = (-lagdot( data1, lag1, data2, lag1, n ));
			GTG[1][1] =   lagdot( data2, 0,    data2, 0,    n );
			GTG[1][2] =   lagdot( data2, 0,    data2, lag1, n );
			GTG[2][2] =   lagdot( data2, lag1, data2, lag1, n );

			GTG[1][0] = GTG[0][1];
			GTG[2][0] = GTG[0][2];
			GTG[2][1] = GTG[1][2];

			GTd[0] = (-lagdot( data1, lag1, data1, 0, n ));
			GTd[1] =   lagdot( data2, 0,    data1, 0, n );
			GTd[2] =   lagdot( data2, lag1, data1, 0, n );

			gauss(GTG,GTd,3,3,1.0e-6,&ierror,TRUE);

			E = 0.0; Ea=0.0; Eb=0.0;
			for( ii=0; ii<n; ii++ ) {
				jj = ii-lag1;
				e = data1[ii];
				e += GTd[1] * (-data2[ii]);
				ea = data1[ii];
				eb = GTd[1] * data2[ii];
				if( (jj>=0) && (jj<n) ) {
					e += GTd[0] * data1[jj] + GTd[2] * (-data2[jj]);
					ea += GTd[0] * data1[jj];
					eb += GTd[2] * data2[jj];
					}
				data3[ii]=(float)ea;
				data4[ii]=(float)eb;
				E += e*e; Ea += ea*ea; Eb += eb*eb;
				}

			R = E/(Ea+Eb);

			if( lag1p==1 ) {
				soln[0]=GTd[0]; soln[1]=GTd[1]; soln[2]=GTd[2];
				Emin=E; Eamin=Ea; Ebmin=Eb; Rmin=R;
				lagmin=lag1;
				}
			else if( R<Rmin ) {
				soln[0]=GTd[0]; soln[1]=GTd[1]; soln[2]=GTd[2];
				Emin=E; Eamin=Ea; Ebmin=Eb; Rmin=R;
				lagmin=lag1;
				}
			}

		lag1=lagmin;
		E = 0.0; Ea=0.0; Eb=0.0;
		for( ii=0; ii<n; ii++ ) {
			jj = ii-lag1;
			e = data1[ii];
			e += soln[1] * (-data2[ii]);
			ea = data1[ii];
			eb = soln[1] * data2[ii];
			if( (jj>=0) && (jj<n) ) {
				e += soln[0] * data1[jj] + soln[2] * (-data2[jj]);
				ea += soln[0] * data1[jj];
				eb += soln[2] * data2[jj];
				}
			data3[ii]=(float)ea;
			data4[ii]=(float)eb;
			E += e*e; Ea += ea*ea; Eb += eb*eb;
			}
		R = E/(Ea+Eb);

		fprintf(stdout,"%f\t%f\t%f\t%f\t%f\t%f\t%e\n", 0.0,
			head1.record.delta*(double)lagmin, soln[1], soln[2], 1.0, soln[0], Rmin);

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
	
		free(data1);free(data2);free(data3);free(data4);

		j++;
		}
	}

double lagdot( float *d1, int l1, float* d2, int l2, int n) {
	int i, j, k;
	double sum;

	sum = 0.0;
	for( i=0; i<n; i++ ) {
		j=i-l1; if( (j<0) || (j>=n) ) continue;
		k=i-l2; if( (k<0) || (k>=n) ) continue;
		sum+=d1[j]*d2[k];
		}

	return(sum);
	}


