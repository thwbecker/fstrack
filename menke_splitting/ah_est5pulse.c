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
	double GTG[9][9], GTd[9], E, e, ea, eb, Ea, Eb, Emin, Eamin, Ebmin, R, Rmin, GTGmax;
	double soln[9], tlag, sn[9], temp, damping;
	float *c[9];
	int i, j, n, ii, jj, kk, ierror, maxlag, clag[9];
	int lag1, lag1p, lag2, lag2p, lag3, lag3p, lag4p, lag4;
	int lag1min, lag2min, lag3min, lag4min, first;

	progname=argv[0];
	damping=1.0e-5;

	if( argc != 5 ) {
		fprintf(stderr,"usage: %s max_lag A(t)_ah_file B(t)_ah_file e(t)_ah_file > results.txt\n", argv[0]);
		fprintf(stderr,"       with max_lag in seconds.\n");
		fprintf(stderr,"\n");
		fprintf(stderr,"purpose: estimate the 5-pulse operator that relates two\n");
		fprintf(stderr,"    timeseries A(t) and B(t)\n");
		fprintf(stderr,"\n");
		fprintf(stderr,"    Suppose A(t) = s(t)*a(t) and B(t) = s(t)*b(t)\n");
		fprintf(stderr,"    where   a(t) = (sum i=0,4) ai delta(t-ti) and\n");
		fprintf(stderr,"                b(t) = (sum i=0,4) bi delta(t-ti)\n");
		fprintf(stderr,"    and where t0=0 and s(t) is unknown.  This program determines\n");
		fprintf(stderr,"    the amplitudes ai and bi (up to a multiplicative constant,\n");
		fprintf(stderr,"    e.g. b0=1) by the Menke & Levin cross-conolution method.  The\n");
		fprintf(stderr,"    idea is to minimize the norm, R, of the prediction\n");
		fprintf(stderr,"    error, e(t) = A(t)*b(t) - B(t)*a(t), where\n");
		fprintf(stderr,"    R = ||e||**2 / [ ||A(t)*b(t)||**2 + ||B(t)*a(t)||**2 ]\n");
		fprintf(stderr,"    The best-fitting ti, ai, bi, R are written to stdout\n");
		fprintf(stderr,"    and A(t)*b(t) and B(t)*a(t) written to the e(t)_ah_file.\n");
		fprintf(stderr,"\n");
		fprintf(stderr,"    When using this routine to to solve for receiver function operator\n");
		fprintf(stderr,"    I recommend that A(t) represent the radial-horizontal component\n");
		fprintf(stderr,"    and B(t) represent the vertical component.\n");
		exit(-1);
		}

	if( sscanf(argv[1],"%le",&tlag) != 1 ) {
		fprintf(stderr,"error: %s: can read max_lag from command line\n", argv[0] );
		exit(-1);
		}

	if( (f1=fopen(argv[2],"r"))==NULL ) {
		fprintf(stderr,"error: %s: cant open A(t)_ah_file <%s> for read\n", argv[0], argv[2]);
		exit(-1);
		}

	if( (f2=fopen(argv[3],"r"))==NULL ) {
		fprintf(stderr,"error: %s: cant open B(t)_ah_file <%s> for read\n", argv[0], argv[3]);
		exit(-1);
		}

	if( (f3=fopen(argv[4],"w"))==NULL ) {
		fprintf(stderr,"error: %s: cant open e(t)_ah_file <%s> for write\n", argv[0], argv[4]);
		exit(-1);
		}

	xdrstdio_create(&in1,  f1,  XDR_DECODE);
	xdrstdio_create(&in2,  f2,  XDR_DECODE);
	xdrstdio_create(&out, f3, XDR_ENCODE);

	fprintf(stdout,"t0\tt1\tt2\tt3\tt4\ta0\ta1\ta2\ta3\ta4\tb0\tb1\tb2\tb3\tb4\tR\n");

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

		first=TRUE;
		for( lag4p=4; lag4p<=maxlag; lag4p++ )
		for( lag3p=3; lag3p<lag4p;   lag3p++ )
		for( lag2p=2; lag2p<lag3p;   lag2p++ )
		for( lag1p=1; lag1p<lag2p;   lag1p++ ) {

			lag1=(lag1p); lag2=(lag2p); lag3=(lag3p); lag4=(lag4p);

			c[0]=data1; c[1]=data1; c[2]=data1; c[3]=data1;
			c[4]=data2; c[5]=data2; c[6]=data2; c[7]=data2; c[8]=data2;
			clag[0]=lag1; clag[1]=lag2; clag[2]=lag3; clag[3]=lag4;
			clag[4]=0; clag[5]=lag1; clag[6]=lag2; clag[7]=lag3; clag[8]=lag4;
			sn[0]=(1.0); sn[1]=(1.0); sn[2]=(1.0); sn[3]=(1.0);
			sn[4]=(-1.0); sn[5]=(-1.0); sn[6]=(-1.0); sn[7]=(-1.0); sn[8]=(-1.0);

			for( ii=0; ii<9; ii++ ) for (jj=ii; jj<9; jj++ ) {
				temp = lagdot( c[ii], clag[ii], c[jj], clag[jj], n );
				GTG[ii][jj] = sn[ii] * sn[jj] * temp;
				GTG[jj][ii] = GTG[ii][jj];
				}

			GTGmax = 0.0;
			for( ii=0; ii<9; ii++ ) for (jj=ii; jj<9; jj++ ) {
				temp = fabs(GTG[ii][jj]);
				if( temp>GTGmax ) GTGmax = temp;
				}

			if( GTGmax <= 0.0 ) GTGmax = 1.0;
			for( ii=0; ii<9; ii++) GTG[ii][ii] += damping*GTGmax;

			for( ii=0; ii<9; ii++ ) {
				GTd[ii] = (-sn[ii]) * lagdot( c[ii], clag[ii], data1, 0, n );
				}

			gauss(GTG,GTd,9,9,1.0e-6,&ierror,TRUE);

			E = 0.0; Ea = 0.0; Eb = 0.0;
			for( ii=0; ii<n; ii++ ) {
				e = data1[ii];
				ea = data1[ii]; eb=0.0;
				for( jj=0; jj<9; jj++ ) {
					kk = ii-clag[jj];
					if((kk>=0)&&(kk<n)) {
						e+=sn[jj]*GTd[jj]*(double)c[jj][kk];
						}
					}
				for( jj=0; jj<4; jj++ ) {
					kk = ii-clag[jj];
					if((kk>=0)&&(kk<n)) {
						ea+=GTd[jj]*(double)c[jj][kk];
						}
					}
				for( jj=4; jj<9; jj++ ) {
					kk = ii-clag[jj];
					if((kk>=0)&&(kk<n)) {
						eb+=GTd[jj]*(double)c[jj][kk];
						}
					}
				data3[ii]=(float)ea;
				data4[ii]=(float)eb;
				E += e*e;  Ea += ea*ea; Eb += eb*eb;
				}

			R = E/(Ea+Eb);

			if( first ) {
				for(ii=0;ii<9;ii++) soln[ii]=GTd[ii];
				Emin=E; Eamin=Ea; Ebmin=Eb; Rmin=R;
				lag1min=lag1; lag2min=lag2; lag3min=lag3; lag4min=lag4;
				first=FALSE;
				}
			else if( R<Rmin ) {
				for(ii=0;ii<9;ii++) soln[ii]=GTd[ii];
				Emin=E; Eamin=Ea; Ebmin=Eb; Rmin=R;
				lag1min=lag1; lag2min=lag2; lag3min=lag3; lag4min=lag4;
				}
			}

		E = 0.0; Ea = 0.0; Eb = 0.0;
		lag1=lag1min; lag2=lag2min; lag3=lag3min; lag4=lag4min;
		clag[0]=lag1; clag[1]=lag2; clag[2]=lag3; clag[3]=lag4;
		clag[4]=0; clag[5]=lag1; clag[6]=lag2; clag[7]=lag3; clag[8]=lag4;
		for( ii=0; ii<n; ii++ ) {
			e = data1[ii];
			ea = data1[ii]; eb=0.0;
			for( jj=0; jj<9; jj++ ) {
				kk = ii-clag[jj];
				if((kk>=0)&&(kk<n)) {
					e+=sn[jj]*soln[jj]*(double)c[jj][kk];
					}
				}
			for( jj=0; jj<4; jj++ ) {
				kk = ii-clag[jj];
				if((kk>=0)&&(kk<n)) {
					ea+=soln[jj]*(double)c[jj][kk];
					}
				}
			for( jj=4; jj<9; jj++ ) {
				kk = ii-clag[jj];
				if((kk>=0)&&(kk<n)) {
					eb+=soln[jj]*(double)c[jj][kk];
					}
				}
			data3[ii]=(float)ea;
			data4[ii]=(float)eb;
			E += e*e;  Ea += ea*ea; Eb += eb*eb;
			}

		R = E/(Ea+Eb);

fprintf(stdout,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%e\n", 0.0,
			head1.record.delta*(double)lag1min,
			head1.record.delta*(double)lag2min,
			head1.record.delta*(double)lag3min,
			head1.record.delta*(double)lag4min,
			soln[4], soln[5], soln[6], soln[7], soln[8],
			1.0, soln[0], soln[1], soln[2], soln[3], Rmin );

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


