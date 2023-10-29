#define TRUE (1)
#define FALSE (0)

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <rpc/rpc.h>
#include "SHahhead.h"

char *progname;

main(int argc, char *argv[])
{
	FILE *f1, *f2, *f3;
	XDR in1, in2, out;
	ahhed head1, head2;
	float *data1, *data2, *data3, *data4;
	int i, j, n, ii, jj;
	double E, e, Emin, ea, eb, Ea, Eb, R, tlag;
	double xx, xy, a0, b0;

	progname=argv[0];

	if( argc != 5 ) {
		fprintf(stderr,"usage: %s max_lag A(t)_ah_file B(t)_ah_file e(t)_ah_file > results.txt\n", argv[0]);
		fprintf(stderr,"       The parameter max_lag (in seconds) is provided for\n");
		fprintf(stderr,"       compatibility with other software and is in fact not used.\n");
		fprintf(stderr,"\n");
		fprintf(stderr,"purpose: estimate trivial 1-pulse operator that relates two\n");
		fprintf(stderr,"    timeseries A(t) and B(t)\n");
		fprintf(stderr,"\n");
		fprintf(stderr,"    Suppose A(t) = s(t)*a(t) and B(t) = s(t)*b(t)\n");
		fprintf(stderr,"    where   a(t) = a0 delta(t-t0) and\n");
		fprintf(stderr,"            b(t) = b0 delta(t-t0)\n");
		fprintf(stderr,"    and where s(t) is unknown.  This program determines\n");
		fprintf(stderr,"    the amplitudes a0 and b0 (up to a multiplicative constant,\n");
		fprintf(stderr,"    e.g. b0=1) by a least-squares fit of a0 and b0. The\n");
		fprintf(stderr,"    time lag t0 is constrained to be zero.  The fundamental\n");
		fprintf(stderr,"    idea is to minimize the norm, R, of the prediction\n");
		fprintf(stderr,"    error, e(t) = A(t)*b(t) - B(t)*a(t), where\n");
		fprintf(stderr,"    Here R = ||e||**2 / [ ||A(t)*b(t)||**2 + ||B(t)*a(t)||**2 ]\n");
		fprintf(stderr,"    The best fitting t0, a0, b0, R are written to stdout\n");
		fprintf(stderr,"    and A(t)*b(t) and B(t)*a(t) written to the e(t)_ah_file.\n");
		fprintf(stderr,"\n");
		fprintf(stderr,"    When using this routine to to solve for receiver function operator\n");
		fprintf(stderr,"    I recommend that A(t) represent the radial-horizontal component\n");
		fprintf(stderr,"    and B(t) represent the vertical component.\n");
		exit(-1);
		}

	
	fprintf(stdout,"t0\ta0\tb0\tR\n");

	if( sscanf(argv[1],"%le",&tlag) != 1 ) {
		fprintf(stderr,"error: %s: can read max_lag from command line\n", argv[0] );
		exit(-1);
		}

	if( (f1=fopen(argv[2],"r"))==NULL ) {
		fprintf(stderr,"error: %s: cant open A(t)_ah_file <%s> for read\n", argv[0], argv[1]);
		exit(-1);
		}

	if( (f2=fopen(argv[3],"r"))==NULL ) {
		fprintf(stderr,"error: %s: cant open B(t)_ah_file <%s> for read\n", argv[0], argv[2]);
		exit(-1);
		}

	if( (f3=fopen(argv[4],"w"))==NULL ) {
		fprintf(stderr,"error: %s: cant open e(t)_ah_file <%s> for write\n", argv[0], argv[3]);
		exit(-1);
		}

	xdrstdio_create(&in1,  f1,  XDR_DECODE);
	xdrstdio_create(&in2,  f2,  XDR_DECODE);
	xdrstdio_create(&out,  f3, XDR_ENCODE);

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

		if( (data3=(float*)malloc(n*sizeof(float))) == NULL ) {
			fprintf(stderr,"error: %s: out of memory!\n", argv[0]);
			exit(-1);
			}
		if( (data4=(float*)malloc(n*sizeof(float))) == NULL ) {
			fprintf(stderr,"error: %s: out of memory!\n", argv[0]);
			exit(-1);
			}

		xx=0.0; xy=0.0;
		for( ii=0; ii<n; ii++ ) {
			xx+=data1[ii]*data1[ii];
			xy+=data1[ii]*data2[ii];
			}
		a0=xx/xy;
		b0=1.0;

		E = 0.0; Ea=0.0; Eb=0.0;
		for( ii=0; ii<n; ii++ ) {
			e = b0*data1[ii]-a0*data2[ii];
			ea = b0*data1[ii];
			eb = a0*data2[ii];
			data3[ii]=(float)ea;
			data4[ii]=(float)eb;
			E += e*e; Ea += ea*ea; Eb += eb*eb;
			}
		R = E/(Ea+Eb);

		fprintf(stdout,"%f\t%f\t%f\t%e\n", 0.0, a0, b0, R);

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

	exit(0);
	}
