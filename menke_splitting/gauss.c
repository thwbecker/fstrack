#include <stdio.h>
#include <malloc.h>
#include <ctype.h>
#include <math.h>

gauss(a,vec,n,nstore,test,ierror,itriag)
double *a, vec[], test;
int n, nstore, *ierror, itriag;
{
 
/* subroutine gauss, by william menke */
/* july 1978 (modified feb 1983, nov 85) */
 
/* a subroutine to solve a system of n linear equations in n unknowns*/
/* where n doesn't exceed MAX_TABLE_COLS */
/* gaussian reduction with partial pivoting is used */
/* 	a		(sent, destroyed)	n by n matrix		*/
/*	vec		(sent, overwritten)	n vector, replaced w/ solution*//* 	nstore		(sent)			dimension of a	*/
/*	test		(sent)			div by zero check number*/
/*	ierror		(returned)		zero on no error*/
/*	itriag		(sent)			matrix triangularized only*/
/*						 on TRUE useful when solving*/
/*						 multiple systems with same a */	static int isub[100], l1;
	int line[100], iet, ieb, i, j, k, l, j2;
	double big, testa, b, sum;
	

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
                                       testa=(double) fabs(
						(double) (*(a+l1*nstore+j)) );
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
 
                       sum=1.0/(*(a+i*nstore+j)); 
				/*reduce matrix towards triangle */
		       for( k=0; k<n; k++ ) {
				if( line[k]==0 ) {
					b=(*(a+k*nstore+j))*sum;
				       	for( l=j+1; l<n; l++ ) {
                                               *(a+k*nstore+l)=
							(*(a+k*nstore+l))
							-b*(*(a+i*nstore+l));
					       } /*end for l*/
                                       *(a+k*nstore+j)=b;
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
 
	for( j=0; j<n-1; j++) { /*transform the vector to match triang. matrix*/               b=vec[isub[j]];
               for( k=0; k<n; k++ ) {
                      if (line[k]>j) {	/* skip elements outside of triangle*/
                                vec[k]=vec[k]-(*(a+k*nstore+j))*b;
				} /*end if*/
			} /*end for k*/
		} /*end for j*/
 
      b = *(a+l1*nstore+(n-1));   /*apex of triangle*/
      if( ((double)fabs( (double) b))<=test) {
		/*check for div by zero in backsolving*/
		ieb=2;
		} /*end if*/
      vec[isub[n-1]]=vec[isub[n-1]]/b;
 
      for( j=n-2; j>=0; j-- ) {	/* backsolve rest of triangle*/
		sum=vec[isub[j]];
		for( j2=j+1; j2<n; j2++ ) {
			sum = sum - vec[isub[j2]] * (*(a+isub[j]*nstore+j2));
			} /*end for j2*/
			b = *(a+isub[j]*nstore+j);
               if( ((double)fabs((double)b))<=test) {
			/* test for div by 0 in backsolving */
			ieb=2;
			} /*end if*/
		vec[isub[j]]=sum/b;   /*solution returned in vec*/
		} /*end for j*/

/*put the solution vector into the proper order*/

      for( i=0; i<n; i++ ) {    /* reorder solution */
		for( k=i; k<n; k++ ) {  /* search for i-th solution element */
			if( line[k]==i ) {
				j=k;
				break;
				} /*end if*/
			} /*end for k*/
               b = vec[j];       /* swap solution and pointer elements*/
               vec[j] = vec[i];
               vec[i] = b;
               line[j] = line[i];
		} /*end for i*/
 
      *ierror = iet + ieb;   /* set final error flag*/
}


