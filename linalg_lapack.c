/*


 linear algebra type routines that rely on lapack


$Id: linalg_lapack.c,v 1.1 2010/12/31 18:59:17 becker Exp $



*/
#include "fstrack.h"

/*

calculate exp(A t) using DGPADM from EXPOKIT for a 3x3

(works for both C and FORTRAN style)

*/
void calc_exp_matrixt(COMP_PRECISION *a,COMP_PRECISION t,
		      COMP_PRECISION *expout)
{
  static int ideg=8;// degre of Pade approximation, six should be ok
  static int m=3,ldh=3;// a(ldh,m) dimensions
  COMP_PRECISION *wsp;// work space
  int *ipiv,iexph,ns;// workspace, output pointer, nr of squareas
  int iflag,lwsp;// exit code, size of workspace
  int i,j,m2;
  if(fabs(t) < EPS_PREC){// return identity matrix in C style
    unity_matrix(expout,m);
  }else{
    /* 
       work space 
    */
    m2 = m*m;
    lwsp = 2*(4*m2+ideg+1);// size of workspace, oversized by factor two
    ipiv = (int *)malloc(sizeof(int)*m);
    my_vecalloc(&wsp,lwsp,"calc_exp_matrixt");
    if(!wsp || ! ipiv)MEMERROR("calc_exp_matrixt");
    //
    // call to expokit routine
    dgpadm(&ideg,&m,&t,a,&ldh,wsp,&lwsp,ipiv,&iexph,&ns,&iflag);
    if(iflag < 0){
      fprintf(stderr,"calc_exp_matrixt: problem in dgpadm: code: %i\n",
	      iflag);
      exit(-1);
    }
    // assign to output
    for(i=0,j=iexph-1;i<m2;i++,j++)
      expout[i] = wsp[j];
    free(wsp);free(ipiv);
  }
}


/*

  test if a m by m matrix a, is singular.  if the nullity, 0 <= nnull
  <= m, is non-zero the en[nnull*m] vector will hold the nnull
  eigenvectors[m] of the nullspace

  a will be destroyed

  wlim is the limit for singular value = 0 as a fraction of 
  the max singular value

*/
void test_null_space(COMP_PRECISION *a,COMP_PRECISION **en,
		     int m, int *nnull,COMP_PRECISION wlim)
{
  COMP_PRECISION *w, *v, *rv1,wmax;
  int i,mm;
  mm = m * m;
  //
  // do SVD decomposition
  // on output, the columns of V with zero singular value
  // will hold the an orthonormal basis of the null space
  //
  w=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*m);
  rv1=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*m);
  v=(COMP_PRECISION *)malloc(sizeof(COMP_PRECISION)*mm);
  if(!w || !rv1 || !v)
    MEMERROR("test_null_space");
  // this assumes that the A matrix is FORTRAN 
  // sorted
  svdcmp(a,&m,&m,&m,&m,w,v,rv1);
  free(rv1);
  *en=NULL;
  for(wmax=0.0,i=0;i<m;i++)	/* find max singular value */
    if(w[i] > wmax)
      wmax = w[i];
  wlim *= wmax;
  
  // search for zero singular values
  for((*nnull)=i=0;i<m;i++)
    if(w[i] < wlim){
      /* nnull-th eigenvector of the
	 null space is the i-th
	 column of the V matrix */
      *en=(COMP_PRECISION *)
	realloc(*en,(*nnull+1)*sizeof(COMP_PRECISION)*m);
      if(! *en)MEMERROR("test_null_space");
      copy_vector((v+i*m),(*en+(*nnull)*m),m);
      (*nnull) ++;
    }
  free(v);free(w);
}


/*

  given a m by m matrix a, will check if the matrix has full rank. if
  not, will return one of the b matrices (and allocate it) such that

  A * B = 0

  if nnull = 0 (not singular), b will not be allocated
  
  wlim should be eps small, say 1e-7

  
  this routine doesn't care if A is C or FORTRAN sorted

*/
void find_ABzero(COMP_PRECISION *a, COMP_PRECISION **b, int m,
		 int *nnull,COMP_PRECISION wlim)
{
  COMP_PRECISION *c,*en;
  int i,mm;
  mm=m * m;
  c=NULL;
  assign_c_for_ab_hs(a,&c,m); /* assemble the C[m*m,m*m] matrix
				 for C . x = 0 */
  en=NULL;test_null_space(c,&en,mm,nnull,wlim);
  free(c);
  if(*nnull){
    my_vecrealloc(b,mm,"find_ABzero");
    for(i=0;i<mm;i++)
      *(*b+i) = 0.0;
    for(i=0;i< *nnull;i++)
      add_a_to_b_vector((en+i*mm),*b,mm);
  }
  free(en);
}
