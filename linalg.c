/*


miscellaneous linear algebra type routines


$Id: linalg.c,v 1.20 2005/10/19 22:52:57 becker Exp becker $



*/
#include "fstrack.h"
#include <math.h>
/*
  
  given real, symmetric a 3 x 3 matrix A where only the RR, RT, RP, TT, TP, and PP
  components have to be filled (the corresponding other ones will be overwritten)

  since this is for symmteric matrices, doesn't matter if we have C or
  FORTRAN storage

  calculate the eigenvalues eval[3]
  and the corresponding eigenvectors, evec is 3x3

  val will have the eigenvalues in ascending order, 
  e1 = val[2],e2 = val[1],e3 = val[0]
  or 
  e1 = val[E1],e2 = val[E2],e3 = val[E3]

  with e1 > e2 > e3

  vec will be each normalized to unity 

  r, theta, phi
  
  vec[6,7,8] is the vector that corresponds to the highest  eigenvalue,     val[2]
  vec[3,4,5] is the vector that corresponds to the intermediate eigenvalue, val[1]
  vec[0,1,2] is the vector that corresponds to the smallest eigenvalue,     val[0]


  WARNING:

  if icalc_vectors is not set to TRUE, will only calculate the eigenvalues

  uses EISPACK routine

  we have also defined macros E1, E2, E3 that refer to the indices
  2, 1, and 0, respectively (see above)

*/
void calc_eigensystem_sym(COMP_PRECISION *a,COMP_PRECISION *eval,
			  COMP_PRECISION *evec,
			  my_boolean icalc_vectors)
{
  static int n=3;// dimension, don't use a define since we want to pass n to 'rs'
  int i,ierr,matz;
  COMP_PRECISION fv1[3],fv2[3],loca[9];
  //
  // assign the symmetric values of the matrix
  // this is unnecessary, do it for safety
  // (will, however, overwrite the presumed symmetric entries)
  //
  a[TR]=a[RT];a[PR]=a[RP];a[PT]=a[TP];
  copy_3dmatrix(a,loca);// save the original A
  //
  // EISPACK matz flag, 0: only eigenvalues, !=0: values + vectors
  matz=(icalc_vectors)?(1):(0);
  // call the appropriate symmetric matrix eigensystem routine
  // from EISPACK
  SYMREAL_ES_ROUTINE(&n,&n,loca,eval,&matz,evec,fv1,fv2,&ierr);
  if(ierr){
    fprintf(stderr,"calc_eigensystem_sym: runtime error %i in routine\n",ierr);
    exit(-1);
  }
  for(i=0;i<3;i++)
    if(!finite(eval[i]))
      fprintf(stderr,"calc_eigensystem_sym: WARNING: eigenvalue %i not finite\n",i);
}
//
// FORTRAN wrapper (ie. to be called from FORTRAN routines) for above routine
void calc_eigensystem_sym_(COMP_PRECISION *a,COMP_PRECISION *eval,
			   COMP_PRECISION *evec,int icalc_vectors)
{
  if(icalc_vectors)
    calc_eigensystem_sym(a,eval,evec,TRUE);
  else
    calc_eigensystem_sym(a,eval,evec,FALSE);
}
/*
  
  given real general 3 x 3 matrix A 
  
  IN FORTRAN STORAGE

  calculate the eigenvalues eval[3]

  modified from EISPACK docu:

  c        evalr  and  evali  contain the real and imaginary parts,
  c        respectively, of the eigenvalues.  complex conjugate
  c        pairs of eigenvalues appear consecutively with the
  c        eigenvalue having the positive imaginary part first.
  c
  c        evec[3x3] 
  c        contains the real and imaginary parts of the eigenvectors
  c        if matz is not zero.  if the j-th eigenvalue is real, the
  c        j-th column of  z  contains its eigenvector.  if the j-th
  c        eigenvalue is complex with positive imaginary part, the
  c        j-th and (j+1)-th columns of  z  contain the real and
  c        imaginary parts of its eigenvector.  the conjugate of this
  c        vector is the eigenvector for the conjugate eigenvalue.
  c
  c        ierr  is an integer output variable set equal to an error
  c           completion code described in the documentation for hqr
  c           and hqr2.  the normal completion code is zero.
  c
  c        iv1  and  fv1  are temporary storage arrays.
  c

*/
void calc_eigensystem(COMP_PRECISION *a,COMP_PRECISION *evalr,
		      COMP_PRECISION *evali,
		      COMP_PRECISION *evec,
		      my_boolean icalc_vectors,
		      my_boolean *imaginary)
{
  static int n=3;// dimension, don't use a define since we want to pass n to 'rs'
  int i,ierr,matz;
  int iv1[3];

  COMP_PRECISION fv1[3],loca[9];
  copy_vector(a,loca,9);// save the original A
  // EISPACK matz flag, 0: only eigenvalues, !=0: values + vectors
  matz=(icalc_vectors)?(1):(0);
  // call the appropriate symmetric matrix eigensystem routine
  // from EISPACK
  REAL_ES_ROUTINE(&n,&n,loca,evalr,evali,&matz,
		  evec,iv1,fv1,&ierr);
  if(ierr){
    fprintf(stderr,"calc_eigensystem: runtime error %i in routine\n",ierr);
    exit(-1);
  }
  *imaginary = FALSE;
  for(i=0;i < n;i++){
    if(fabs(evali[i]) >= EPS_COMP_PREC)
      *imaginary = TRUE;
    if(!finite(hypot(evalr[i],evali[i])))
      fprintf(stderr,"calc_eigensystem: WARNING: |eigenvalue| %i not finite\n",i);
  }
}
/*

  calculate the eigenvectors and eigenvalues of a 2-D symmetric matrix
  S given as s11, s12, and s22 output: eigenvalues e1 > e2 and azi,
  the angle in degrees clockwise from 2-dir (North)

*/
void eigen2d(COMP_PRECISION s11,COMP_PRECISION s12, 
	     COMP_PRECISION s22,COMP_PRECISION *e1,
	     COMP_PRECISION *e2, COMP_PRECISION *azi)
{
  COMP_PRECISION x1,x2,r;
  x1 = (s11 + s22)/2.0;
  x2 = (s11 - s22)/2.0;
  r = x2 * x2 + s12 * s12;
  if(r > 0.0){
    r = sqrt(r);
    *e1 = x1 + r;
    *e2 = x1 - r;
  }else{
    *e1 = *e2 = x1;
  }
  if(fabs(x2)>EPS_PREC)
    *azi = 28.64788975654116044 * atan2(s12,x2);// factor is 90/pi
  else if(s12 <= 0.0)
    *azi = -45.0;
  else
    *azi =  45.0;
  *azi = 90.0 - *azi;
}


/*

  calculate the sqrt of a symmetric tensor by doing a EV decomposition like 
  A = X \Lambda X^-1 = X \Lambda X^T and taking the sqrt of it's 
  eigenvalues in \Lambda and reassembling 

  this function uses the save_sqrt function, ie. small entries
  will be thrown away

*/
void calc_sqrt_sym_matrix(COMP_PRECISION *l2, COMP_PRECISION *l)
{
  COMP_PRECISION eval[3],evec[9];
  int i;
  calc_eigensystem_sym(l2,eval,evec,TRUE);
  for(i=0;i<3;i++)
    eval[i] = save_sqrt(eval[i]);
  assemble_sym_matrix(evec,eval,l);
}

/*

  assemble a symmetric matrix based on the orthogonal and normalized
  eigenvectors and the eigenvalues from a decomposition

  (for symmetric matrices, X^-1 == X^T, where X is the matrix that holds the 
  eigenvectors)

  output is the matrix A

  used to obtain L from L^2 

*/
void assemble_sym_matrix(COMP_PRECISION *evec,COMP_PRECISION *eval,
			 COMP_PRECISION *a)
{
  COMP_PRECISION x[9];
  int i,j,k,os1,os2;
  // do C  = Lambda dot X^T
  for(i=os1=0;i<3;i++,os1+=3)
    for(j=0;j<3;j++){
      x[os1+j]=0.0;
      for(k=0;k<3;k++)
	x[os1+j] += ((i==k)?(eval[i]):0.0) * evec[k*3+j];
    }
  // do X dot C
  for(i=os1=0;i<3;i++,os1+=3)
    for(j=0;j<3;j++){
      a[os1+j]=0.0;
      for(k=os2=0;k<3;k++,os2+=3)
	a[os1+j] += evec[os2+i] * x[os2+j];
    }
}
/*

  calculate 

  B = A dot A^T

  for 3x3 matrices stored in [9] format

*/
void calc_a_dot_at(COMP_PRECISION *a, COMP_PRECISION *b)
{
  int j,k,i3,j3,index;
  for(i3=0;i3<9;i3+=3){
    for(j3=j=0;j < 3;j++,j3+=3){
      index=i3+j;
      for(b[index]=0.0,k=0;k < 3;k++)
	b[index] += a[i3+k]*a[j3+k];
    }
  }
}
/*

  calculate 

  B = A dot A^T

  for 3x3 matrices stored in [3][3] format

*/
void calc_a_dot_at_3x3(COMP_PRECISION a[3][3], 
		       COMP_PRECISION b[3][3])
{
  int i,j,k;
  for(i=0;i<3;i++)
    for(j=0;j<3;j++){
      b[i][j] = 0.0;
      for(k=0;k<3;k++)
	b[i][j] += a[i][k] * a[j][k];
    }
}
/*

  calculate 

  B = A^T dot A 

  for 3x3 matrices
*/
void calc_at_dot_a(COMP_PRECISION *a, COMP_PRECISION *b)
{
  int index,i,j,k;
  for(index=i=0;i<3;i++)
    for(j=0;j < 3;j++,index++){
      for(b[index]=0.0,k=0;k < 9;k+=3)
	b[index] += a[k+i]*a[k+j];
    }
}
/*

  calculate the symmetric part of a cartesian 
  decomposition of a matrix
  b = 0.5 (a + a^T)

  (C style)
*/
void calc_cd_symm_part(COMP_PRECISION *a, COMP_PRECISION *b)
{
  int i,i3,j,j3,index;
  for(i=i3=0;i<3;i++,i3+=3)
    for(j=j3=0;j<3;j++,j3+=3){
      index = i3 + j;
      b[index] = 0.5 * (a[index] + a[j3+i]);
    }
}

/*

  calculate the asymmetric part of a cartesian 
  decomposition of a matrix
  b = 0.5 (a - a^T)

  (C style)
*/
void calc_cd_asymm_part(COMP_PRECISION *a, COMP_PRECISION *b)
{
  int i,i3,j,j3,index;
  for(i=i3=0;i<3;i++,i3+=3)
    for(j=j3=0;j<3;j++,j3+=3){
      index = i3 + j;
      b[index] = 0.5 * (a[index] - a[j3+i]);
    }
}
/* 
  
   return the determinant of a 3x3 matrix that is input in vector storage
   0 1 2  rr rt rp
   3 4 5  tr tt tp
   6 7 8  pr pt pp
*/
COMP_PRECISION det3x3(COMP_PRECISION *a)
{
  
  return 
      a[0]*a[4]*a[8]+a[1]*a[5]*a[6]
    + a[2]*a[3]*a[7]-a[0]*a[5]*a[7]
    - a[1]*a[3]*a[8]-a[2]*a[4]*a[6];
}
/* 

   compute the trace of a 3x3 matrix

*/
COMP_PRECISION trace3x3(COMP_PRECISION *a)
{
  return a[RR]+a[TT]+a[PP];
}
/*

  return the second invariant of a symmetric 3x3 matrix in C style

  I_2 =    a_xx a_yy + a_xx a_zz + a_yy a_zz 
         - a_xy^2    - a_xz^2    - a_yz^2

*/
COMP_PRECISION sec_inv3x3(COMP_PRECISION *a)
{
  return a[RR]*a[TT] + a[RR]*a[PP] + a[TT]*a[PP] -
    a[RT]*a[RT] - a[RP]*a[RP] - a[TP]*a[TP];
}
/*

  returns the double dot product

  A : B = a_11 b_11 + a_12 b_21 + a_13 b_31 + 
          a_21 b_21 + a_22 b_22 + a_23 b_23 +
	  a_31 b_31 + a_32 b_32 + a_33 b_33
*/
COMP_PRECISION double_dot(COMP_PRECISION *a,COMP_PRECISION *b)
{
  COMP_PRECISION tmp;
  int i,j,i3,j3;
  tmp = 0.0;
  for(i=i3=0;i<3;i++,i3+=3)
    for(j=j3=0;j<3;j++,j3+=3)
      tmp += a[i3+j] * b[j3+i];
  return tmp;
}
/*


calculate C = A . B for 3x3 matrices stored in 9 format

*/
void calc_a_times_b_mat(COMP_PRECISION *a, COMP_PRECISION *b,
			COMP_PRECISION *c)
{
  int i3,j,k,k3;
  for(i3=0;i3<9;i3+=3)
    for(j=0;j<3;j++){
      c[i3+j] = 0.0;
      for(k=k3=0;k<3;k++,k3+=3)
	c[i3+j] += a[i3+k] * b[k3+j];
    }
}
/*


calculate C = A . B for 3x3 matrices stored in 3x3 format

*/
void calc_a_times_b_3x3mat(COMP_PRECISION a[3][3],
			   COMP_PRECISION b[3][3],
			   COMP_PRECISION c[3][3])
{
  int i,j,k;
  for(i=0;i<3;i++)
    for(j=0;j<3;j++){
      c[i][j] = 0.0;
      for(k=0;k<3;k++)
	c[i][j] += a[i][k] * b[k][j];
    }
}




/*

low-level, BLAS type calculations

*/
COMP_PRECISION vector_diff(COMP_PRECISION *x,COMP_PRECISION *y,
			   int n)
{
  int i;
  COMP_PRECISION sum,tmp;
  sum = 0.0;
  for(i=0;i<n;i++){
    tmp = x[i]-y[i];
    sum += tmp * tmp;
  }
  if(sum > 0)
    return sqrt(sum);
  else 
    return 0.0;
}
void zero_vector(COMP_PRECISION *x,int n)
{
  int i;
  for(i=0;i<n;i++)
    x[i]=0.0;
}
void scale_vector(COMP_PRECISION *x,COMP_PRECISION fac,int n)
{
  int i;
  for(i=0;i<n;i++)
    x[i] *= fac;
}
void scale_vector3d(COMP_PRECISION *x,COMP_PRECISION fac)
{
  x[0] *= fac;  x[1] *= fac;  x[2] *= fac;
}
// b = a 
void copy_vector(COMP_PRECISION *a, COMP_PRECISION *b, int n)
{
  a_equals_b_vector(b,a,n);
}
void copy_3dvector(COMP_PRECISION *a, COMP_PRECISION *b)
{
  a_equals_b_vector3d(b,a);
}
void copy_3dmatrix(COMP_PRECISION *a, COMP_PRECISION *b)
{
  a_equals_b_vector(b,a,9);
}
void add_a_to_b_vector(COMP_PRECISION *a, COMP_PRECISION *b,
		       int n)
{
  int i;
  for(i=0;i<n;i++)
    b[i] += a[i];
}
void sub_a_from_b_vector(COMP_PRECISION *a, COMP_PRECISION *b,
			 int n)
{
  int i;
  for(i=0;i<n;i++)
    b[i] -= a[i];
}
void a_equals_b_minus_c_vector(COMP_PRECISION *a, COMP_PRECISION *b, COMP_PRECISION *c,
			       int n)
{
  int i;
  for(i=0;i<n;i++)
    a[i] = b[i] - c[i];
}
//
// 3-D vector dot product
//
COMP_PRECISION vec3ddotp(COMP_PRECISION *a,COMP_PRECISION *b)
{
  COMP_PRECISION tmp;
  tmp  = a[0] * b[0];
  tmp += a[1] * b[1];
  tmp += a[2] * b[2];
  return tmp;
}
//
// 2-D vector dot product
//
COMP_PRECISION vec2ddotp(COMP_PRECISION *a,COMP_PRECISION *b)
{
  return a[0] * b[0] + a[1] * b[1];
}
COMP_PRECISION norm3d(COMP_PRECISION *a)
{
  COMP_PRECISION xl;
  xl = vec3ddotp(a,a);
  if(xl <= 0.0)
    return 0.0;
  else
    return sqrt(xl);
}
/* normalize a vector */
void normalize3d(COMP_PRECISION *x)
{
  COMP_PRECISION length;
  length = norm3d(x);
  scale_vector3d(x,1.0/length);
}
/* return the sum of the elements of a vector */
COMP_PRECISION vec_sum(COMP_PRECISION *x, int n)
{
  COMP_PRECISION tmp;
  int i;
  for(tmp = 0.0, i=0;i<n;i++)
    tmp += x[i];
  return tmp;
}
/*

  given a n by n matrix

  sorting, will allocate and create a c matrix which is n*n by n*n so
  that the solution x[n*n] of

  C . x = 0 

  is the matrix B that solves

  A * B = 0

 */
void assign_c_for_ab_hs(COMP_PRECISION *a,COMP_PRECISION **c, 
			int n)
{
  int i,j,k,irow,icol,m;
  m = n * n;
  *c=(COMP_PRECISION *)realloc(*c,sizeof(COMP_PRECISION)*m*m);
  for(i=0;i<m;i++)
    *(*c+i) = 0.0;
  // assign to C
  irow = 0;
  for(i=0;i < n;i++){
    for(j=0;j < n;j++){
      for(k=0;k < n;k++){
	icol = j * n + k;// this is B(k,j)
	*(*c+ icol*m + irow) = a[i + k*n];// this is A(i,k)
      }      
      irow++;
    }
  }
}
/* 
   using a symmetric 3x3 matrix stored in [9] format, go 
   to [3][3] format 
*/
void symmat9to3x3(COMP_PRECISION *l, COMP_PRECISION b[3][3])
{
  b[FSTRACK_R][FSTRACK_R] = l[RR];
  b[FSTRACK_R][FSTRACK_THETA] = b[FSTRACK_THETA][FSTRACK_R] = l[RT];
  b[FSTRACK_R][FSTRACK_PHI] = b[FSTRACK_PHI][FSTRACK_R] = l[RP];
  b[FSTRACK_THETA][FSTRACK_THETA] = l[TT];
  b[FSTRACK_THETA][FSTRACK_PHI] = b[FSTRACK_PHI][FSTRACK_THETA] = l[TP];
  b[FSTRACK_PHI][FSTRACK_PHI] = l[PP];
}
/* convert a matrix stored in 0...8 format to [3][3] */
void mat9to3x3(COMP_PRECISION *a,COMP_PRECISION b[3][3])
{
  b[FSTRACK_R][FSTRACK_R] =     a[RR];    b[FSTRACK_R][FSTRACK_THETA] = a[RT]; b[FSTRACK_R][FSTRACK_PHI]     = a[RP];
  b[FSTRACK_THETA][FSTRACK_R] = a[TR];b[FSTRACK_THETA][FSTRACK_THETA] = a[TT]; b[FSTRACK_THETA][FSTRACK_PHI] = a[TP];
  b[FSTRACK_PHI][FSTRACK_R] =   a[PR];  b[FSTRACK_PHI][FSTRACK_THETA] = a[PT]; b[FSTRACK_PHI][FSTRACK_PHI]   = a[PP];
}
/* go the other way */
void mat3x3to9(COMP_PRECISION a[3][3],COMP_PRECISION *b)
{
  b[RR]=a[FSTRACK_R][FSTRACK_R];    b[RT]=a[FSTRACK_R][FSTRACK_THETA];      b[RP] = a[FSTRACK_R][FSTRACK_PHI];
  b[TR]=a[FSTRACK_THETA][FSTRACK_R];b[TT]=a[FSTRACK_THETA][FSTRACK_THETA] ; b[TP] = a[FSTRACK_THETA][FSTRACK_PHI] ;
  b[PR]=a[FSTRACK_PHI][FSTRACK_R]; b[PT]=a[FSTRACK_PHI][FSTRACK_THETA] ;   b[PP] = a[FSTRACK_PHI][FSTRACK_PHI];
}
/*
  find the largest element of a vector

*/
COMP_PRECISION max_vec(COMP_PRECISION *x, int n)
{
  int i;
  COMP_PRECISION max;
  max = x[0];
  for(i=1;i<n;i++)
    if(x[i] > max)
      max = x[i];
  return max;
}

/*

find the maximum of the absolute values of a vector

*/

COMP_PRECISION max_abs_vec(COMP_PRECISION *x, int n)
{
  int i;
  COMP_PRECISION max=0.0,tmp;
  for(i=0;i<n;i++){
    tmp = fabs(x[i]);
    if(tmp > max)
      max = tmp;
  }
  return max;
}

//  \vec{a} = \vec{b}
void a_equals_b_vector(COMP_PRECISION *a,COMP_PRECISION *b,
		       int n)
{
  memcpy(a,b,n*sizeof(COMP_PRECISION));
}
void a_equals_b_vector3d(COMP_PRECISION *a,COMP_PRECISION *b)
{
  memcpy(a,b,3*sizeof(COMP_PRECISION));
}
/* 
   
given a x[n][n] C style matrix, return the unity matrix

*/
void unity_matrix(COMP_PRECISION *x,int n)
{
  int i,j;
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      x[i*n+j] = (i==j)?(1.0):(0.0);
}

/* 
   compute the transpose of B and assign it to A

   A(3,3) = B(3,3)^T

   a and b are passed as 9 element vectors
 */
void a_is_b_transpose_3by3(COMP_PRECISION *a,COMP_PRECISION *b)
{
  COMP_PRECISION c[9];
  c[0] = b[0];c[1] = b[3];c[2] = b[6];
  c[3] = b[1];c[4] = b[4];c[5] = b[7];
  c[6] = b[2];c[7] = b[5];c[8] = b[8];
  a_equals_b_vector(a,c,9);	/* to allow calling as (x,x) */
}
void transpose_3by3_inplace(COMP_PRECISION *a)
{
  COMP_PRECISION c[9];
  c[0] = a[0];c[1] = a[3];c[2] = a[6];
  c[3] = a[1];c[4] = a[4];c[5] = a[7];
  c[6] = a[2];c[7] = a[5];c[8] = a[8];
  a_equals_b_vector(a,c,9);	/* to allow calling as (x,x) */
}

void a_is_b_transpose_nbyn(COMP_PRECISION *a,COMP_PRECISION *b,int n)
{
  COMP_PRECISION *c;
  int nn,j,i;
  nn = n * n;
  my_vecalloc(&c,nn,"a_is_b_transposenbyn");
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      c[j*n+i] = b[i*n+j];
  a_equals_b_vector(a,c,nn);	/* to allow calling as (x,x) */
  free(c);
}
/* cross product in 2-D  */
COMP_PRECISION vec2d_crossp(COMP_PRECISION *a,COMP_PRECISION *b)
{
  return a[0]*b[1] -a[1]*b[0];
}
void transpose_nbyn_inplace(COMP_PRECISION *a,int n)
{
  COMP_PRECISION *c;
  int nn,j,i;
  nn = n * n;
  my_vecalloc(&c,nn,"transpose_nbyn_inplace");
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      c[j*n+i] = a[i*n+j];
  a_equals_b_vector(a,c,nn);	
  free(c);
}


/* norm */
COMP_PRECISION norm(COMP_PRECISION *x, int n)
{
  return sqrt(vecdotp(x,x,n));
}
/* dot product */
COMP_PRECISION vecdotp(COMP_PRECISION *a,COMP_PRECISION *b,
		       int n)
{
  int i;
  COMP_PRECISION x;
  x = 0.0;
  for(i=0;i < n;i++)
    x += a[i] * b[i];
  return x;
}
void zero_small_entries(COMP_PRECISION *x,int n)
{
  int i;
#ifndef SINGLE_PRECISION
  for(i=0;i<n;i++){
    if(fabs(x[i]) < 5e-14)
      x[i] = 0.0;
  }
#else
  for(i=0;i<n;i++){
    if(fabs(x[i]) < 5e-7)
      x[i] = 0.0;
  }
#endif
}
/* remove the trace of a 3 x 3 matrix given as 9 vector */
COMP_PRECISION remove_trace(COMP_PRECISION *a)
{
  COMP_PRECISION div;
  div = (a[XX] + a[YY] + a[ZZ])/3.0;
  a[XX] -= div;    
  a[YY] -= div;    
  a[ZZ] -= div;
  return div;
}
