c
c aniosotropic reflectivity package from Vera Schulte-Pelkum in March 2005
c based on code by Peter Shearer and various others, see individual READMEs 
c and source code comments for details and description
c  
c minor modifications by Thorsten Becker
c
c $Id: utilmat.f,v 1.1 2005/03/23 20:29:35 becker Exp $
c
c
cc ZERO66 fills a 6x6 matrix with zeros
      SUBROUTINE ZERO66(c)
      implicit double precision (a-h,o-z)
      dimension c(6,6)
      do 20 i=1,6
      do 10 j=1,6
      c(i,j)=0.0d0
10    continue
20    continue
      return
      end
c
c
c DISP66 displays a 6x6 matrix
      SUBROUTINE DISP66(c)
      implicit double precision (a-h,o-z)
      dimension c(6,6)
      do 10 j=1,6
10       print 20,c(j,1),c(j,2),c(j,3),c(j,4),c(j,5),c(j,6)
20    format (6f10.4)
      return
      end
c
c
c DISP33 displays a 3x3 matrix
      SUBROUTINE DISP33(c)
      implicit double precision (a-h,o-z)
      dimension c(3,3)
      do 10 j=1,3
10       print 20,c(j,1),c(j,2),c(j,3)
20    format (3f10.4)
      return
      end
c
c
c subroutine MMULT multiplies two square matrices
      SUBROUTINE MMULT(a,b,c,n)
      implicit double precision (a-h,o-z)
      dimension a(n,n),b(n,n),c(n,n)
      do 40 i=1,n
      do 30 j=1,n
         c(i,j)=0.0d0
         do 20 k=1,n
20          c(i,j)=c(i,j)+a(i,k)*b(k,j)
30    continue
40    continue
      return
      end
c
c
c
c
      SUBROUTINE INVERT(A, N, DETERM)
C$$$$$ CALLS NO OTHER ROUTINES
C  MODIFIED VERSION OF UCSD MATINV,  TO INVERT SQUARE MATRIX  A  IN PLACE.
C  THE ARRAY  A  MUST BE N BY N CLOSE PACKED, WITH N.LE.30,  DETERMINANT
C  APPEARS IN  DETERM. SINGULAR MATRICES ARE DETECTED BY ZERO DETERMINANT
      implicit double precision (a-h,o-z)
      DIMENSION A(N,N),IPIVOT(30),INDEX(30,2),PIVOT(30)
      EQUIVALENCE(IROW,JROW),(ICOLUM,JCOLUM),(AMAX,T,SWAP)
      DETERM=1.0d0
      DO 20 J=1,N
   20 IPIVOT(J)=0.0d0
      DO 550 I=1,N
      AMAX=0.0d0
      DO 105 J=1,N
      IF(IPIVOT(J)-1)60,105,60
   60 DO 100 K=1,N
      IF(IPIVOT(K)-1)80,100,740
   80 IF(ABS(AMAX)-ABS(A(J,K)))85,83,83
 83   IF (.NOT.(AMAX.EQ.0.0d0 .AND. A(J,K).EQ.0.0d0)) GO TO 100
   85 IROW=J
      ICOLUM=K
      AMAX=A(J,K)
  100 CONTINUE
  105 CONTINUE
      IPIVOT(ICOLUM)=IPIVOT(ICOLUM)+1
      IF(IROW-ICOLUM)140,260,140
  140 DETERM=-DETERM
      DO 200 L=1,N
      SWAP=A(IROW,L)
      A(IROW,L)=A(ICOLUM,L)
  200 A(ICOLUM,L)=SWAP
  260 INDEX(I,1)=IROW
      INDEX(I,2)=ICOLUM
      PIVOT(I)=A(ICOLUM,ICOLUM)
      DETERM=DETERM*PIVOT(I)
      IF (DETERM.EQ.0.0d0) RETURN
      A(ICOLUM,ICOLUM)=1.0d0
      DO 350 L=1,N
  350 A(ICOLUM,L)=A(ICOLUM,L)/PIVOT(I)
      DO 550 L1=1,N
      IF(L1-ICOLUM)400,550,400
  400 T=A(L1,ICOLUM)
      A(L1,ICOLUM)=0.0d0
      DO 450 L=1,N
  450 A(L1,L)=A(L1,L)-A(ICOLUM,L)*T
  550 CONTINUE
      DO 710 I=1,N
      L=N+1-I
      IF(INDEX(L,1)-INDEX(L,2))630,710,630
  630 JROW=INDEX(L,1)
      JCOLUM=INDEX(L,2)
      DO 705 K=1,N
      SWAP=A(K,JROW)
      A(K,JROW)=A(K,JCOLUM)
      A(K,JCOLUM)=SWAP
  705 CONTINUE
  710 CONTINUE
  740 RETURN
      END                                                                !INVERT
