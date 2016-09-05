c Subroutine SLOWNESS solves for the slowness, polarization,
c and group velocity vectors for a specified phase direction.
c For the case of general anisotropy, this is solved as an
c eigenvalue problem.
c  Given:
c    c(3,3,3,3) = elastic tensor
c    qhat(3)    = phase (slowness) direction.  This is a unit vector.
c  Returned:
c    qmag(3) = slowness magnitudes (eigenvalues)
c    q(3,3)  = slowness vectors
c    a(3,3)  = polarization vectors (eigenvectors)
c    umag(3) = group velocity magnitudes
c    u(3,3)  = group velocity vectors
c All eigenvalues and eigenvectors are sorted as
c      1) qP
c      2) qSV
c      3) qSH
c Requires: JACOBI
c
      subroutine SLOWNESS(c,qhat,qmag,q,a,umag,u)
      real c(3,3,3,3),qhat(3),qmag(3),q(3,3),a(3,3),umag(3),u(3,3)
      real m(3,3),h(3,3),eigen(3,3)
      integer index(3)
c
      do 40 i=1,3
         do 30 k=1,3
            m(i,k)=0.
            do 20 j=1,3
               do 10 l=1,3
10                m(i,k)=m(i,k)+c(i,j,k,l)*qhat(j)*qhat(l)
20          continue
            h(i,k)=m(i,k)
30       continue
40    continue
c      call DISP33(h)
c
      err=.00001
      call JACOBI(h,eigen,3,3,err)
c
c Now sort by polarization
c qP is largest eigenvalue (smallest slowness)
      emax=0.
      do 50 i=1,3
         if (h(i,i).gt.emax) then
            emax=h(i,i)
            index(1)=i
         end if
50    continue
c qSV has most vertical (001) polarization
      amax=0.
      do 60 i=1,3
         if (i.eq.index(1)) go to 60
         if (abs(eigen(3,i)).ge.amax) then
            amax=abs(eigen(3,i))
            index(2)=i
         end if
60    continue
c qSH is the remaining polarization
      do 70 i=1,3
         if (i.ne.index(1).and.i.ne.index(2)) then
            index(3)=i
         end if
70    continue
c
      do 90 i=1,3
         scr=h(index(i),index(i))
         if (scr.lt.0.) then
            print *,'We are going to bomb!!!!!'
            print *,'i,index= ',i,index
            print *,'h= ',h
         end if
         if (scr.ne.0.) then
            qmag(i)=1./sqrt(h(index(i),index(i)))
         else
            qmag(i)=0.
         end if
         do 80 j=1,3
            q(j,i)=qhat(j)*qmag(i)
            a(j,i)=eigen(j,index(i))
80       continue
90    continue
c
c Now find group velocity vectors
      do 150 i=1,3
         dotcor=0.
         do 130 j=1,3
            u(j,i)=0.
            do 120 k=1,3
            do 110 l=1,3
            do 100 n=1,3
               u(j,i)=u(j,i)+c(j,k,l,n)*a(k,i)*q(l,i)*a(n,i)
100         continue
110         continue
120         continue
            dotcor=dotcor+u(j,i)*q(j,i)
130      continue
         umag(i)=0.
         do 140 j=1,3
            if (dotcor.ne.0.) u(j,i)=u(j,i)/dotcor
            umag(i)=umag(i)+u(j,i)**2
140      continue
         umag(i)=sqrt(umag(i))
150   continue
c
      return
      end
c
c
c
c
      SUBROUTINE JACOBI(H,EIGEN,NN,MM,ERR)
C$$$$$ CALLS NO OTHER ROUTINES
C  JACOBI DIAGONALIZATION OF SYMMETRIC MATRICES.  SEE WILKINSON - THE
C  ALGEBRAIC EIGENVALUE PROBLEM, OUP 1964.
C  A VERSION FOR DEALING WITH HERMITIAN MATRICES IS AVAILABLE
C  H  IS A SYMMETRIC MATRIX OF ORDER  NN, DIMENSIONED H(MM,MM) IN
C  THE CALLING PROGRAM (MM.GE.NN).  ON RETURN  EIGEN  IS A REAL MATRIX
C  DIMENSIONED LIKE H, CONTAINING NORMALIZED EIGENVECTORS AS COLUMNS.
C  THE EIGENVALUES APPEAR ON THE DIAGONAL OF  H  , WHICH HAS BEEN OVER
C -WRITTEN DURING THE CALCULATION.
C  ERR  IS THE ERROR CRITERION. THE LARGEST OFF-DIAGONAL TERM WILL BE
C  AT LEAST  ERR  TIMES SMALLER THAN THE SUM OF MAGNITUDES OF EIGENVALS.
      DIMENSION H(MM,MM),EIGEN(MM,MM),HH(2)
      INTEGER P,Q,ROTS
C
      N=NN
      ROTS=0
C  EIGENVECTORS ARE FOUND AS THE PRODUCT MATRIX OF THE SEQUENCE OF
C  ELEMENTARY ROTATIONS THAT BRING  H  INTO DIAGONAL FORM.
      DO 110 J=1,N
      DO 100 I=1,N
100   EIGEN(I,J)=0.0
110   EIGEN(J,J)=1.0
      IF (N.EQ.1) RETURN
C  FIND THE  LARGEST  OFF-DIAGONAL ELEMNT AND PUT ITS POSITION IN P,Q
C  ALSO SUM THE DIAGONAL ELEMENTS  ABSOLUTE VALUES
 200  BIG=0.0
      DIAG=0.0
      DO 210 J=1,N
      DO 210 I=1,J
      HIJ=H(I,J)
      IF (I.EQ.J) GO TO 205
      HABS=ABS(HIJ)
      IF (BIG.GT.HABS) GO TO 210
      BIG=HABS
      P=I
      Q=J
      GO TO 210
205   DIAG=DIAG+ABS(HIJ)
210   CONTINUE
C  CHECKS ERROR CRITERION
      IF (BIG.GT.DIAG*ERR) GO TO 300
      RETURN
C
C  SETS UP VARIOUS CONSTANTS TO TRANSFORM CURRENT H
300   CONTINUE
      ROTS=ROTS+1
      HPQ=H(P,Q)
      HPP=H(P,P)
      HQQ=H(Q,Q)
      X=HPQ*SIGN(2.0,HPP-HQQ)
      Y=ABS(HPP-HQQ)
      D=0.5/SQRT(X*X+Y*Y)
      C2=0.5+Y*D
      C=SQRT(C2)
      CS=X*D
      S=CS/C
      S2=S*S
C
C  APPLIES ELEMENTARY UNITARY TRANSFORMS TO EIGEN, AND SIMILARITY TRANS
C -FORM TO H
      DO 310 I=1,N
      EIGNIP=EIGEN(I,P)
      EIGNIQ=EIGEN(I,Q)
      EIGEN(I,P)=EIGNIP*C+EIGNIQ*S
      EIGEN(I,Q)=EIGNIQ*C-EIGNIP*S
      IF (I.EQ.Q .OR. I.EQ.P) GO TO 310
      HIP=H(I,P)
      HIQ=H(I,Q)
      H(I,P)=HIP*C+HIQ*S
      H(P,I)=H(I,P)
      H(I,Q)=HIQ*C-HIP*S
      H(Q,I)=H(I,Q)
310   CONTINUE
C
C  DEALS WITH THE SPECIAL ELEMENTS AT P,P Q,Q AND P,Q
      H(P,P)=HPP*C2+HQQ*S2+2.0*CS*HPQ
      H(Q,Q)=HPP*S2+HQQ*C2-2.0*CS*HPQ
      H(P,Q)=0.0
      H(Q,P)=0.0
C  REPEATS ITERATION
      GO TO 200
      END
