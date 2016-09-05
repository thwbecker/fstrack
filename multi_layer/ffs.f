c
c aniosotropic reflectivity package from Vera Schulte-Pelkum in March 2005
c based on code by Peter Shearer and various others, see individual READMEs 
c and source code comments for details and description
c  
c minor modifications by Thorsten Becker
c
c $Id: ffs.f,v 1.2 2005/03/24 23:38:11 becker Exp twb $
c
c
cC-----------------------------------------------------------------------
C SUBROUTINE:  FFS
C FAST FOURIER SYNTHESIS SUBROUTINE
C RADIX 8-4-2
C-----------------------------------------------------------------------
C
      SUBROUTINE FFS (b , nfft , conj , norm)                          
C...INPUT:
C........   nfft: INTEGER*4 number of points in the output time
C........         series; a power of 2 .LE. 16k.
C........
C........b(nfft): REAL*4 array containing the complex spectrum to be
C........         inverse transformed. The spectral component at
C........         frequency k is COMPLX( b(2*k+1) , b(2*k+2) ), except
C........         for the (real) Nyquist component which is packed
C........         into the (zero) imaginary part of the DC component.
C........         Thus, the spectrum consists of nfft reals rather
C........         than nfft+2.
C........
C........   conj: LOGICAL flag defining conjugation operation. If TRUE,
C........         the spectrum is conjugated before inverse-transforming;
C........         c.f. notes below on the complex exponent convention.
C........
C........   norm: LOGICAL flag defining normalization operation. If TRUE,
C........         the output time series is normaliized by 1./nfft
C........
C...OUTPUT:
C........b(nfft): the real-valued inverse transform of the input spectrum.
C...
C....FFS calls the following procedures:
C...
C........R2TR   ===========================================================
C........R4SYN
C........R8SYN  (in //RITTER/USERS/USERLIB/SOURCE/OBS.DIR/FFA_FFS_SUBS.FTN)
C........ORD1
C........ORD2   ===========================================================
C........VEC_$MULT_CONSTANT
C........VEC_$MULT_CONSTANT_I
C...
C...     r.g.adair   5 Sep 83
C...............................................................................
c% include '/sys/ins/vec.ins.ftn'
      real b
      DIMENSION B(1)
      CHARACTER bell*1
      LOGICAL   conj,norm
      PARAMETER (k16=16*1024)
C
C
      bell=CHAR(7)
      IF(nfft.GT.k16) THEN
         PRINT *,bell//'(S.R. FFS) Your series is of length ',nfft
         PRINT *,'The input series must be .LE. 16k points long'
         STOP
      ENDIF
C
C
      N = 1
      DO 10 I=1,15
        M = I
        N = N*2
        IF (N.EQ.NFFT) GO TO 20
  10  CONTINUE
      PRINT *,bell//'ERROR: your series length (',nfft,
     : ') is not a power of two'
      STOP
  20  CONTINUE
      N8POW = M/3
C
c      IF(conj) THEN
c         CALL VEC_$MULT_CONSTANT_I(b(4),2,nfft/2-1,-1.,b(4),2)
c      ENDIF
C
C REORDER THE INPUT FOURIER COEFFICIENTS
C
      CALL ORD2(M, B)
      CALL ORD1(M, B)
C
      IF (N8POW.EQ.0) GO TO 60
C
C PERFORM THE RADIX 8 ITERATIONS
C
      NN = N
      DO 50 IT=1,N8POW
        INT = N/NN
        CALL R8SYN(INT, NN, B, B(INT+1), B(2*INT+1), B(3*INT+1),
     *      B(4*INT+1), B(5*INT+1), B(6*INT+1), B(7*INT+1), B(1),
     *      B(INT+1), B(2*INT+1), B(3*INT+1), B(4*INT+1), B(5*INT+1),
     *      B(6*INT+1), B(7*INT+1))
        NN = NN/8
  50  CONTINUE
C
C DO A RADIX 2 OR RADIX 4 ITERATION IF ONE IS REQUIRED
C
  60  IF (M-N8POW*3-1) 90, 80, 70
  70  INT = N/4
      CALL R4SYN(INT, B(1), B(INT+1), B(2*INT+1), B(3*INT+1))
      GO TO 90
  80  INT = N/2
      CALL R2TR(INT, B(1), B(INT+1))
C
90    continue
c90    IF(norm) THEN
c         rfft=1./FLOAT(nfft)
c         CALL VEC_$MULT_CONSTANT(b(1),nfft,rfft,b(1))
c      ENDIF
C
      RETURN
      END                                                               


