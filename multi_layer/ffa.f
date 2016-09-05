c
c aniosotropic reflectivity package from Vera Schulte-Pelkum in March 2005
c based on code by Peter Shearer and various others, see individual READMEs 
c and source code comments for details and description
c  
c minor modifications by Thorsten Becker
c
c $Id: ffa.f,v 1.1 2005/03/23 20:29:35 becker Exp twb $
c
c
cC-----------------------------------------------------------------      FFA
C SUBROUTINE:  FFA
C FAST FOURIER ANALYSIS SUBROUTINE
C------------------------------------------------------------------
C
      SUBROUTINE FFA(b , nfft , conj , norm)
C...
C...INPUT:
C........   nfft: INTEGER*4 number of points in the input time
C........         series; a power of 2 .LE. 16k.
C........
C........b(nfft): REAL*4  real-valued data series to be transformed.
C........
C........   conj: LOGICAL flag defining conjugation operation. If TRUE,
C........         the output spectrum is conjugated; c.f. notes below
C........         on the complex exponent convention.
C........
C........   norm: LOGICAL flag defining normalization operation. If TRUE,
C........         the output spectrum is normaliized by 1./nfft
C........
C...OUTPUT:
C........b(nfft): REAL*4 array containing the complex spectrum to be
C........         inverse transformed. The spectral component at
C........         frequency k is COMPLX( b(2*k+1) , b(2*k+2) ), except
C........         for the (real) Nyquist component which is packed
C........         into the (zero) imaginary part of the DC component.
C........         Thus, the spectrum consists of nfft reals rather
C........         than nfft+2.
C........
C....FFA calls the following procedures:
C...
C........R2TR   ===========================================================
C........R4TR
C........R8TR   (in //RITTER/USERS/USERLIB/SOURCE/OBS.DIR/FFA_FFS_SUBS.FTN)
C........ORD1
C........ORD2   ===========================================================
C........VEC_$MULT_CONSTANT
C........VEC_$MULT_CONSTANT_I
C...
C...     r.g.adair   5 Sep 83
C..........................................................................
c% include '/sys/ins/vec.ins.ftn'        apollo specific functions
      DIMENSION B(1)
      CHARACTER bell*1
      LOGICAL   conj,norm
      PARAMETER (k16=16*1024)
C
C
      bell=CHAR(7)
      IF(nfft.GT.k16) THEN
         PRINT *,bell//'(S.R. FFA) Your series is of length ',nfft
         PRINT *,'The input series must be .LE. 16k points long'
         STOP
      ENDIF
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
C DO A RADIX 2 OR RADIX 4 ITERATION FIRST IF ONE IS REQUIRED
C
      IF (M-N8POW*3-1) 50, 40, 30
  30  NN = 4
      INT = N/NN
      CALL R4TR(INT, B(1), B(INT+1), B(2*INT+1), B(3*INT+1))
      GO TO 60
  40  NN = 2
      INT = N/NN
      CALL R2TR(INT, B(1), B(INT+1))
      GO TO 60
  50  NN = 1
C
C PERFORM RADIX 8 ITERATIONS
C
  60  IF (N8POW) 90, 90, 70
  70  DO 80 IT=1,N8POW
        NN = NN*8
        INT = N/NN
        CALL R8TR(INT, NN, B(1), B(INT+1), B(2*INT+1), B(3*INT+1),
     *      B(4*INT+1), B(5*INT+1), B(6*INT+1), B(7*INT+1), B(1),
     *      B(INT+1), B(2*INT+1), B(3*INT+1), B(4*INT+1), B(5*INT+1),
     *      B(6*INT+1), B(7*INT+1))
  80  CONTINUE
C
C PERFORM IN-PLACE REORDERING
C
  90  CALL ORD1(M, B)
      CALL ORD2(M, B)
C
c apollo specific functions follow
c      IF(conj) THEN
c         CALL VEC_$MULT_CONSTANT_I(b(4),2,nfft/2 -1,-1.,b(4),2)
c      ENDIF
C
c      IF(norm) THEN
c         rfft=1./FLOAT(nfft)
c         CALL VEC_$MULT_CONSTANT(b(1),nfft,rfft,b(1))
c      ENDIF
C
      RETURN
      END                                                               
C
