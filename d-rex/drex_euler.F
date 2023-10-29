c
c
c these routines are from ccp4
c
c
C     POLARRFN
C     Copyright (C) 1986 Wolfgang Kabsch, Phil Evans
C
C     This code is distributed under the terms and conditions of the
C     CCP4 licence agreement as `Part 2' (Annex 2) software.
C     A copy of the CCP4 licence can be obtained by writing to the
C     CCP4 Secretary, Daresbury Laboratory, Warrington WA4 4AD, UK.
C
CC
      SUBROUTINE DCOSFD(ROT,DIRCOS,KAPPA)
C     ===================================
C******************************************************
C     Given a rotation matrix ROTN, expressing a rotation
C     through an angle KAPPA right-handedly about an axis
C     with direction cosines DIRCOS(I),
C     DCOSFD determines DIRCOS() and KAPPA
C     KAPPA is returned in degrees in the range 0 to 180
C
C     Expression for rotation matrix is (see J & J, "Mathl. Phys.",
C     p. 122):-
C
C  cw + n(1)n(1)(1-cw)    n(1)n(2)(1-cw)-n(3)sw  n(1)n(3)(1-cw)+n(3)sw
C  n(1)n(2)(1-cw)+n(3)sw  cw + n(2)n(2)(1-cw)    n(2)n(3)(1-cw)-n(1)sw
C  n(3)n(1)(1-cw)-n(2)sw  n(3)n(2)(1-cw)+n(1)sw  cw + n(3)n(3)(1-cw)
C
C     where cw = cos(KAPPA), sw = sin(KAPPA),
C     & n(1), n(2), n(3) are the direction cosines.
C
C*******************************************************}
C
      REAL ROT(3,3),DIRCOS(3),S(3),P(3),D(3)
      REAL KAPPA
      REAL DIFF,PI,RTD,TRACE,R2,R3,SINKAP,COSKAP,TERM
      INTEGER NS(3),NP(3),NCS(3),ND(3),KNTRL,J,JMX
      PARAMETER(PI=3.141592654D0, RTD=180./PI)
C
C
      DIFF=1E-04
      KNTRL=0
C     Trace = 1 + 2.cos w
      TRACE=ROT(1,1)+ROT(2,2)+ROT(3,3)
C     S(1)=2.n(1).sin w
C     S(2)=2.n(2).sin w
C     S(3)=2.n(3).sin w
      S(1)=ROT(3,2)-ROT(2,3)
      S(2)=ROT(1,3)-ROT(3,1)
      S(3)=ROT(2,1)-ROT(1,2)
      CALL ORDR3(KNTRL,S,NS)
C--Is biggest S() zero (i.e. is sin w = 0)?
      IF (ABS(S(NS(1))).GT.DIFF) THEN
        R2=S(NS(2))/S(NS(1))
        R3=S(NS(3))/S(NS(1))
        DIRCOS(NS(1))=1./SQRT(1. + R2**2 + R3**2)
        DIRCOS(NS(2))=R2*DIRCOS(NS(1))
        DIRCOS(NS(3))=R3*DIRCOS(NS(1))
        GOTO 123
      ENDIF
C--Calculation when sin w = 0 (esp. when w = 180').
C     P(1)=n(2).n(3).{0(W=0)or2(w=180)}
C     P(2)=n(3).n(1).{0(W=0)or2(w=180)}
C     P(3)=n(1).n(2).{0(W=0)or2(w=180)}
      P(1)=ROT(3,2)
      P(2)=ROT(1,3)
      P(3)=ROT(2,1)
      CALL ORDR3(KNTRL,P,NP)
C--   Is biggest P() zero (i.e. are all off-diag. terms zero)?
      IF (ABS(P(NP(1))).LT.DIFF) THEN
        IF (TRACE.GT.0.) THEN
C--Matrix is unit matrix.
          KAPPA=0.0
          DIRCOS(1)=0.
          DIRCOS(2)=0.
          DIRCOS(3)=1.
          RETURN
        ENDIF
C--Trace -ve, so dyad about x,y or z.
        KAPPA=PI
        DO J=1,3
          D(J)=ROT(J,J)+1.
        ENDDO
        CALL ORDR3(KNTRL,D,ND)
        DIRCOS(ND(1))=1.
        DIRCOS(ND(2))=0.
        DIRCOS(ND(3))=0.
        RETURN
      ENDIF
      IF (ABS(P(NP(2))).LT.DIFF) THEN
C--One d.c. is zero.
        DIRCOS(NP(1))=0.
        DIRCOS(NP(2))=SQRT(MAX(.5*(1.+ROT(NP(2),NP(2))),0.))
        DIRCOS(NP(3))=SQRT(MAX(.5*(1.+ROT(NP(3),NP(3))),0.))
        IF (P(NP(1)).LT.0.) DIRCOS(NP(2))=-DIRCOS(NP(2))
        GOTO 123
      ENDIF
C--All d.c's are nonzero.
C     R2=n(NP(1))/n(NP(2))
C     R3=n(NP(1))/n(NP(3))
      R2=P(NP(2))/P(NP(1))
      R3=P(NP(3))/P(NP(1))
      DIRCOS(NP(2))=1./SQRT(1. + R2**2 + (R2/R3)**2)
      DIRCOS(NP(3))=1./SQRT(1. + R3**2 + (R3/R2)**2)
C--Check: don't take sqrt of negative values
      TERM=1.-DIRCOS(NP(2))**2-DIRCOS(NP(3))**2
      IF (TERM.LT.0.) THEN
        IF (ABS(TERM).GT.DIFF)
     +  CALL CCPERR (1,'Negative value in DCOSFD')
        TERM=0.
      ENDIF
      DIRCOS(NP(1))=SQRT(TERM)
C--Adjust signs of DIRCOS().
      JMX=0
      DO J=1,3
        IF (P(J).GT.0.) THEN
          IF (JMX.GT.0) JMX=-1
          IF (JMX.EQ.0) JMX=J
        ENDIF
      ENDDO
      IF (JMX.GT.0) DIRCOS(JMX)=-DIRCOS(JMX)
      IF (JMX.EQ.0) WRITE(6,1358)
1358  FORMAT(' All P()s. are negative or zero')
C--Given d.c's, calculate KAPPA.
123   KNTRL=1
      CALL ORDR3(KNTRL,DIRCOS,NCS)
C     Find KAPPA.
      SINKAP=ROT(NCS(3),NCS(2))-ROT(NCS(2),NCS(3))
      SINKAP=.5*SINKAP/DIRCOS(NCS(1))
      COSKAP=.5*(TRACE-1.)
      KAPPA=RTD*ATAN2(SINKAP,COSKAP)
C  If kappa negative, negate kappa and invert DIRCOS
      IF (KAPPA.LT.0.) THEN
        KAPPA = -KAPPA
        DIRCOS(1) = -DIRCOS(1)
        DIRCOS(2) = -DIRCOS(2)
        DIRCOS(3) = -DIRCOS(3)
      ENDIF
      END
C
C
C
C
      SUBROUTINE POLMAT(ROTN,ANGLES)
C     ==============================
C
C**************************************************************
C       Sets rotation matrix ROTN expressing a rotation through an
C       angle AKAPPA right-handedly about an axis
C       with polar angles OMEGA, PHI (OMEGA with Z-axis, projection
C       giving angle PHI with X-axis).
C       These angles give direction cosines DIRCOS(I).
C***************************************************************
C
      REAL ROTN(3,3),DIRCOS(3),ANGLES(3)
      PARAMETER(PI=3.141592654D0, DTR=PI/180.)
C
      SNOM=SIN(DTR*ANGLES(1))
      CSOM=COS(DTR*ANGLES(1))
      SNPH=SIN(DTR*ANGLES(2))
      CSPH=COS(DTR*ANGLES(2))
      SNKA=SIN(DTR*ANGLES(3))
      CSKA=COS(DTR*ANGLES(3))
      DIRCOS(1)=SNOM*CSPH
      DIRCOS(2)=SNOM*SNPH
      DIRCOS(3)=CSOM
      DO I=1,3
        DO J=1,3
          K1=6-I-J
          K=K1
          IF (K1.LT.1 .OR. K1.GT.3) K=3
          EPSIJK=((I-J)*(J-K)*(K-I))/2
          ROTN(I,J)=DIRCOS(I)*DIRCOS(J)*(1.-CSKA)-EPSIJK*DIRCOS(K)*SNKA
          IF (I.EQ.J) ROTN(I,J)=ROTN(I,J)+CSKA
        ENDDO
      ENDDO
      END
