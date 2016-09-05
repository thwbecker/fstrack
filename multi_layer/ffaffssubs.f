c
c aniosotropic reflectivity package from Vera Schulte-Pelkum in March 2005
c based on code by Peter Shearer and various others, see individual READMEs 
c and source code comments for details and description
c  
c minor modifications by Thorsten Becker
c
c $Id: ffaffssubs.f,v 1.1 2005/03/23 20:29:35 becker Exp twb $
c
c
cC-----------------------------------------------------------------------
C SUBROUTINE:  R2TR                                                     R2TR
C RADIX 2 ITERATION SUBROUTINE
C-----------------------------------------------------------------------
C
C
      SUBROUTINE R2TR(INT, B0, B1)
      DIMENSION B0(1), B1(1)
      DO 10 K=1,INT
        T = B0(K) + B1(K)
        B1(K) = B0(K) - B1(K)
        B0(K) = T
  10  CONTINUE
      RETURN
      END                                                               
C
C-----------------------------------------------------------------------
C SUBROUTINE:  R4TR                                                     R4TR
C RADIX 4 ITERATION SUBROUTINE
C-----------------------------------------------------------------------
C
      SUBROUTINE R4TR(INT, B0, B1, B2, B3)
      DIMENSION B0(1), B1(1), B2(1), B3(1)
      DO 10 K=1,INT
        R0 = B0(K) + B2(K)
        R1 = B1(K) + B3(K)
        B2(K) = B0(K) - B2(K)
        B3(K) = B1(K) - B3(K)
        B0(K) = R0 + R1
        B1(K) = R0 - R1
  10  CONTINUE
      RETURN
      END                                                              
C
C-----------------------------------------------------------------------
C SUBROUTINE: R8TR                                                      R8TR
C RADIX 8 ITERATION SUBROUTINE
C-----------------------------------------------------------------------
C
      SUBROUTINE R8TR(INT, NN, BR0, BR1, BR2, BR3, BR4, BR5, BR6, BR7,
     *    BI0, BI1, BI2, BI3, BI4, BI5, BI6, BI7)
      DIMENSION L(15), BR0(1), BR1(1), BR2(1), BR3(1), BR4(1), BR5(1),
     *    BR6(1), BR7(1), BI0(1), BI1(1), BI2(1), BI3(1), BI4(1),
     *    BI5(1), BI6(1), BI7(1)
C.........................................................................
      REAL*4 pii,pi2,pi8,c22,s22,p7,p7two
      PARAMETER (pii=3.1415927,pi2=pii*2,pi8=pii/8,p7=0.7071067)
      PARAMETER (p7two=2*p7,c22=0.9238795,s22=0.3826834)
C.........................................................................
      EQUIVALENCE (L15,L(1)), (L14,L(2)), (L13,L(3)), (L12,L(4)),
     *    (L11,L(5)), (L10,L(6)), (L9,L(7)), (L8,L(8)), (L7,L(9)),
     *    (L6,L(10)), (L5,L(11)), (L4,L(12)), (L3,L(13)), (L2,L(14)),
     *    (L1,L(15))
C
C SET UP COUNTERS SUCH THAT JTHET STEPS THROUGH THE ARGUMENTS
C OF W, JR STEPS THROUGH STARTING LOCATIONS FOR THE REAL PART OF THE
C INTERMEDIATE RESULTS AND JI STEPS THROUGH STARTING LOCATIONS
C OF THE IMAGINARY PART OF THE INTERMEDIATE RESULTS.
C
      L(1) = NN/8
      DO 40 K=2,15
        IF (L(K-1)-2) 10, 20, 30
  10    L(K-1) = 2
  20    L(K) = 2
        GO TO 40
  30    L(K) = L(K-1)/2
  40  CONTINUE
      PIOVN = PII/FLOAT(NN)
      JI = 3
      JL = 2
      JR = 2
      DO 120 J1=2,L1,2
      DO 120 J2=J1,L2,L1
      DO 120 J3=J2,L3,L2
      DO 120 J4=J3,L4,L3
      DO 120 J5=J4,L5,L4
      DO 120 J6=J5,L6,L5
      DO 120 J7=J6,L7,L6
      DO 120 J8=J7,L8,L7
      DO 120 J9=J8,L9,L8
      DO 120 J10=J9,L10,L9
      DO 120 J11=J10,L11,L10
      DO 120 J12=J11,L12,L11
      DO 120 J13=J12,L13,L12
      DO 120 J14=J13,L14,L13
      DO 120 JTHET=J14,L15,L14
        TH2 = JTHET - 2
        IF (TH2) 50, 50, 90
  50    DO 60 K=1,INT
          T0 = BR0(K) + BR4(K)
          T1 = BR1(K) + BR5(K)
          T2 = BR2(K) + BR6(K)
          T3 = BR3(K) + BR7(K)
          T4 = BR0(K) - BR4(K)
          T5 = BR1(K) - BR5(K)
          T6 = BR2(K) - BR6(K)
          T7 = BR3(K) - BR7(K)
          BR2(K) = T0 - T2
          BR3(K) = T1 - T3
          T0 = T0 + T2
          T1 = T1 + T3
          BR0(K) = T0 + T1
          BR1(K) = T0 - T1
          PR = P7*(T5-T7)
          PI = P7*(T5+T7)
          BR4(K) = T4 + PR
          BR7(K) = T6 + PI
          BR6(K) = T4 - PR
          BR5(K) = PI - T6
  60    CONTINUE
        IF (NN-8) 120, 120, 70
  70    K0 = INT*8 + 1
        KL = K0 + INT - 1
        DO 80 K=K0,KL
          PR = P7*(BI2(K)-BI6(K))
          PI = P7*(BI2(K)+BI6(K))
          TR0 = BI0(K) + PR
          TI0 = BI4(K) + PI
          TR2 = BI0(K) - PR
          TI2 = BI4(K) - PI
          PR = P7*(BI3(K)-BI7(K))
          PI = P7*(BI3(K)+BI7(K))
          TR1 = BI1(K) + PR
          TI1 = BI5(K) + PI
          TR3 = BI1(K) - PR
          TI3 = BI5(K) - PI
          PR = TR1*C22 - TI1*S22
          PI = TI1*C22 + TR1*S22
          BI0(K) = TR0 + PR
          BI6(K) = TR0 - PR
          BI7(K) = TI0 + PI
          BI1(K) = PI - TI0
          PR = -TR3*S22 - TI3*C22
          PI = TR3*C22 - TI3*S22
          BI2(K) = TR2 + PR
          BI4(K) = TR2 - PR
          BI5(K) = TI2 + PI
          BI3(K) = PI - TI2
  80    CONTINUE
        GO TO 120
  90    ARG = TH2*PIOVN
        C1 = COS(ARG)
        S1 = SIN(ARG)
        C2 = C1**2 - S1**2
        S2 = C1*S1 + C1*S1
        C3 = C1*C2 - S1*S2
        S3 = C2*S1 + S2*C1
        C4 = C2**2 - S2**2
        S4 = C2*S2 + C2*S2
        C5 = C2*C3 - S2*S3
        S5 = C3*S2 + S3*C2
        C6 = C3**2 - S3**2
        S6 = C3*S3 + C3*S3
        C7 = C3*C4 - S3*S4
        S7 = C4*S3 + S4*C3
        INT8 = INT*8
        J0 = JR*INT8 + 1
        K0 = JI*INT8 + 1
        JLAST = J0 + INT - 1
        DO 100 J=J0,JLAST
          K = K0 + J - J0
          TR1 = BR1(J)*C1 - BI1(K)*S1
          TI1 = BR1(J)*S1 + BI1(K)*C1
          TR2 = BR2(J)*C2 - BI2(K)*S2
          TI2 = BR2(J)*S2 + BI2(K)*C2
          TR3 = BR3(J)*C3 - BI3(K)*S3
          TI3 = BR3(J)*S3 + BI3(K)*C3
          TR4 = BR4(J)*C4 - BI4(K)*S4
          TI4 = BR4(J)*S4 + BI4(K)*C4
          TR5 = BR5(J)*C5 - BI5(K)*S5
          TI5 = BR5(J)*S5 + BI5(K)*C5
          TR6 = BR6(J)*C6 - BI6(K)*S6
          TI6 = BR6(J)*S6 + BI6(K)*C6
          TR7 = BR7(J)*C7 - BI7(K)*S7
          TI7 = BR7(J)*S7 + BI7(K)*C7
C
          T0 = BR0(J) + TR4
          T1 = BI0(K) + TI4
          TR4 = BR0(J) - TR4
          TI4 = BI0(K) - TI4
          T2 = TR1 + TR5
          T3 = TI1 + TI5
          TR5 = TR1 - TR5
          TI5 = TI1 - TI5
          T4 = TR2 + TR6
          T5 = TI2 + TI6
          TR6 = TR2 - TR6
          TI6 = TI2 - TI6
          T6 = TR3 + TR7
          T7 = TI3 + TI7
          TR7 = TR3 - TR7
          TI7 = TI3 - TI7
C
          TR0 = T0 + T4
          TI0 = T1 + T5
          TR2 = T0 - T4
          TI2 = T1 - T5
          TR1 = T2 + T6
          TI1 = T3 + T7
          TR3 = T2 - T6
          TI3 = T3 - T7
          T0 = TR4 - TI6
          T1 = TI4 + TR6
          T4 = TR4 + TI6
          T5 = TI4 - TR6
          T2 = TR5 - TI7
          T3 = TI5 + TR7
          T6 = TR5 + TI7
          T7 = TI5 - TR7
          BR0(J) = TR0 + TR1
          BI7(K) = TI0 + TI1
          BI6(K) = TR0 - TR1
          BR1(J) = TI1 - TI0
          BR2(J) = TR2 - TI3
          BI5(K) = TI2 + TR3
          BI4(K) = TR2 + TI3
          BR3(J) = TR3 - TI2
          PR = P7*(T2-T3)
          PI = P7*(T2+T3)
          BR4(J) = T0 + PR
          BI3(K) = T1 + PI
          BI2(K) = T0 - PR
          BR5(J) = PI - T1
          PR = -P7*(T6+T7)
          PI = P7*(T6-T7)
          BR6(J) = T4 + PR
          BI1(K) = T5 + PI
          BI0(K) = T4 - PR
          BR7(J) = PI - T5
 100    CONTINUE
        JR = JR + 2
        JI = JI - 2
        IF (JI-JL) 110, 110, 120
 110    JI = 2*JR - 1
        JL = JR
 120  CONTINUE
      RETURN
      END                                                              
C
C
C-----------------------------------------------------------------------
C SUBROUTINE:  R4SYN                                                    R4SYN
C RADIX 4 SYNTHESIS
C-----------------------------------------------------------------------
C
      SUBROUTINE R4SYN(INT, B0, B1, B2, B3)
      DIMENSION B0(1), B1(1), B2(1), B3(1)
      DO 10 K=1,INT
        T0 = B0(K) + B1(K)
        T1 = B0(K) - B1(K)
        T2 = B2(K) + B2(K)
        T3 = B3(K) + B3(K)
        B0(K) = T0 + T2
        B2(K) = T0 - T2
        B1(K) = T1 + T3
        B3(K) = T1 - T3
  10  CONTINUE
      RETURN
      END                                                               
C
C-----------------------------------------------------------------------
C SUBROUTINE:  R8SYN                                                    R8SYN
C RADIX 8 SYNTHESIS SUBROUTINE
C-----------------------------------------------------------------------
C
      SUBROUTINE R8SYN(INT, NN, BR0, BR1, BR2, BR3, BR4, BR5, BR6, BR7,
     *    BI0, BI1, BI2, BI3, BI4, BI5, BI6, BI7)
      DIMENSION L(15), BR0(1), BR1(1), BR2(1), BR3(1), BR4(1), BR5(1),
     *    BR6(1), BR7(1), BI0(1), BI1(1), BI2(1), BI3(1), BI4(1),
     *    BI5(1), BI6(1), BI7(1)
C.........................................................................
      REAL*4 pii,pi2,pi8,c22,s22,p7,p7two
      PARAMETER (pii=3.1415927,pi2=pii*2,pi8=pii/8,p7=0.7071067)
      PARAMETER (p7two=2*p7,c22=0.9238795,s22=0.3826834)
C.........................................................................
      EQUIVALENCE (L15,L(1)), (L14,L(2)), (L13,L(3)), (L12,L(4)),
     *    (L11,L(5)), (L10,L(6)), (L9,L(7)), (L8,L(8)), (L7,L(9)),
     *    (L6,L(10)), (L5,L(11)), (L4,L(12)), (L3,L(13)), (L2,L(14)),
     *    (L1,L(15))
      L(1) = NN/8
      DO 40 K=2,15
        IF (L(K-1)-2) 10, 20, 30
  10    L(K-1) = 2
  20    L(K) = 2
        GO TO 40
  30    L(K) = L(K-1)/2
  40  CONTINUE
      PIOVN = PII/FLOAT(NN)
      JI = 3
      JL = 2
      JR = 2
C
      DO 120 J1=2,L1,2
      DO 120 J2=J1,L2,L1
      DO 120 J3=J2,L3,L2
      DO 120 J4=J3,L4,L3
      DO 120 J5=J4,L5,L4
      DO 120 J6=J5,L6,L5
      DO 120 J7=J6,L7,L6
      DO 120 J8=J7,L8,L7
      DO 120 J9=J8,L9,L8
      DO 120 J10=J9,L10,L9
      DO 120 J11=J10,L11,L10
      DO 120 J12=J11,L12,L11
      DO 120 J13=J12,L13,L12
      DO 120 J14=J13,L14,L13
      DO 120 JTHET=J14,L15,L14
        TH2 = JTHET - 2
        IF (TH2) 50, 50, 90
  50    DO 60 K=1,INT
          T0 = BR0(K) + BR1(K)
          T1 = BR0(K) - BR1(K)
          T2 = BR2(K) + BR2(K)
          T3 = BR3(K) + BR3(K)
          T4 = BR4(K) + BR6(K)
          T6 = BR7(K) - BR5(K)
          T5 = BR4(K) - BR6(K)
          T7 = BR7(K) + BR5(K)
          PR = P7*(T7+T5)
          PI = P7*(T7-T5)
          TT0 = T0 + T2
          TT1 = T1 + T3
          T2 = T0 - T2
          T3 = T1 - T3
          T4 = T4 + T4
          T5 = PR + PR
          T6 = T6 + T6
          T7 = PI + PI
          BR0(K) = TT0 + T4
          BR1(K) = TT1 + T5
          BR2(K) = T2 + T6
          BR3(K) = T3 + T7
          BR4(K) = TT0 - T4
          BR5(K) = TT1 - T5
          BR6(K) = T2 - T6
          BR7(K) = T3 - T7
  60    CONTINUE
        IF (NN-8) 120, 120, 70
  70    K0 = INT*8 + 1
        KL = K0 + INT - 1
        DO 80 K=K0,KL
          T1 = BI0(K) + BI6(K)
          T2 = BI7(K) - BI1(K)
          T3 = BI0(K) - BI6(K)
          T4 = BI7(K) + BI1(K)
          PR = T3*C22 + T4*S22
          PI = T4*C22 - T3*S22
          T5 = BI2(K) + BI4(K)
          T6 = BI5(K) - BI3(K)
          T7 = BI2(K) - BI4(K)
          T8 = BI5(K) + BI3(K)
          RR = T8*C22 - T7*S22
          RI = -T8*S22 - T7*C22
          BI0(K) = (T1+T5) + (T1+T5)
          BI4(K) = (T2+T6) + (T2+T6)
          BI1(K) = (PR+RR) + (PR+RR)
          BI5(K) = (PI+RI) + (PI+RI)
          T5 = T1 - T5
          T6 = T2 - T6
          BI2(K) = P7TWO*(T6+T5)
          BI6(K) = P7TWO*(T6-T5)
          RR = PR - RR
          RI = PI - RI
          BI3(K) = P7TWO*(RI+RR)
          BI7(K) = P7TWO*(RI-RR)
  80    CONTINUE
        GO TO 120
  90    ARG = TH2*PIOVN
        C1 = COS(ARG)
        S1 = -SIN(ARG)
        C2 = C1**2 - S1**2
        S2 = C1*S1 + C1*S1
        C3 = C1*C2 - S1*S2
        S3 = C2*S1 + S2*C1
        C4 = C2**2 - S2**2
        S4 = C2*S2 + C2*S2
        C5 = C2*C3 - S2*S3
        S5 = C3*S2 + S3*C2
        C6 = C3**2 - S3**2
        S6 = C3*S3 + C3*S3
        C7 = C3*C4 - S3*S4
        S7 = C4*S3 + S4*C3
        INT8 = INT*8
        J0 = JR*INT8 + 1
        K0 = JI*INT8 + 1
        JLAST = J0 + INT - 1
        DO 100 J=J0,JLAST
          K = K0 + J - J0
          TR0 = BR0(J) + BI6(K)
          TI0 = BI7(K) - BR1(J)
          TR1 = BR0(J) - BI6(K)
          TI1 = BI7(K) + BR1(J)
          TR2 = BR2(J) + BI4(K)
          TI2 = BI5(K) - BR3(J)
          TR3 = BI5(K) + BR3(J)
          TI3 = BI4(K) - BR2(J)
          TR4 = BR4(J) + BI2(K)
          TI4 = BI3(K) - BR5(J)
          T0 = BR4(J) - BI2(K)
          T1 = BI3(K) + BR5(J)
          TR5 = P7*(T0+T1)
          TI5 = P7*(T1-T0)
          TR6 = BR6(J) + BI0(K)
          TI6 = BI1(K) - BR7(J)
          T0 = BR6(J) - BI0(K)
          T1 = BI1(K) + BR7(J)
          TR7 = -P7*(T0-T1)
          TI7 = -P7*(T1+T0)
          T0 = TR0 + TR2
          T1 = TI0 + TI2
          T2 = TR1 + TR3
          T3 = TI1 + TI3
          TR2 = TR0 - TR2
          TI2 = TI0 - TI2
          TR3 = TR1 - TR3
          TI3 = TI1 - TI3
          T4 = TR4 + TR6
          T5 = TI4 + TI6
          T6 = TR5 + TR7
          T7 = TI5 + TI7
          TTR6 = TI4 - TI6
          TI6 = TR6 - TR4
          TTR7 = TI5 - TI7
          TI7 = TR7 - TR5
          BR0(J) = T0 + T4
          BI0(K) = T1 + T5
          BR1(J) = C1*(T2+T6) - S1*(T3+T7)
          BI1(K) = C1*(T3+T7) + S1*(T2+T6)
          BR2(J) = C2*(TR2+TTR6) - S2*(TI2+TI6)
          BI2(K) = C2*(TI2+TI6) + S2*(TR2+TTR6)
          BR3(J) = C3*(TR3+TTR7) - S3*(TI3+TI7)
          BI3(K) = C3*(TI3+TI7) + S3*(TR3+TTR7)
          BR4(J) = C4*(T0-T4) - S4*(T1-T5)
          BI4(K) = C4*(T1-T5) + S4*(T0-T4)
          BR5(J) = C5*(T2-T6) - S5*(T3-T7)
          BI5(K) = C5*(T3-T7) + S5*(T2-T6)
          BR6(J) = C6*(TR2-TTR6) - S6*(TI2-TI6)
          BI6(K) = C6*(TI2-TI6) + S6*(TR2-TTR6)
          BR7(J) = C7*(TR3-TTR7) - S7*(TI3-TI7)
          BI7(K) = C7*(TI3-TI7) + S7*(TR3-TTR7)
 100    CONTINUE
        JR = JR + 2
        JI = JI - 2
        IF (JI-JL) 110, 110, 120
 110    JI = 2*JR - 1
        JL = JR
 120  CONTINUE
      RETURN
      END                                                              
C
C-----------------------------------------------------------------------
C SUBROUTINE:  ORD1                                                     ORD1
C IN-PLACE REORDERING SUBROUTINE
C-----------------------------------------------------------------------
C
      SUBROUTINE ORD1(M, B)
      DIMENSION B(1)
C
      K = 4
      KL = 2
      N = 2**M
      DO 40 J=4,N,2
        IF (K-J) 20, 20, 10
  10    T = B(J)
        B(J) = B(K)
        B(K) = T
  20    K = K - 2
        IF (K-KL) 30, 30, 40
  30    K = 2*J
        KL = J
  40  CONTINUE
      RETURN
      END                                                               
C
C-----------------------------------------------------------------------
C SUBROUTINE:  ORD2                                                     ORD2
C IN-PLACE REORDERING SUBROUTINE
C-----------------------------------------------------------------------
C
      SUBROUTINE ORD2(M, B)
      DIMENSION L(15), B(1)
      EQUIVALENCE (L15,L(1)), (L14,L(2)), (L13,L(3)), (L12,L(4)),
     *    (L11,L(5)), (L10,L(6)), (L9,L(7)), (L8,L(8)), (L7,L(9)),
     *    (L6,L(10)), (L5,L(11)), (L4,L(12)), (L3,L(13)), (L2,L(14)),
     *    (L1,L(15))
      N = 2**M
      L(1) = N
      DO 10 K=2,M
        L(K) = L(K-1)/2
  10  CONTINUE
      DO 20 K=M,14
        L(K+1) = 2
  20  CONTINUE
      IJ = 2
      DO 40 J1=2,L1,2
      DO 40 J2=J1,L2,L1
      DO 40 J3=J2,L3,L2
      DO 40 J4=J3,L4,L3
      DO 40 J5=J4,L5,L4
      DO 40 J6=J5,L6,L5
      DO 40 J7=J6,L7,L6
      DO 40 J8=J7,L8,L7
      DO 40 J9=J8,L9,L8
      DO 40 J10=J9,L10,L9
      DO 40 J11=J10,L11,L10
      DO 40 J12=J11,L12,L11
      DO 40 J13=J12,L13,L12
      DO 40 J14=J13,L14,L13
      DO 40 JI=J14,L15,L14
        IF (IJ-JI) 30, 40, 40
  30    T = B(IJ-1)
        B(IJ-1) = B(JI-1)
        B(JI-1) = T
        T = B(IJ)
        B(IJ) = B(JI)
        B(JI) = T
  40    IJ = IJ + 2
      RETURN
      END                                                               
