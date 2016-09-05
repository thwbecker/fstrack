c
c aniosotropic reflectivity package from Vera Schulte-Pelkum in March 2005
c based on code by Peter Shearer and various others, see individual READMEs 
c and source code comments for details and description
c  
c minor modifications by Thorsten Becker
c
c $Id: sgeev.f,v 1.2 2005/03/23 20:29:35 becker Exp $
c
c
cc This contains SGEEV and all necessary subroutines.  These subroutines are:
c    balanc,balbak,orthes,ortran,hqr,hqr2,scopy,scopym,xerror,cdiv,xerrwv,
c    fdump,i1mach,j4save,xerabt,xerctl,xerprt,xersav,xgetua,s88fmt
c
c
c
*DECK SGEEV
      SUBROUTINE SGEEV(A,LDA,N,E,V,LDV,WORK,JOB,INFO)
C***BEGIN PROLOGUE  SGEEV
C***REVISION DATE  811015   (YYMMDD)
C***CATEGORY NO.  F2C3
C***KEYWORDS  EIGENVALUE,EIGENVECTOR,REAL,GENERAL MATRIX
C***AUTHOR  KAHANER D.K.(NBS), MOLER C.B.(UNM), STEWART G.W.(UMD)
C***DATE WRITTEN  AUGUST 8, 1980
C***PURPOSE
C   TO COMPUTE THE EIGENVALUES AND, OPTIONALLY, THE EIGENVECTORS
C   OF A GENERAL REAL MATRIX.
C***DESCRIPTION
C
C     LICEPACK.    THIS VERSION DATED 08/08/80.
C     DAVID KAHANER, CLEVE MOLER, G. W. STEWART,
C       N.B.S.         U.N.M.      N.B.S./U.MD.
C
C     ABSTRACT
C      SGEEV COMPUTES THE EIGENVALUES AND, OPTIONALLY,
C      THE EIGENVECTORS OF A GENERAL REAL MATRIX.
C
C     CALL SEQUENCE PARAMETERS-
C       (THE VALUES OF PARAMETERS MARKED WITH * (STAR) WILL BE CHANGED
C         BY SGEEV.)
C
C        A*      REAL(LDA,N)
C                REAL NONSYMMETRIC INPUT MATRIX.
C
C        LDA     INTEGER
C                SET BY THE USER TO
C                THE LEADING DIMENSION OF THE REAL ARRAY A.
C
C        N       INTEGER
C                SET BY THE USER TO
C                THE ORDER OF THE MATRICES A AND V, AND
C                THE NUMBER OF ELEMENTS IN E.
C
C        E*      COMPLEX(N)
C                ON RETURN FROM SGEEV E CONTAINS THE EIGENVALUES OF A.
C                SEE ALSO INFO BELOW.
C
C        V*      COMPLEX(LDV,N)
C                ON RETURN FROM SGEEV IF THE USER HAS SET JOB
C                = 0        V IS NOT REFERENCED.
C                = NONZERO  THE N EIGENVECTORS OF A ARE STORED IN THE
C                FIRST N COLUMNS OF V.  SEE ALSO INFO BELOW.
C                (NOTE THAT IF THE INPUT MATRIX A IS NEARLY DEGENERATE,
C                 V MAY BE BADLY CONDITIONED, I.E. MAY HAVE NEARLY
C                 DEPENDENT COLUMNS.)
C
C        LDV     INTEGER
C                SET BY THE USER TO
C                THE LEADING DIMENSION OF THE ARRAY V IF JOB IS ALSO
C                SET NONZERO.  IN THAT CASE N MUST BE .LE. LDV.
C                IF JOB IS SET TO ZERO LDV IS NOT REFERENCED.
C
C        WORK*   REAL(2N)
C                TEMPORARY STORAGE VECTOR.  CONTENTS CHANGED BY SGEEV.
C
C        JOB     INTEGER
C                SET BY THE USER TO
C                = 0        EIGENVALUES ONLY TO BE CALCULATED BY SGEEV.
C                           NEITHER V NOR LDV ARE REFERENCED.
C                = NONZERO  EIGENVALUES AND VECTORS TO BE CALCULATED.
C                           IN THIS CASE A AND V MUST BE DISTINCT ARRAYS
C                           ALSO IF LDA.GT.LDV SGEEV CHANGES ALL THE
C                           ELEMENTS OF A THRU COLUMN N.  IF LDA.LT.LDV
C                           SGEEV CHANGES ALL THE ELEMENTS OF V THROUGH
C                           COLUMN N. IF LDA.EQ.LDV ONLY A(I,J) AND V(I,
C                           J) FOR I,J = 1,...,N ARE CHANGED BY SGEEV.
C
C        INFO*   INTEGER
C                ON RETURN FROM SGEEV THE VALUE OF INFO IS
C                = 0  NORMAL RETURN, CALCULATION SUCCESSFUL.
C                = K  IF THE EIGENVALUE ITERATION FAILS TO CONVERGE,
C                     EIGENVALUES K+1 THROUGH N ARE CORRECT, BUT
C                     NO EIGENVECTORS WERE COMPUTED EVEN IF THEY WERE
C                     REQUESTED (JOB NONZERO).
C
C      ERROR MESSAGES
C           NO. 1  RECOVERABLE  N IS GREATER THAN LDA
C           NO. 2  RECOVERABLE  N IS LESS THAN ONE.
C           NO. 3  RECOVERABLE  JOB IS NONZERO AND N IS GREATER THAN LDV
C           NO. 4  WARNING      LDA.GT.LDV,  ELEMENTS OF A OTHER THAN TH
C                               N BY N INPUT ELEMENTS HAVE BEEN CHANGED.
C           NO. 5  WARNING      LDA.LT.LDV,  ELEMENTS OF V OTHER THAN TH
C                               N BY N OUTPUT ELEMENTS HAVE BEEN CHANGED
C
C
C     SUBROUTINES USED
C
C     EISPACK-  BALANC,BALBAK, ORTHES, ORTRAN, HQR, HQR2
C     BLAS-  SCOPY, SCOPYM
C     SLATEC- XERROR
C***REFERENCES
C***ROUTINES CALLED BALANC,BALBAK,ORTHES,ORTRAN,HQR,HQR2,SCOPY,
C  SCOPYM,XERROR
C***END PROLOGUE  SGEEV
      implicit double precision (a-h,o-z)
      INTEGER I,IHI,ILO,INFO,J,JB,JOB,K,KM,KP,L,LDA,LDV,
     1        MDIM,MIN0,N
      REAL*8 A(1),E(1),WORK(1),V(1)
C***FIRST EXECUTABLE STATEMENT  SGEEV
      IF(N .GT. LDA) CALL XERROR(17HSGEEV-N .GT. LDA.,17,1,1)
      IF (N .GT. LDA) RETURN
      IF(N .LT. 1) CALL XERROR(14HSGEEV-N .LT. 1,14,2,1)
      IF(N .LT. 1) RETURN
      IF(N .EQ. 1 .AND. JOB .EQ. 0) GO TO 35
      MDIM = LDA
      IF(JOB .EQ. 0) GO TO 5
      IF(N .GT. LDV) CALL XERROR(33HSGEEV-JOB .NE. 0, AND N .GT. LDV.,
     1  33,3,1)
      IF(N .GT. LDV) RETURN
      IF(N .EQ. 1) GO TO 35
C
C       REARRANGE A IF NECESSARY WHEN LDA.GT.LDV AND JOB .NE.0
C
      MDIM = MIN0(LDA,LDV)
C***COMMENT BY PETER SHEARER
C***FOR SOME MYSTERIOUS REASON, I COULDN'T GET THESE NEXT THREE LINES
C***TO COMPILE.  THUS I REWROTE THE CODE.
C      IF(LDA.LT.LDV) CALL XERROR
C     1(89HSGEEV-LDA.LT.LDV,  ELEMENTS OF V OTHER THAN THE N BY N OUTPUT
C     2ELEMENTS HAVE BEEN CHANGED.,89,5,0)
      IF (LDA.LT.LDV) THEN
         print *,'error: ELEMENTS OF V OTHER THAN THE...'
         stop
c         CALL XERROR(89HSGEEV-LDA.LT.LDV,  ELEMENTS OF V OTHER THAN THE 
c     &        N BY N OUTPUT ELEMENTS HAVE BEEN CHANGED.,89,5,0)
      END IF
C***END PETER SHEARER MODIFICATION
      IF(LDA.LE.LDV) GO TO 5
c      CALL XERROR(87HSGEEV-LDA.GT.LDV, ELEMENTS OF A OTHER THAN THE N BY 
c     &     N INPUT ELEMENTS HAVE BEEN CHANGED.,87,4,0)

      print *,'error: ELEMENTS OF A OTHER THAN THE N B...'
      stop
      L = N - 1
      DO 4 J=1,L
         M = 1+J*LDV
         K = 1+J*LDA
         CALL SCOPY(N,A(K),1,A(M),1)
    4 CONTINUE
    5 CONTINUE
C
C     SCALE AND ORTHOGONAL REDUCTION TO HESSENBERG.
C
      CALL BALANC(MDIM,N,A,ILO,IHI,WORK(1))
      CALL ORTHES(MDIM,N,ILO,IHI,A,WORK(N+1))
      IF(JOB .NE. 0) GO TO 10
C
C     EIGENVALUES ONLY
C
      CALL HQR(LDA,N,ILO,IHI,A,E(1),E(N+1),INFO)
      GO TO 30
C
C     EIGENVALUES AND EIGENVECTORS.
C
   10 CALL ORTRAN(MDIM,N,ILO,IHI,A,WORK(N+1),V)
      CALL HQR2(MDIM,N,ILO,IHI,A,E(1),E(N+1),V,INFO)
      IF (INFO .NE. 0) GO TO 30
      CALL BALBAK(MDIM,N,ILO,IHI,WORK(1),N,V)
C
C     CONVERT EIGENVECTORS TO COMPLEX STORAGE.
C
      DO 20 JB = 1,N
         J=N+1-JB
         I=N+J
         K=(J-1)*MDIM+1
         KP=K+MDIM
         KM=K-MDIM
         IF(E(I).GE.0.0D0) CALL SCOPY(N,V(K),1,WORK(1),2)
         IF(E(I).LT.0.0D0) CALL SCOPY(N,V(KM),1,WORK(1),2)
         IF(E(I).EQ.0.0D0) CALL SCOPY(N,0.0D0,0,WORK(2),2)
         IF(E(I).GT.0.0D0) CALL SCOPY(N,V(KP),1,WORK(2),2)
         IF(E(I).LT.0.0D0) CALL SCOPYM(N,V(K),1,WORK(2),2)
         L=2*(J-1)*LDV+1
         CALL SCOPY(2*N,WORK(1),1,V(L),1)
   20 CONTINUE
C
C     CONVERT EIGENVALUES TO COMPLEX STORAGE.
C
   30 CALL SCOPY(N,E(1),1,WORK(1),1)
      CALL SCOPY(N,E(N+1),1,E(2),2)
      CALL SCOPY(N,WORK(1),1,E(1),2)
      RETURN
C
C     TAKE CARE OF N=1 CASE
C
   35 E(1) = A(1)
      E(2) = 0.D0
      INFO = 0
      IF(JOB .EQ. 0) RETURN
      V(1) = A(1)
      V(2) = 0.D0
      RETURN
      END
c
c
c
*DECK BALANC
      SUBROUTINE BALANC(NM,N,A,LOW,IGH,SCALE)
C***BEGIN PROLOGUE  BALANC
C***REFER TO  EISDOC
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE BALANCE,
C     NUM. MATH. 13, 293-304(1969) BY PARLETT AND REINSCH.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 315-326(1971).
C
C     THIS SUBROUTINE BALANCES A REAL MATRIX AND ISOLATES
C     EIGENVALUES WHENEVER POSSIBLE.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        A CONTAINS THE INPUT MATRIX TO BE BALANCED.
C
C     ON OUTPUT
C
C        A CONTAINS THE BALANCED MATRIX.
C
C        LOW AND IGH ARE TWO INTEGERS SUCH THAT A(I,J)
C          IS EQUAL TO ZERO IF
C           (1) I IS GREATER THAN J AND
C           (2) J=1,...,LOW-1 OR I=IGH+1,...,N.
C
C        SCALE CONTAINS INFORMATION DETERMINING THE
C           PERMUTATIONS AND SCALING FACTORS USED.
C
C     SUPPOSE THAT THE PRINCIPAL SUBMATRIX IN ROWS LOW THROUGH IGH
C     HAS BEEN BALANCED, THAT P(J) DENOTES THE INDEX INTERCHANGED
C     WITH J DURING THE PERMUTATION STEP, AND THAT THE ELEMENTS
C     OF THE DIAGONAL MATRIX USED ARE DENOTED BY D(I,J).  THEN
C        SCALE(J) = P(J),    FOR J = 1,...,LOW-1
C                 = D(J,J),      J = LOW,...,IGH
C                 = P(J)         J = IGH+1,...,N.
C     THE ORDER IN WHICH THE INTERCHANGES ARE MADE IS N TO IGH+1,
C     THEN 1 TO LOW-1.
C
C     NOTE THAT 1 IS RETURNED FOR IGH IF IGH IS ZERO FORMALLY.
C
C     THE ALGOL PROCEDURE EXC CONTAINED IN BALANCE APPEARS IN
C     BALANC  IN LINE.  (NOTE THAT THE ALGOL ROLES OF IDENTIFIERS
C     K,L HAVE BEEN REVERSED.)
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  BALANC
C
      implicit double precision (a-h,o-z)
      INTEGER I,J,K,L,M,N,JJ,NM,IGH,LOW,IEXC
      REAL*8 A(NM,N),SCALE(N)
      REAL*8 C,F,G,R,S,B2,RADIX
      LOGICAL NOCONV
C
C***FIRST EXECUTABLE STATEMENT  BALANC
      RADIX = 16
C
      B2 = RADIX * RADIX
      K = 1
      L = N
      GO TO 100
C     .......... IN-LINE PROCEDURE FOR ROW AND
C                COLUMN EXCHANGE ..........
   20 SCALE(M) = J
      IF (J .EQ. M) GO TO 50
C
      DO 30 I = 1, L
         F = A(I,J)
         A(I,J) = A(I,M)
         A(I,M) = F
   30 CONTINUE
C
      DO 40 I = K, N
         F = A(J,I)
         A(J,I) = A(M,I)
         A(M,I) = F
   40 CONTINUE
C
   50 GO TO (80,130), IEXC
C     .......... SEARCH FOR ROWS ISOLATING AN EIGENVALUE
C                AND PUSH THEM DOWN ..........
   80 IF (L .EQ. 1) GO TO 280
      L = L - 1
C     .......... FOR J=L STEP -1 UNTIL 1 DO -- ..........
  100 DO 120 JJ = 1, L
         J = L + 1 - JJ
C
         DO 110 I = 1, L
            IF (I .EQ. J) GO TO 110
            IF (A(J,I) .NE. 0.0D0) GO TO 120
  110    CONTINUE
C
         M = L
         IEXC = 1
         GO TO 20
  120 CONTINUE
C
      GO TO 140
C     .......... SEARCH FOR COLUMNS ISOLATING AN EIGENVALUE
C                AND PUSH THEM LEFT ..........
  130 K = K + 1
C
  140 DO 170 J = K, L
C
         DO 150 I = K, L
            IF (I .EQ. J) GO TO 150
            IF (A(I,J) .NE. 0.0D0) GO TO 170
  150    CONTINUE
C
         M = K
         IEXC = 2
         GO TO 20
  170 CONTINUE
C     .......... NOW BALANCE THE SUBMATRIX IN ROWS K TO L ..........
      DO 180 I = K, L
  180 SCALE(I) = 1.0D0
C     .......... ITERATIVE LOOP FOR NORM REDUCTION ..........
  190 NOCONV = .FALSE.
C
      DO 270 I = K, L
         C = 0.0D0
         R = 0.0D0
C
         DO 200 J = K, L
            IF (J .EQ. I) GO TO 200
            C = C + ABS(A(J,I))
            R = R + ABS(A(I,J))
  200    CONTINUE
C     .......... GUARD AGAINST ZERO C OR R DUE TO UNDERFLOW ..........
         IF (C .EQ. 0.0D0 .OR. R .EQ. 0.0D0) GO TO 270
         G = R / RADIX
         F = 1.0D0
         S = C + R
  210    IF (C .GE. G) GO TO 220
         F = F * RADIX
         C = C * B2
         GO TO 210
  220    G = R * RADIX
  230    IF (C .LT. G) GO TO 240
         F = F / RADIX
         C = C / B2
         GO TO 230
C     .......... NOW BALANCE ..........
  240    IF ((C + R) / F .GE. 0.95D0 * S) GO TO 270
         G = 1.0D0 / F
         SCALE(I) = SCALE(I) * F
         NOCONV = .TRUE.
C
         DO 250 J = K, N
  250    A(I,J) = A(I,J) * G
C
         DO 260 J = 1, L
  260    A(J,I) = A(J,I) * F
C
  270 CONTINUE
C
      IF (NOCONV) GO TO 190
C
  280 LOW = K
      IGH = L
      RETURN
      END
c
c
c
*DECK BALBAK
      SUBROUTINE BALBAK(NM,N,LOW,IGH,SCALE,M,Z)
C***BEGIN PROLOGUE  BALBAK
C***REFER TO  EISDOC
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE BALBAK,
C     NUM. MATH. 13, 293-304(1969) BY PARLETT AND REINSCH.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 315-326(1971).
C
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A REAL GENERAL
C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
C     BALANCED MATRIX DETERMINED BY  BALANC.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY  BALANC.
C
C        SCALE CONTAINS INFORMATION DETERMINING THE PERMUTATIONS
C          AND SCALING FACTORS USED BY  BALANC.
C
C        M IS THE NUMBER OF COLUMNS OF Z TO BE BACK TRANSFORMED.
C
C        Z CONTAINS THE REAL AND IMAGINARY PARTS OF THE EIGEN-
C          VECTORS TO BE BACK TRANSFORMED IN ITS FIRST M COLUMNS.
C
C     ON OUTPUT
C
C        Z CONTAINS THE REAL AND IMAGINARY PARTS OF THE
C          TRANSFORMED EIGENVECTORS IN ITS FIRST M COLUMNS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  BALBAK
C
      implicit double precision (a-h,o-z)
      INTEGER I,J,K,M,N,II,NM,IGH,LOW
      REAL*8 SCALE(N),Z(NM,M)
      REAL*8 S
C
C***FIRST EXECUTABLE STATEMENT  BALBAK
      IF (M .EQ. 0) GO TO 200
      IF (IGH .EQ. LOW) GO TO 120
C
      DO 110 I = LOW, IGH
         S = SCALE(I)
C     .......... LEFT HAND EIGENVECTORS ARE BACK TRANSFORMED
C                IF THE FOREGOING STATEMENT IS REPLACED BY
C                S=1.0D0/SCALE(I). ..........
         DO 100 J = 1, M
  100    Z(I,J) = Z(I,J) * S
C
  110 CONTINUE
C     ......... FOR I=LOW-1 STEP -1 UNTIL 1,
C               IGH+1 STEP 1 UNTIL N DO -- ..........
  120 DO 140 II = 1, N
         I = II
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 140
         IF (I .LT. LOW) I = LOW - II
         K = SCALE(I)
         IF (K .EQ. I) GO TO 140
C
         DO 130 J = 1, M
            S = Z(I,J)
            Z(I,J) = Z(K,J)
            Z(K,J) = S
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
      END
c
c
c
*DECK ORTHES
      SUBROUTINE ORTHES(NM,N,LOW,IGH,A,ORT)
C***BEGIN PROLOGUE  ORTHES
C***REFER TO  EISDOC
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE ORTHES,
C     NUM. MATH. 12, 349-368(1968) BY MARTIN AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
C
C     GIVEN A REAL GENERAL MATRIX, THIS SUBROUTINE
C     REDUCES A SUBMATRIX SITUATED IN ROWS AND COLUMNS
C     LOW THROUGH IGH TO UPPER HESSENBERG FORM BY
C     ORTHOGONAL SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        A CONTAINS THE INPUT MATRIX.
C
C     ON OUTPUT
C
C        A CONTAINS THE HESSENBERG MATRIX.  INFORMATION ABOUT
C          THE ORTHOGONAL TRANSFORMATIONS USED IN THE REDUCTION
C          IS STORED IN THE REMAINING TRIANGLE UNDER THE
C          HESSENBERG MATRIX.
C
C        ORT CONTAINS FURTHER INFORMATION ABOUT THE TRANSFORMATIONS.
C          ONLY ELEMENTS LOW THROUGH IGH ARE USED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  ORTHES
C
      implicit double precision (a-h,o-z)
      INTEGER I,J,M,N,II,JJ,LA,MP,NM,IGH,KP1,LOW
      REAL*8 A(NM,N),ORT(IGH)
      REAL*8 F,G,H,SCALE
C
C***FIRST EXECUTABLE STATEMENT  ORTHES
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA .LT. KP1) GO TO 200
C
      DO 180 M = KP1, LA
         H = 0.0D0
         ORT(M) = 0.0D0
         SCALE = 0.0D0
C     .......... SCALE COLUMN (ALGOL TOL THEN NOT NEEDED) ..........
         DO 90 I = M, IGH
   90    SCALE = SCALE + ABS(A(I,M-1))
C
         IF (SCALE .EQ. 0.0D0) GO TO 180
         MP = M + IGH
C     .......... FOR I=IGH STEP -1 UNTIL M DO -- ..........
         DO 100 II = M, IGH
            I = MP - II
            ORT(I) = A(I,M-1) / SCALE
            H = H + ORT(I) * ORT(I)
  100    CONTINUE
C
         G = -SIGN(SQRT(H),ORT(M))
         H = H - ORT(M) * G
         ORT(M) = ORT(M) - G
C     .......... FORM (I-(U*UT)/H) * A ..........
         DO 130 J = M, N
            F = 0.0D0
C     .......... FOR I=IGH STEP -1 UNTIL M DO -- ..........
            DO 110 II = M, IGH
               I = MP - II
               F = F + ORT(I) * A(I,J)
  110       CONTINUE
C
            F = F / H
C
            DO 120 I = M, IGH
  120       A(I,J) = A(I,J) - F * ORT(I)
C
  130    CONTINUE
C     .......... FORM (I-(U*UT)/H)*A*(I-(U*UT)/H) ..........
         DO 160 I = 1, IGH
            F = 0.0D0
C     .......... FOR J=IGH STEP -1 UNTIL M DO -- ..........
            DO 140 JJ = M, IGH
               J = MP - JJ
               F = F + ORT(J) * A(I,J)
  140       CONTINUE
C
            F = F / H
C
            DO 150 J = M, IGH
  150       A(I,J) = A(I,J) - F * ORT(J)
C
  160    CONTINUE
C
         ORT(M) = SCALE * ORT(M)
         A(M,M-1) = SCALE * G
  180 CONTINUE
C
  200 RETURN
      END
c
c
c
*DECK ORTRAN
      SUBROUTINE ORTRAN(NM,N,LOW,IGH,A,ORT,Z)
C***BEGIN PROLOGUE  ORTRAN
C***REFER TO  EISDOC
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE ORTRANS,
C     NUM. MATH. 16, 181-204(1970) BY PETERS AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
C
C     THIS SUBROUTINE ACCUMULATES THE ORTHOGONAL SIMILARITY
C     TRANSFORMATIONS USED IN THE REDUCTION OF A REAL GENERAL
C     MATRIX TO UPPER HESSENBERG FORM BY  ORTHES.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANS-
C          FORMATIONS USED IN THE REDUCTION BY  ORTHES
C          IN ITS STRICT LOWER TRIANGLE.
C
C        ORT CONTAINS FURTHER INFORMATION ABOUT THE TRANS-
C          FORMATIONS USED IN THE REDUCTION BY  ORTHES.
C          ONLY ELEMENTS LOW THROUGH IGH ARE USED.
C
C     ON OUTPUT
C
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C          REDUCTION BY  ORTHES.
C
C        ORT HAS BEEN ALTERED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     .......... INITIALIZE Z TO IDENTITY MATRIX ..........
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  ORTRAN
C
      implicit double precision (a-h,o-z)
      INTEGER I,J,N,KL,MM,MP,NM,IGH,LOW,MP1
      REAL*8 A(NM,IGH),ORT(IGH),Z(NM,N)
      REAL*8 G
C
C***FIRST EXECUTABLE STATEMENT  ORTRAN
C     .......... INITIALIZE Z TO IDENTITY MATRIX ..........
      DO 80 I = 1, N
C
         DO 60 J = 1, N
   60    Z(I,J) = 0.0D0
C
         Z(I,I) = 1.0D0
   80 CONTINUE
C
      KL = IGH - LOW - 1
      IF (KL .LT. 1) GO TO 200
C     .......... FOR MP=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
      DO 140 MM = 1, KL
         MP = IGH - MM
         IF (A(MP,MP-1) .EQ. 0.0D0) GO TO 140
         MP1 = MP + 1
C
         DO 100 I = MP1, IGH
  100    ORT(I) = A(I,MP-1)
C
         DO 130 J = MP, IGH
            G = 0.0D0
C
            DO 110 I = MP, IGH
  110       G = G + ORT(I) * Z(I,J)
C     .......... DIVISOR BELOW IS NEGATIVE OF H FORMED IN ORTHES.
C                DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
            G = (G / ORT(MP)) / A(MP,MP-1)
C
            DO 120 I = MP, IGH
  120       Z(I,J) = Z(I,J) + G * ORT(I)
C
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
      END
c
c
c
*DECK HQR
      SUBROUTINE HQR(NM,N,LOW,IGH,H,WR,WI,IERR)
C***BEGIN PROLOGUE  HQR
C***REFER TO  EISDOC
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE HQR,
C     NUM. MATH. 14, 219-231(1970) BY MARTIN, PETERS, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 359-371(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A REAL
C     UPPER HESSENBERG MATRIX BY THE QR METHOD.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        H CONTAINS THE UPPER HESSENBERG MATRIX.  INFORMATION ABOUT
C          THE TRANSFORMATIONS USED IN THE REDUCTION TO HESSENBERG
C          FORM BY  ELMHES  OR  ORTHES, IF PERFORMED, IS STORED
C          IN THE REMAINING TRIANGLE UNDER THE HESSENBERG MATRIX.
C
C     ON OUTPUT
C
C        H HAS BEEN DESTROYED.  THEREFORE, IT MUST BE SAVED
C          BEFORE CALLING  HQR  IF SUBSEQUENT CALCULATION AND
C          BACK TRANSFORMATION OF EIGENVECTORS IS TO BE PERFORMED.
C
C        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVALUES.  THE EIGENVALUES
C          ARE UNORDERED EXCEPT THAT COMPLEX CONJUGATE PAIRS
C          OF VALUES APPEAR CONSECUTIVELY WITH THE EIGENVALUE
C          HAVING THE POSITIVE IMAGINARY PART FIRST.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES SHOULD BE CORRECT
C          FOR INDICES IERR+1,...,N.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER A TOTAL OF 30*N ITERATIONS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  HQR
C
      implicit double precision (a-h,o-z)
      INTEGER I,J,K,L,M,N,EN,LL,MM,NA,NM,IGH,ITN,ITS,LOW,MP2,ENM2,IERR
      REAL*8 H(NM,N),WR(N),WI(N)
      REAL*8 P,Q,R,S,T,W,X,Y,ZZ,NORM,S1,S2
      LOGICAL NOTLAS
C
C***FIRST EXECUTABLE STATEMENT  HQR
      IERR = 0
      NORM = 0.0D0
      K = 1
C     .......... STORE ROOTS ISOLATED BY BALANC
C                AND COMPUTE MATRIX NORM ..........
      DO 50 I = 1, N
C
         DO 40 J = K, N
   40    NORM = NORM + ABS(H(I,J))
C
         K = I
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 50
         WR(I) = H(I,I)
         WI(I) = 0.0D0
   50 CONTINUE
C
      EN = IGH
      T = 0.0D0
      ITN = 30*N
C     .......... SEARCH FOR NEXT EIGENVALUES ..........
   60 IF (EN .LT. LOW) GO TO 1001
      ITS = 0
      NA = EN - 1
      ENM2 = NA - 1
C     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
C                FOR L=EN STEP -1 UNTIL LOW DO -- ..........
   70 DO 80 LL = LOW, EN
         L = EN + LOW - LL
         IF (L .EQ. LOW) GO TO 100
         S = ABS(H(L-1,L-1)) + ABS(H(L,L))
         IF (S .EQ. 0.0D0) S = NORM
         S2 = S + ABS(H(L,L-1))
         IF (S2 .EQ. S) GO TO 100
   80 CONTINUE
C     .......... FORM SHIFT ..........
  100 X = H(EN,EN)
      IF (L .EQ. EN) GO TO 270
      Y = H(NA,NA)
      W = H(EN,NA) * H(NA,EN)
      IF (L .EQ. NA) GO TO 280
      IF (ITN .EQ. 0) GO TO 1000
      IF (ITS .NE. 10 .AND. ITS .NE. 20) GO TO 130
C     .......... FORM EXCEPTIONAL SHIFT ..........
      T = T + X
C
      DO 120 I = LOW, EN
  120 H(I,I) = H(I,I) - X
C
      S = ABS(H(EN,NA)) + ABS(H(NA,ENM2))
      X = 0.75D0 * S
      Y = X
      W = -0.4375D0 * S * S
  130 ITS = ITS + 1
      ITN = ITN - 1
C     .......... LOOK FOR TWO CONSECUTIVE SMALL
C                SUB-DIAGONAL ELEMENTS.
C                FOR M=EN-2 STEP -1 UNTIL L DO -- ..........
      DO 140 MM = L, ENM2
         M = ENM2 + L - MM
         ZZ = H(M,M)
         R = X - ZZ
         S = Y - ZZ
         P = (R * S - W) / H(M+1,M) + H(M,M+1)
         Q = H(M+1,M+1) - ZZ - R - S
         R = H(M+2,M+1)
         S = ABS(P) + ABS(Q) + ABS(R)
         P = P / S
         Q = Q / S
         R = R / S
         IF (M .EQ. L) GO TO 150
         S1 = ABS(P) * (ABS(H(M-1,M-1)) + ABS(ZZ) + ABS(H(M+1,M+1)))
         S2 = S1 + ABS(H(M,M-1)) * (ABS(Q) + ABS(R))
         IF (S2 .EQ. S1) GO TO 150
  140 CONTINUE
C
  150 MP2 = M + 2
C
      DO 160 I = MP2, EN
         H(I,I-2) = 0.0D0
         IF (I .EQ. MP2) GO TO 160
         H(I,I-3) = 0.0D0
  160 CONTINUE
C     .......... DOUBLE QR STEP INVOLVING ROWS L TO EN AND
C                COLUMNS M TO EN ..........
      DO 260 K = M, NA
         NOTLAS = K .NE. NA
         IF (K .EQ. M) GO TO 170
         P = H(K,K-1)
         Q = H(K+1,K-1)
         R = 0.0D0
         IF (NOTLAS) R = H(K+2,K-1)
         X = ABS(P) + ABS(Q) + ABS(R)
         IF (X .EQ. 0.0D0) GO TO 260
         P = P / X
         Q = Q / X
         R = R / X
  170    S = SIGN(SQRT(P*P+Q*Q+R*R),P)
         IF (K .EQ. M) GO TO 180
         H(K,K-1) = -S * X
         GO TO 190
  180    IF (L .NE. M) H(K,K-1) = -H(K,K-1)
  190    P = P + S
         X = P / S
         Y = Q / S
         ZZ = R / S
         Q = Q / P
         R = R / P
C     .......... ROW MODIFICATION ..........
         DO 210 J = K, EN
            P = H(K,J) + Q * H(K+1,J)
            IF (.NOT. NOTLAS) GO TO 200
            P = P + R * H(K+2,J)
            H(K+2,J) = H(K+2,J) - P * ZZ
  200       H(K+1,J) = H(K+1,J) - P * Y
            H(K,J) = H(K,J) - P * X
  210    CONTINUE
C
         J = MIN0(EN,K+3)
C     .......... COLUMN MODIFICATION ..........
         DO 230 I = L, J
            P = X * H(I,K) + Y * H(I,K+1)
            IF (.NOT. NOTLAS) GO TO 220
            P = P + ZZ * H(I,K+2)
            H(I,K+2) = H(I,K+2) - P * R
  220       H(I,K+1) = H(I,K+1) - P * Q
            H(I,K) = H(I,K) - P
  230    CONTINUE
C
  260 CONTINUE
C
      GO TO 70
C     .......... ONE ROOT FOUND ..........
  270 WR(EN) = X + T
      WI(EN) = 0.0D0
      EN = NA
      GO TO 60
C     .......... TWO ROOTS FOUND ..........
  280 P = (Y - X) / 2.0D0
      Q = P * P + W
      ZZ = SQRT(ABS(Q))
      X = X + T
      IF (Q .LT. 0.0D0) GO TO 320
C     .......... REAL PAIR ..........
      ZZ = P + SIGN(ZZ,P)
      WR(NA) = X + ZZ
      WR(EN) = WR(NA)
      IF (ZZ .NE. 0.0D0) WR(EN) = X - W / ZZ
      WI(NA) = 0.0D0
      WI(EN) = 0.0D0
      GO TO 330
C     .......... COMPLEX PAIR ..........
  320 WR(NA) = X + P
      WR(EN) = X + P
      WI(NA) = ZZ
      WI(EN) = -ZZ
  330 EN = ENM2
      GO TO 60
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30*N ITERATIONS ..........
 1000 IERR = EN
 1001 RETURN
      END
c
c
c
*DECK HQR2
      SUBROUTINE HQR2(NM,N,LOW,IGH,H,WR,WI,Z,IERR)
C***BEGIN PROLOGUE  HQR2
C***REFER TO  EISDOC
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE HQR2,
C     NUM. MATH. 16, 181-204(1970) BY PETERS AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
C     OF A REAL UPPER HESSENBERG MATRIX BY THE QR METHOD.  THE
C     EIGENVECTORS OF A REAL GENERAL MATRIX CAN ALSO BE FOUND
C     IF  ELMHES  AND  ELTRAN  OR  ORTHES  AND  ORTRAN  HAVE
C     BEEN USED TO REDUCE THIS GENERAL MATRIX TO HESSENBERG FORM
C     AND TO ACCUMULATE THE SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        H CONTAINS THE UPPER HESSENBERG MATRIX.
C
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED BY  ELTRAN
C          AFTER THE REDUCTION BY  ELMHES, OR BY  ORTRAN  AFTER THE
C          REDUCTION BY  ORTHES, IF PERFORMED.  IF THE EIGENVECTORS
C          OF THE HESSENBERG MATRIX ARE DESIRED, Z MUST CONTAIN THE
C          IDENTITY MATRIX.
C
C     ON OUTPUT
C
C        H HAS BEEN DESTROYED.
C
C        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVALUES.  THE EIGENVALUES
C          ARE UNORDERED EXCEPT THAT COMPLEX CONJUGATE PAIRS
C          OF VALUES APPEAR CONSECUTIVELY WITH THE EIGENVALUE
C          HAVING THE POSITIVE IMAGINARY PART FIRST.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES SHOULD BE CORRECT
C          FOR INDICES IERR+1,...,N.
C
C        Z CONTAINS THE REAL AND IMAGINARY PARTS OF THE EIGENVECTORS.
C          IF THE I-TH EIGENVALUE IS REAL, THE I-TH COLUMN OF Z
C          CONTAINS ITS EIGENVECTOR.  IF THE I-TH EIGENVALUE IS COMPLEX
C          WITH POSITIVE IMAGINARY PART, THE I-TH AND (I+1)-TH
C          COLUMNS OF Z CONTAIN THE REAL AND IMAGINARY PARTS OF ITS
C          EIGENVECTOR.  THE EIGENVECTORS ARE UNNORMALIZED.  IF AN
C          ERROR EXIT IS MADE, NONE OF THE EIGENVECTORS HAS BEEN FOUND.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER A TOTAL OF 30*N ITERATIONS.
C
C     CALLS CDIV FOR COMPLEX DIVISION.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C***ROUTINES CALLED  CDIV
C***END PROLOGUE  HQR2
C
      implicit double precision (a-h,o-z)
      INTEGER I,J,K,L,M,N,EN,II,JJ,LL,MM,NA,NM,NN
      INTEGER IGH,ITN,ITS,LOW,MP2,ENM2,IERR
      REAL*8 H(NM,N),WR(N),WI(N),Z(NM,N)
      REAL*8 P,Q,R,S,T,W,X,Y,RA,SA,VI,VR,ZZ,NORM,S1
      LOGICAL NOTLAS
C
C***FIRST EXECUTABLE STATEMENT  HQR2
      IERR = 0
      NORM = 0.0D0
      K = 1
C     .......... STORE ROOTS ISOLATED BY BALANC
C                AND COMPUTE MATRIX NORM ..........
      DO 50 I = 1, N
C
         DO 40 J = K, N
   40    NORM = NORM + ABS(H(I,J))
C
         K = I
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 50
         WR(I) = H(I,I)
         WI(I) = 0.0D0
   50 CONTINUE
C
      EN = IGH
      T = 0.0D0
      ITN = 30*N
C     .......... SEARCH FOR NEXT EIGENVALUES ..........
   60 IF (EN .LT. LOW) GO TO 340
      ITS = 0
      NA = EN - 1
      ENM2 = NA - 1
C     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
C                FOR L=EN STEP -1 UNTIL LOW DO -- ..........
   70 DO 80 LL = LOW, EN
         L = EN + LOW - LL
         IF (L .EQ. LOW) GO TO 100
         S = ABS(H(L-1,L-1)) + ABS(H(L,L))
         IF (S .EQ. 0.0D0) S = NORM
         S2 = S + ABS(H(L,L-1))
         IF (S2 .EQ. S) GO TO 100
   80 CONTINUE
C     .......... FORM SHIFT ..........
  100 X = H(EN,EN)
      IF (L .EQ. EN) GO TO 270
      Y = H(NA,NA)
      W = H(EN,NA) * H(NA,EN)
      IF (L .EQ. NA) GO TO 280
      IF (ITN .EQ. 0) GO TO 1000
      IF (ITS .NE. 10 .AND. ITS .NE. 20) GO TO 130
C     .......... FORM EXCEPTIONAL SHIFT ..........
      T = T + X
C
      DO 120 I = LOW, EN
  120 H(I,I) = H(I,I) - X
C
      S = ABS(H(EN,NA)) + ABS(H(NA,ENM2))
      X = 0.75D0 * S
      Y = X
      W = -0.4375D0 * S * S
  130 ITS = ITS + 1
      ITN = ITN - 1
C     .......... LOOK FOR TWO CONSECUTIVE SMALL
C                SUB-DIAGONAL ELEMENTS.
C                FOR M=EN-2 STEP -1 UNTIL L DO -- ..........
      DO 140 MM = L, ENM2
         M = ENM2 + L - MM
         ZZ = H(M,M)
         R = X - ZZ
         S = Y - ZZ
         P = (R * S - W) / H(M+1,M) + H(M,M+1)
         Q = H(M+1,M+1) - ZZ - R - S
         R = H(M+2,M+1)
         S = ABS(P) + ABS(Q) + ABS(R)
         P = P / S
         Q = Q / S
         R = R / S
         IF (M .EQ. L) GO TO 150
         S1 = ABS(P) * (ABS(H(M-1,M-1)) + ABS(ZZ) + ABS(H(M+1,M+1)))
         S2 = S1 + ABS(H(M,M-1)) * (ABS(Q) + ABS(R))
         IF (S2 .EQ. S1) GO TO 150
  140 CONTINUE
C
  150 MP2 = M + 2
C
      DO 160 I = MP2, EN
         H(I,I-2) = 0.0D0
         IF (I .EQ. MP2) GO TO 160
         H(I,I-3) = 0.0D0
  160 CONTINUE
C     .......... DOUBLE QR STEP INVOLVING ROWS L TO EN AND
C                COLUMNS M TO EN ..........
      DO 260 K = M, NA
         NOTLAS = K .NE. NA
         IF (K .EQ. M) GO TO 170
         P = H(K,K-1)
         Q = H(K+1,K-1)
         R = 0.0D0
         IF (NOTLAS) R = H(K+2,K-1)
         X = ABS(P) + ABS(Q) + ABS(R)
         IF (X .EQ. 0.0D0) GO TO 260
         P = P / X
         Q = Q / X
         R = R / X
  170    S = SIGN(SQRT(P*P+Q*Q+R*R),P)
         IF (K .EQ. M) GO TO 180
         H(K,K-1) = -S * X
         GO TO 190
  180    IF (L .NE. M) H(K,K-1) = -H(K,K-1)
  190    P = P + S
         X = P / S
         Y = Q / S
         ZZ = R / S
         Q = Q / P
         R = R / P
C     .......... ROW MODIFICATION ..........
         DO 210 J = K, N
            P = H(K,J) + Q * H(K+1,J)
            IF (.NOT. NOTLAS) GO TO 200
            P = P + R * H(K+2,J)
            H(K+2,J) = H(K+2,J) - P * ZZ
  200       H(K+1,J) = H(K+1,J) - P * Y
            H(K,J) = H(K,J) - P * X
  210    CONTINUE
C
         J = MIN0(EN,K+3)
C     .......... COLUMN MODIFICATION ..........
         DO 230 I = 1, J
            P = X * H(I,K) + Y * H(I,K+1)
            IF (.NOT. NOTLAS) GO TO 220
            P = P + ZZ * H(I,K+2)
            H(I,K+2) = H(I,K+2) - P * R
  220       H(I,K+1) = H(I,K+1) - P * Q
            H(I,K) = H(I,K) - P
  230    CONTINUE
C     .......... ACCUMULATE TRANSFORMATIONS ..........
         DO 250 I = LOW, IGH
            P = X * Z(I,K) + Y * Z(I,K+1)
            IF (.NOT. NOTLAS) GO TO 240
            P = P + ZZ * Z(I,K+2)
            Z(I,K+2) = Z(I,K+2) - P * R
  240       Z(I,K+1) = Z(I,K+1) - P * Q
            Z(I,K) = Z(I,K) - P
  250    CONTINUE
C
  260 CONTINUE
C
      GO TO 70
C     .......... ONE ROOT FOUND ..........
  270 H(EN,EN) = X + T
      WR(EN) = H(EN,EN)
      WI(EN) = 0.0D0
      EN = NA
      GO TO 60
C     .......... TWO ROOTS FOUND ..........
  280 P = (Y - X) / 2.0D0
      Q = P * P + W
      ZZ = SQRT(ABS(Q))
      H(EN,EN) = X + T
      X = H(EN,EN)
      H(NA,NA) = Y + T
      IF (Q .LT. 0.0D0) GO TO 320
C     .......... REAL PAIR ..........
      ZZ = P + SIGN(ZZ,P)
      WR(NA) = X + ZZ
      WR(EN) = WR(NA)
      IF (ZZ .NE. 0.0D0) WR(EN) = X - W / ZZ
      WI(NA) = 0.0D0
      WI(EN) = 0.0D0
      X = H(EN,NA)
      S = ABS(X) + ABS(ZZ)
      P = X / S
      Q = ZZ / S
      R = SQRT(P*P+Q*Q)
      P = P / R
      Q = Q / R
C     .......... ROW MODIFICATION ..........
      DO 290 J = NA, N
         ZZ = H(NA,J)
         H(NA,J) = Q * ZZ + P * H(EN,J)
         H(EN,J) = Q * H(EN,J) - P * ZZ
  290 CONTINUE
C     .......... COLUMN MODIFICATION ..........
      DO 300 I = 1, EN
         ZZ = H(I,NA)
         H(I,NA) = Q * ZZ + P * H(I,EN)
         H(I,EN) = Q * H(I,EN) - P * ZZ
  300 CONTINUE
C     .......... ACCUMULATE TRANSFORMATIONS ..........
      DO 310 I = LOW, IGH
         ZZ = Z(I,NA)
         Z(I,NA) = Q * ZZ + P * Z(I,EN)
         Z(I,EN) = Q * Z(I,EN) - P * ZZ
  310 CONTINUE
C
      GO TO 330
C     .......... COMPLEX PAIR ..........
  320 WR(NA) = X + P
      WR(EN) = X + P
      WI(NA) = ZZ
      WI(EN) = -ZZ
  330 EN = ENM2
      GO TO 60
C     .......... ALL ROOTS FOUND.  BACKSUBSTITUTE TO FIND
C                VECTORS OF UPPER TRIANGULAR FORM ..........
  340 IF (NORM .EQ. 0.0D0) GO TO 1001
C     .......... FOR EN=N STEP -1 UNTIL 1 DO -- ..........
      DO 800 NN = 1, N
         EN = N + 1 - NN
         P = WR(EN)
         Q = WI(EN)
         NA = EN - 1
         IF (Q) 710, 600, 800
C     .......... REAL VECTOR ..........
  600    M = EN
         H(EN,EN) = 1.0D0
         IF (NA .EQ. 0) GO TO 800
C     .......... FOR I=EN-1 STEP -1 UNTIL 1 DO -- ..........
         DO 700 II = 1, NA
            I = EN - II
            W = H(I,I) - P
            R = H(I,EN)
            IF (M .GT. NA) GO TO 620
C
            DO 610 J = M, NA
  610       R = R + H(I,J) * H(J,EN)
C
  620       IF (WI(I) .GE. 0.0D0) GO TO 630
            ZZ = W
            S = R
            GO TO 700
  630       M = I
            IF (WI(I) .NE. 0.0D0) GO TO 640
            T = W
            IF (T .NE. 0.0D0) GO TO 635
            T = NORM
  632       T = 0.5D0*T
            IF (NORM + T .GT. NORM) GO TO 632
            T = 2.0D0*T
  635       H(I,EN) = -R / T
            GO TO 700
C     .......... SOLVE REAL EQUATIONS ..........
  640       X = H(I,I+1)
            Y = H(I+1,I)
            Q = (WR(I) - P) * (WR(I) - P) + WI(I) * WI(I)
            T = (X * S - ZZ * R) / Q
            H(I,EN) = T
            IF (ABS(X) .LE. ABS(ZZ)) GO TO 650
            H(I+1,EN) = (-R - W * T) / X
            GO TO 700
  650       H(I+1,EN) = (-S - Y * T) / ZZ
  700    CONTINUE
C     .......... END REAL VECTOR ..........
         GO TO 800
C     .......... COMPLEX VECTOR ..........
  710    M = NA
C     .......... LAST VECTOR COMPONENT CHOSEN IMAGINARY SO THAT
C                EIGENVECTOR MATRIX IS TRIANGULAR ..........
         IF (ABS(H(EN,NA)) .LE. ABS(H(NA,EN))) GO TO 720
         H(NA,NA) = Q / H(EN,NA)
         H(NA,EN) = -(H(EN,EN) - P) / H(EN,NA)
         GO TO 730
  720    CALL CDIV(0.0D0,-H(NA,EN),H(NA,NA)-P,Q,H(NA,NA),H(NA,EN))
  730    H(EN,NA) = 0.0D0
         H(EN,EN) = 1.0D0
         ENM2 = NA - 1
         IF (ENM2 .EQ. 0) GO TO 800
C     .......... FOR I=EN-2 STEP -1 UNTIL 1 DO -- ..........
         DO 790 II = 1, ENM2
            I = NA - II
            W = H(I,I) - P
            RA = 0.0D0
            SA = H(I,EN)
C
            DO 760 J = M, NA
               RA = RA + H(I,J) * H(J,NA)
               SA = SA + H(I,J) * H(J,EN)
  760       CONTINUE
C
            IF (WI(I) .GE. 0.0D0) GO TO 770
            ZZ = W
            R = RA
            S = SA
            GO TO 790
  770       M = I
            IF (WI(I) .NE. 0.0D0) GO TO 780
            CALL CDIV(-RA,-SA,W,Q,H(I,NA),H(I,EN))
            GO TO 790
C     .......... SOLVE COMPLEX EQUATIONS ..........
  780       X = H(I,I+1)
            Y = H(I+1,I)
            VR = (WR(I) - P) * (WR(I) - P) + WI(I) * WI(I) - Q * Q
            VI = (WR(I) - P) * 2.0D0 * Q
            IF (VR .NE. 0.0D0 .OR. VI .NE. 0.0D0) GO TO 783
            S1 = NORM * (ABS(W)+ABS(Q)+ABS(X)+ABS(Y)+ABS(ZZ))
            VR = S1
  782       VR = 0.5D0*VR
            IF (S1 + VR .GT. S1) GO TO 782
            VR = 2.0D0*VR
  783       CALL CDIV(X*R-ZZ*RA+Q*SA,X*S-ZZ*SA-Q*RA,VR,VI,
     X                H(I,NA),H(I,EN))
            IF (ABS(X) .LE. ABS(ZZ) + ABS(Q)) GO TO 785
            H(I+1,NA) = (-RA - W * H(I,NA) + Q * H(I,EN)) / X
            H(I+1,EN) = (-SA - W * H(I,EN) - Q * H(I,NA)) / X
            GO TO 790
  785       CALL CDIV(-R-Y*H(I,NA),-S-Y*H(I,EN),ZZ,Q,
     X                H(I+1,NA),H(I+1,EN))
  790    CONTINUE
C     .......... END COMPLEX VECTOR ..........
  800 CONTINUE
C     .......... END BACK SUBSTITUTION.
C                VECTORS OF ISOLATED ROOTS ..........
      DO 840 I = 1, N
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 840
C
         DO 820 J = I, N
  820    Z(I,J) = H(I,J)
C
  840 CONTINUE
C     .......... MULTIPLY BY TRANSFORMATION MATRIX TO GIVE
C                VECTORS OF ORIGINAL FULL MATRIX.
C                FOR J=N STEP -1 UNTIL LOW DO -- ..........
      DO 880 JJ = LOW, N
         J = N + LOW - JJ
         M = MIN0(J,IGH)
C
         DO 880 I = LOW, IGH
            ZZ = 0.0D0
C
            DO 860 K = LOW, M
  860       ZZ = ZZ + Z(I,K) * H(K,J)
C
            Z(I,J) = ZZ
  880 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30*N ITERATIONS ..........
 1000 IERR = EN
 1001 RETURN
      END
c
c
c
*DECK SCOPY
      SUBROUTINE  SCOPY(N,SX,INCX,SY,INCY)
C***BEGIN PROLOGUE  SCOPY
C***REVISION DATE  811015   (YYMMDD)
C***CATEGORY NO.  F1A
C***KEYWORDS  BLAS,VECTOR,COPY
C***DATE WRITTEN  OCTOBER 1979
C***AUTHOR LAWSON C. (JPL),HANSON R. (SLA),
C                            KINCAID D. (U TEXAS), KROGH F. (JPL)
C***PURPOSE
C   COPY S.P. VECTOR Y = X
C***DESCRIPTION
C                B L A S  SUBPROGRAM
C    DESCRIPTION OF PARAMETERS
C
C     --INPUT--
C        N  NUMBER OF ELEMENTS IN INPUT VECTOR(S)
C       SX  SINGLE PRECISION VECTOR WITH N ELEMENTS
C     INCX  STORAGE SPACING BETWEEN ELEMENTS OF SX
C       SY  SINGLE PRECISION VECTOR WITH N ELEMENTS
C     INCY  STORAGE SPACING BETWEEN ELEMENTS OF SY
C
C     --OUTPUT--
C       SY  COPY OF VECTOR SX (UNCHANGED IF N.LE.0)
C
C     COPY SINGLE PRECISION SX TO SINGLE PRECISION SY.
C     FOR I = 0 TO N-1, COPY  SX(LX+I*INCX) TO SY(LY+I*INCY),
C     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS
C     DEFINED IN A SIMILAR WAY USING INCY.
C
C***REFERENCES
C  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C   *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C  ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C  VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  SCOPY
C
      implicit double precision (a-h,o-z)
      REAL*8 SX(1),SY(1)
C***FIRST EXECUTABLE STATEMENT  SCOPY
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
C
C        CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        SY(IY) = SX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 7.
C
   20 M = MOD(N,7)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SY(I) = SX(I)
   30 CONTINUE
      IF( N .LT. 7 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        SY(I) = SX(I)
        SY(I + 1) = SX(I + 1)
        SY(I + 2) = SX(I + 2)
        SY(I + 3) = SX(I + 3)
        SY(I + 4) = SX(I + 4)
        SY(I + 5) = SX(I + 5)
        SY(I + 6) = SX(I + 6)
   50 CONTINUE
      RETURN
C
C        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
C
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          SY(I) = SX(I)
   70     CONTINUE
      RETURN
      END
c
c
c
*DECK SCOPYM
      SUBROUTINE SCOPYM(N,SX,INCX,SY,INCY)
C***BEGIN PROLOGUE  SCOPYM
C***DATE WRITTEN   801001   (YYMMDD)
C***REVISION DATE  811015   (YYMMDD)
C***CATEGORY NO.  F1A
C***KEYWORDS  BLAS,VECTOR,COPY
C***AUTHOR  KAHANER,DAVID(NBS)
C***PURPOSE  COPY NEGATIVE OF REAL SX TO REAL SY.
C***DESCRIPTION
C       DESCRIPTION OF PARAMETERS
C           THE * FLAGS OUTPUT VARIABLES
C
C       N   NUMBER OF ELEMENTS IN VECTOR(S)
C      SX   REAL VECTOR WITH N ELEMENTS
C    INCX   STORAGE SPACING BETWEEN ELEMENTS OF SX
C      SY*  REAL NEGATIVE COPY OF SX
C    INCY   STORAGE SPACING BETWEEN ELEMENTS OF SY
C
C      ***  NOTE THAT SY = -SX  ***
C
C  COPY NEGATIVE OF REAL SX TO REAL SY.  FOR I=0 TO N-1,
C   COPY  -SX(LX+I*INCX) TO SY(LY+I*INCY), WHERE LX=1 IF
C   INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS DEFINED
C   IN A SIMILAR WAY USING INCY.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  SCOPYM
      implicit double precision (a-h,o-z)
      REAL*8 SX(1),SY(1)
C***FIRST EXECUTABLE STATEMENT  SCOPYM
      IF(N.LE.0) RETURN
       IF(INCX.EQ.INCY)  IF(INCX-1) 5,20,60
    5 CONTINUE
C
C         CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS
C
      IX=1
      IY=1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I=1,N
        SY(IY) = -SX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN UP LOOP SO REMAINING VECTOR LENGTH IS MULTIPLE OF 7
C
   20 M = MOD(N,7)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I=1,M
        SY(I) = -SX(I)
   30 CONTINUE
      IF( N .LT. 7 ) RETURN
   40 MP1 = M + 1
      DO 50 I= MP1,N,7
        SY(I) = -SX(I)
        SY(I + 1) = -SX(I + 1)
        SY(I + 2) = -SX(I + 2)
        SY(I + 3) = -SX(I + 3)
        SY(I + 4) = -SX(I + 4)
        SY(I + 5) = -SX(I + 5)
        SY(I + 6) = -SX(I + 6)
   50 CONTINUE
      RETURN
C
C          CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS
C
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          SY(I) = -SX(I)
   70     CONTINUE
      RETURN
      END
c
c
c
*DECK XERROR
      SUBROUTINE XERROR(MESSG,NMESSG,NERR,LEVEL)
C***BEGIN PROLOGUE  XERROR
C***REVISION DATE  811015   (YYMMDD)
C***CATEGORY NO.  Q
C***KEYWORDS  XERROR,ERROR
C***DATE WRITTEN  AUGUST 1979
C***AUTHOR  JONES R.E. (SLA)
C***PURPOSE
C PROCESSES AN ERROR (DIAGNOSTIC) MESSAGE
C***DESCRIPTION
C     ABSTRACT
C        XERROR PROCESSES A DIAGNOSTIC MESSAGE, IN A MANNER
C        DETERMINED BY THE VALUE OF LEVEL AND THE CURRENT VALUE
C        OF THE LIBRARY ERROR CONTROL FLAG, KONTRL.
C        (SEE SUBROUTINE XSETF FOR DETAILS.)
C
C     DESCRIPTION OF PARAMETERS
C      --INPUT--
C        MESSG - THE HOLLERITH MESSAGE TO BE PROCESSED, CONTAINING
C                NO MORE THAN 72 CHARACTERS.
C        NMESSG- THE ACTUAL NUMBER OF CHARACTERS IN MESSG.
C        NERR  - THE ERROR NUMBER ASSOCIATED WITH THIS MESSAGE.
C                NERR MUST NOT BE ZERO.
C        LEVEL - ERROR CATEGORY.
C                =2 MEANS THIS IS AN UNCONDITIONALLY FATAL ERROR.
C                =1 MEANS THIS IS A RECOVERABLE ERROR.  (I.E., IT IS
C                   NON-FATAL IF XSETF HAS BEEN APPROPRIATELY CALLED.)
C                =0 MEANS THIS IS A WARNING MESSAGE ONLY.
C                =-1 MEANS THIS IS A WARNING MESSAGE WHICH IS TO BE
C                   PRINTED AT MOST ONCE, REGARDLESS OF HOW MANY
C                   TIMES THIS CALL IS EXECUTED.
C
C     EXAMPLES
C        CALL XERROR(23HSMOOTH -- NUM WAS ZERO.,23,1,2)
C        CALL XERROR(43HINTEG  -- LESS THAN FULL ACCURACY ACHIEVED.,
C                    43,2,1)
C        CALL XERROR(65HROOTER -- ACTUAL ZERO OF F FOUND BEFORE INTERVAL
C    1 FULLY COLLAPSED.,65,3,0)
C        CALL XERROR(39HEXP    -- UNDERFLOWS BEING SET TO ZERO.,39,1,-1)
C
C     WRITTEN BY RON JONES, WITH SLATEC COMMON MATH LIBRARY SUBCOMMITTEE
C     LATEST REVISION ---  7 FEB 1979
C
C***REFERENCES
C   JONES R.E., *SLATEC COMMON MATHEMATICAL LIBRARY ERROR HANDLING
C    PACKAGE*, SAND78-1189, SANDIA LABORATORIES, 1978.
C***ROUTINES CALLED  XERRWV
C***END PROLOGUE  XERROR
C
      implicit double precision (a-h,o-z)
      DIMENSION MESSG(NMESSG)
C***FIRST EXECUTABLE STATEMENT  XERROR
      CALL XERRWV(MESSG,NMESSG,NERR,LEVEL,0,0,0,0,0.,0.)
      RETURN
      END
c
c
c
*DECK CDIV
      SUBROUTINE CDIV(AR,AI,BR,BI,CR,CI)
C***BEGIN PROLOGUE  CDIV
C***REFER TO  EISDOC
C     COMPLEX DIVISION, (CR,CI) = (AR,AI)/(BR,BI)
C
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  CDIV
      implicit double precision (a-h,o-z)
      REAL*8 AR,AI,BR,BI,CR,CI
C
      REAL*8 S,ARS,AIS,BRS,BIS
C***FIRST EXECUTABLE STATEMENT  CDIV
      S = ABS(BR) + ABS(BI)
      ARS = AR/S
      AIS = AI/S
      BRS = BR/S
      BIS = BI/S
      S = BRS**2 + BIS**2
      CR = (ARS*BRS + AIS*BIS)/S
      CI = (AIS*BRS - ARS*BIS)/S
      RETURN
      END
c
c
c
*DECK XERRWV
      SUBROUTINE XERRWV(MESSG,NMESSG,NERR,LEVEL,NI,I1,I2,NR,R1,R2)
C***BEGIN PROLOGUE  XERRWV
C***REVISION DATE  811015   (YYMMDD)
C***CATEGORY NO.  Q
C***KEYWORDS  XERROR,ERROR,PRINT WITH VALUES
C***DATE WRITTEN  MARCH 19, 1980
C***AUTHOR  JONES R.E. (SLA)
C***PURPOSE
C PROCESSES ERROR MESSAGE ALLOWING 2 INTEGER AND TWO REAL VALUES
C  TO BE INCLUDED IN THE MESSAGE.
C***DESCRIPTION
C     ABSTRACT
C        XERRWV PROCESSES A DIAGNOSTIC MESSAGE, IN A MANNER
C        DETERMINED BY THE VALUE OF LEVEL AND THE CURRENT VALUE
C        OF THE LIBRARY ERROR CONTROL FLAG, KONTRL.
C        (SEE SUBROUTINE XSETF FOR DETAILS.)
C        IN ADDITION, UP TO TWO INTEGER VALUES AND TWO REAL
C        VALUES MAY BE PRINTED ALONG WITH THE MESSAGE.
C
C     DESCRIPTION OF PARAMETERS
C      --INPUT--
C        MESSG - THE HOLLERITH MESSAGE TO BE PROCESSED.
C        NMESSG- THE ACTUAL NUMBER OF CHARACTERS IN MESSG.
C        NERR  - THE ERROR NUMBER ASSOCIATED WITH THIS MESSAGE.
C                NERR MUST NOT BE ZERO.
C        LEVEL - ERROR CATEGORY.
C                =2 MEANS THIS IS AN UNCONDITIONALLY FATAL ERROR.
C                =1 MEANS THIS IS A RECOVERABLE ERROR.  (I.E., IT IS
C                   NON-FATAL IF XSETF HAS BEEN APPROPRIATELY CALLED.)
C                =0 MEANS THIS IS A WARNING MESSAGE ONLY.
C                =-1 MEANS THIS IS A WARNING MESSAGE WHICH IS TO BE
C                   PRINTED AT MOST ONCE, REGARDLESS OF HOW MANY
C                   TIMES THIS CALL IS EXECUTED.
C        NI    - NUMBER OF INTEGER VALUES TO BE PRINTED. (O TO 2)
C        I1    - FIRST INTEGER VALUE.
C        I2    - SECOND INTEGER VALUE.
C        NR    - NUMBER OF REAL VALUES TO BE PRINTED. (0 TO 2)
C        R1    - FIRST REAL VALUE.
C        R2    - SECOND REAL VALUE.
C
C     EXAMPLES
C        CALL XERROR(29HSMOOTH -- NUM (=I1) WAS ZERO.,29,1,2,
C    1   1,NUM,0,0,0.,0.)
C        CALL XERRWV(54HQUADXY -- REQUESTED ERROR (R1) LESS THAN MINIMUM
C    1 (R2).,54,77,1,0,0,0,2,ERRREQ,ERRMIN)
C
C     WRITTEN BY RON JONES, WITH SLATEC COMMON MATH LIBRARY SUBCOMMITTEE
C     LATEST REVISION ---  19 MAR 1980
C
C***REFERENCES
C   JONES R.E., *SLATEC COMMON MATHEMATICAL LIBRARY ERROR HANDLING
C    PACKAGE*, SAND78-1189, SANDIA LABORATORIES, 1978.
C***ROUTINES CALLED  FDUMP,I1MACH,J4SAVE,XERABT,XERCTL,XERPRT
C                      J4SAVE
C***END PROLOGUE  XERRWV
C
      implicit double precision (a-h,o-z)
      DIMENSION MESSG(NMESSG),LUN(5)
C***FIRST EXECUTABLE STATEMENT  XERRWV
      LKNTRL = J4SAVE(2,0,.FALSE.)
      MAXMES = J4SAVE(4,0,.FALSE.)
C     CHECK FOR VALID INPUT
      IF ((NMESSG.GT.0).AND.(NERR.NE.0).AND.
     1    (LEVEL.GE.(-1)).AND.(LEVEL.LE.2)) GO TO 10
         IF (LKNTRL.GT.0) CALL XERPRT(17HFATAL ERROR IN...,17)
         CALL XERPRT(23HXERROR -- INVALID INPUT,23)
         IF (LKNTRL.GT.0) CALL FDUMP
         IF (LKNTRL.GT.0) CALL XERPRT(29HJOB ABORT DUE TO FATAL ERROR.,
     1   29)
         IF (LKNTRL.GT.0) CALL XERSAV(1H ,0,0,0,KDUMMY)
         CALL XERABT(23HXERROR -- INVALID INPUT,23)
         RETURN
   10 CONTINUE
C     RECORD MESSAGE
      JUNK = J4SAVE(1,NERR,.TRUE.)
      CALL XERSAV(MESSG,NMESSG,NERR,LEVEL,KOUNT)
C     LET USER OVERRIDE
      LFIRST = MESSG(1)
      LMESSG = NMESSG
      LERR = NERR
      LLEVEL = LEVEL
      CALL XERCTL(LFIRST,LMESSG,LERR,LLEVEL,LKNTRL)
C     RESET TO ORIGINAL VALUES
      LMESSG = NMESSG
      LERR = NERR
      LLEVEL = LEVEL
      LKNTRL = MAX0(-2,MIN0(2,LKNTRL))
      MKNTRL = IABS(LKNTRL)
C     DECIDE WHETHER TO PRINT MESSAGE
      IF ((LLEVEL.LT.2).AND.(LKNTRL.EQ.0)) GO TO 100
      IF (((LLEVEL.EQ.(-1)).AND.(KOUNT.GT.MIN0(1,MAXMES)))
     1.OR.((LLEVEL.EQ.0)   .AND.(KOUNT.GT.MAXMES))
     2.OR.((LLEVEL.EQ.1)   .AND.(KOUNT.GT.MAXMES).AND.(MKNTRL.EQ.1))
     3.OR.((LLEVEL.EQ.2)   .AND.(KOUNT.GT.MAX0(1,MAXMES)))) GO TO 100
         IF (LKNTRL.LE.0) GO TO 20
            CALL XERPRT(1H ,1)
C           INTRODUCTION
            IF (LLEVEL.EQ.(-1)) CALL XERPRT
     1(57HWARNING MESSAGE...THIS MESSAGE WILL ONLY BE PRINTED ONCE.,57)
            IF (LLEVEL.EQ.0) CALL XERPRT(13HWARNING IN...,13)
            IF (LLEVEL.EQ.1) CALL XERPRT
     1      (23HRECOVERABLE ERROR IN...,23)
            IF (LLEVEL.EQ.2) CALL XERPRT(17HFATAL ERROR IN...,17)
   20    CONTINUE
C        MESSAGE
         CALL XERPRT(MESSG,LMESSG)
         CALL XGETUA(LUN,NUNIT)
         DO 50 KUNIT=1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
            IF (NI.GE.1) WRITE (IUNIT,22) I1
            IF (NI.GE.2) WRITE (IUNIT,23) I2
            IF (NR.GE.1) WRITE (IUNIT,24) R1
            IF (NR.GE.2) WRITE (IUNIT,25) R2
   22       FORMAT (11X,21HIN ABOVE MESSAGE, I1=,I10)
   23       FORMAT (11X,21HIN ABOVE MESSAGE, I2=,I10)
   24       FORMAT (11X,21HIN ABOVE MESSAGE, R1=,E20.10)
   25       FORMAT (11X,21HIN ABOVE MESSAGE, R2=,E20.10)
            IF (LKNTRL.LE.0) GO TO 40
C              ERROR NUMBER
               WRITE (IUNIT,30) LERR
   30          FORMAT (15H ERROR NUMBER =,I10)
   40       CONTINUE
   50    CONTINUE
C        TRACE-BACK
         CALL FDUMP
  100 CONTINUE
      IFATAL = 0
      IF ((LLEVEL.EQ.2).OR.((LLEVEL.EQ.1).AND.(MKNTRL.EQ.2)))
     1IFATAL = 1
C     QUIT HERE IF MESSAGE IS NOT FATAL
      IF (IFATAL.LE.0) RETURN
      IF ((LKNTRL.LE.0).OR.(KOUNT.GT.MAX0(1,MAXMES))) GO TO 120
C        PRINT REASON FOR ABORT
         IF (LLEVEL.EQ.1) CALL XERPRT
     1   (35HJOB ABORT DUE TO UNRECOVERED ERROR.,35)
         IF (LLEVEL.EQ.2) CALL XERPRT
     1   (29HJOB ABORT DUE TO FATAL ERROR.,29)
C        PRINT ERROR SUMMARY
         CALL XERSAV(1H ,-1,0,0,KDUMMY)
  120 CONTINUE
C     ABORT
      IF ((LLEVEL.EQ.2).AND.(KOUNT.GT.MAX0(1,MAXMES))) LMESSG = 0
      CALL XERABT(MESSG,LMESSG)
      RETURN
      END
c
c
c
*DECK FDUMP
      SUBROUTINE FDUMP
C***BEGIN PROLOGUE  FDUMP
C***REVISION DATE  811015   (YYMMDD)
C***CATEGORY NO.  Q2,N2
C***KEYWORDS  SYMBOLIC DUMP,DUMP
C***DATE WRITTEN  8/79
C***AUTHOR  JONES R.E. (SLA)
C***PURPOSE
C  SYMBOLIC DUMP (SHOULD BE LOCALLY WRITTEN)
C***DESCRIPTION
C        ***NOTE*** MACHINE DEPENDENT ROUTINE
C        FDUMP IS INTENDED TO BE REPLACED BY A LOCALLY WRITTEN
C        VERSION WHICH PRODUCES A SYMBOLIC DUMP.  FAILING THIS,
C        IT SHOULD BE REPLACED BY A VERSION WHICH PRINTS THE
C        SUBPROGRAM NESTING LIST.  NOTE THAT THIS DUMP MUST BE
C        PRINTED ON EACH OF UP TO FIVE FILES, AS INDICATED BY THE
C        XGETUA ROUTINE.  SEE XSETUA AND XGETUA FOR DETAILS.
C
C     WRITTEN BY RON JONES, WITH SLATEC COMMON MATH LIBRARY SUBCOMMITTEE
C     LATEST REVISION ---  23 MAY 1979
C
C***REFERENCES
C   JONES R.E., *SLATEC COMMON MATHEMATICAL LIBRARY ERROR HANDLING
C    PACKAGE*, SAND78-1189, SANDIA LABORATORIES, 1978.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  FDUMP
C
C***FIRST EXECUTABLE STATEMENT  FDUMP
      RETURN
      END
c
c
c
*DECK I1MACH
      INTEGER FUNCTION I1MACH(I)
C***BEGIN PROLOGUE  I1MACH
C***REVISION DATE  811015   (YYMMDD)
C***CATEGORY NO.  Q
C***KEYWORDS  MACHINE CONSTANTS,INTEGER
C***DATE WRITTEN   1975
C***AUTHOR FOX P.A.,HALL A.D.,SCHRYER N.L. (BELL LABS)
C***PURPOSE
C RETURNS INTEGER MACHINE DEPENDENT CONSTANTS
C***DESCRIPTION
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C   THESE MACHINE CONSTANT ROUTINES MUST BE ACTIVATED FOR
C   A PARTICULAR ENVIRONMENT.
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     I1MACH CAN BE USED TO OBTAIN MACHINE-DEPENDENT PARAMETERS
C     FOR THE LOCAL MACHINE ENVIRONMENT.  IT IS A FUNCTION
C     SUBROUTINE WITH ONE (INPUT) ARGUMENT, AND CAN BE CALLED
C     AS FOLLOWS, FOR EXAMPLE
C
C          K = I1MACH(I)
C
C     WHERE I=1,...,16.  THE (OUTPUT) VALUE OF K ABOVE IS
C     DETERMINED BY THE (INPUT) VALUE OF I.  THE RESULTS FOR
C     VARIOUS VALUES OF I ARE DISCUSSED BELOW.
C
C  I/O UNIT NUMBERS.
C    I1MACH( 1) = THE STANDARD INPUT UNIT.
C    I1MACH( 2) = THE STANDARD OUTPUT UNIT.
C    I1MACH( 3) = THE STANDARD PUNCH UNIT.
C    I1MACH( 4) = THE STANDARD ERROR MESSAGE UNIT.
C
C  WORDS.
C    I1MACH( 5) = THE NUMBER OF BITS PER INTEGER STORAGE UNIT.
C    I1MACH( 6) = THE NUMBER OF CHARACTERS PER INTEGER STORAGE UNIT.
C
C  INTEGERS.
C    ASSUME INTEGERS ARE REPRESENTED IN THE S-DIGIT, BASE-A FORM
C
C               SIGN ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
C
C               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,S-1.
C    I1MACH( 7) = A, THE BASE.
C    I1MACH( 8) = S, THE NUMBER OF BASE-A DIGITS.
C    I1MACH( 9) = A**S - 1, THE LARGEST MAGNITUDE.
C
C  FLOATING-POINT NUMBERS.
C    ASSUME FLOATING-POINT NUMBERS ARE REPRESENTED IN THE T-DIGIT,
C    BASE-B FORM
C               SIGN (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C               WHERE 0 .LE. X(I) .LT. B FOR I=1,...,T,
C               0 .LT. X(1), AND EMIN .LE. E .LE. EMAX.
C    I1MACH(10) = B, THE BASE.
C
C  SINGLE-PRECISION
C    I1MACH(11) = T, THE NUMBER OF BASE-B DIGITS.
C    I1MACH(12) = EMIN, THE SMALLEST EXPONENT E.
C    I1MACH(13) = EMAX, THE LARGEST EXPONENT E.
C
C  DOUBLE-PRECISION
C    I1MACH(14) = T, THE NUMBER OF BASE-B DIGITS.
C    I1MACH(15) = EMIN, THE SMALLEST EXPONENT E.
C    I1MACH(16) = EMAX, THE LARGEST EXPONENT E.
C
C  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
C  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY
C  REMOVING THE C FROM COLUMN 1.  ALSO, THE VALUES OF
C  I1MACH(1) - I1MACH(4) SHOULD BE CHECKED FOR CONSISTENCY
C  WITH THE LOCAL OPERATING SYSTEM.
C
C***REFERENCES
C  FOX P.A., HALL A.D., SCHRYER N.L.,*FRAMEWORK FOR A PORTABLE LIBRARY*,
C  ACM TRANSACTION ON MATHEMATICAL SOFTWARE, VOL. 4, NO. 2,
C  JUNE 1978, PP. 177-188.
C***ROUTINES CALLED  XERROR
C***END PROLOGUE  I1MACH
C
      INTEGER IMACH(16),OUTPUT
C
C     EQUIVALENCE (IMACH(4),OUTPUT)
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
C
C     DATA IMACH( 1) /    7 /
C     DATA IMACH( 2) /    2 /
C     DATA IMACH( 3) /    2 /
C     DATA IMACH( 4) /    2 /
C     DATA IMACH( 5) /   36 /
C     DATA IMACH( 6) /    4 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   33 /
C     DATA IMACH( 9) / Z1FFFFFFFF /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   24 /
C     DATA IMACH(12) / -256 /
C     DATA IMACH(13) /  255 /
C     DATA IMACH(14) /   60 /
C     DATA IMACH(15) / -256 /
C     DATA IMACH(16) /  255 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.
C
C     DATA IMACH( 1) /   5 /
C     DATA IMACH( 2) /   6 /
C     DATA IMACH( 3) /   7 /
C     DATA IMACH( 4) /   6 /
C     DATA IMACH( 5) /  48 /
C     DATA IMACH( 6) /   6 /
C     DATA IMACH( 7) /   2 /
C     DATA IMACH( 8) /  39 /
C     DATA IMACH( 9) / O0007777777777777 /
C     DATA IMACH(10) /   8 /
C     DATA IMACH(11) /  13 /
C     DATA IMACH(12) / -50 /
C     DATA IMACH(13) /  76 /
C     DATA IMACH(14) /  26 /
C     DATA IMACH(15) / -50 /
C     DATA IMACH(16) /  76 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.
C
C     DATA IMACH( 1) /   5 /
C     DATA IMACH( 2) /   6 /
C     DATA IMACH( 3) /   7 /
C     DATA IMACH( 4) /   6 /
C     DATA IMACH( 5) /  48 /
C     DATA IMACH( 6) /   6 /
C     DATA IMACH( 7) /   2 /
C     DATA IMACH( 8) /  39 /
C     DATA IMACH( 9) / O0007777777777777 /
C     DATA IMACH(10) /   8 /
C     DATA IMACH(11) /  13 /
C     DATA IMACH(12) / -50 /
C     DATA IMACH(13) /  76 /
C     DATA IMACH(14) /  26 /
C     DATA IMACH(15) / -32754 /
C     DATA IMACH(16) /  32780 /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    7 /
C     DATA IMACH( 4) /6LOUTPUT/
C     DATA IMACH( 5) /   60 /
C     DATA IMACH( 6) /   10 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   48 /
C     DATA IMACH( 9) / 00007777777777777777B /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   48 /
C     DATA IMACH(12) / -974 /
C     DATA IMACH(13) / 1070 /
C     DATA IMACH(14) /   96 /
C     DATA IMACH(15) / -927 /
C     DATA IMACH(16) / 1070 /
C
C     MACHINE CONSTANTS FOR THE CRAY 1
C
C     DATA IMACH( 1) /   100 /
C     DATA IMACH( 2) /   101 /
C     DATA IMACH( 3) /   102 /
C     DATA IMACH( 4) /   101 /
C     DATA IMACH( 5) /    64 /
C     DATA IMACH( 6) /     8 /
C     DATA IMACH( 7) /     2 /
C     DATA IMACH( 8) /    63 /
C     DATA IMACH( 9) /  777777777777777777777B /
C     DATA IMACH(10) /     2 /
C     DATA IMACH(11) /    48 /
C     DATA IMACH(12) / -8192 /
C     DATA IMACH(13) /  8191 /
C     DATA IMACH(14) /    96 /
C     DATA IMACH(15) / -8192 /
C     DATA IMACH(16) /  8191 /
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
C
C     DATA IMACH( 1) /   11 /
C     DATA IMACH( 2) /   12 /
C     DATA IMACH( 3) /    8 /
C     DATA IMACH( 4) /   10 /
C     DATA IMACH( 5) /   16 /
C     DATA IMACH( 6) /    2 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   15 /
C     DATA IMACH( 9) /32767 /
C     DATA IMACH(10) /   16 /
C     DATA IMACH(11) /    6 /
C     DATA IMACH(12) /  -64 /
C     DATA IMACH(13) /   63 /
C     DATA IMACH(14) /   14 /
C     DATA IMACH(15) /  -64 /
C     DATA IMACH(16) /   63 /
C
C     MACHINE CONSTANTS FOR THE HARRIS 220
C
C     DATA IMACH( 1) /       5 /
C     DATA IMACH( 2) /       6 /
C     DATA IMACH( 3) /       0 /
C     DATA IMACH( 4) /       6 /
C     DATA IMACH( 5) /      24 /
C     DATA IMACH( 6) /       3 /
C     DATA IMACH( 7) /       2 /
C     DATA IMACH( 8) /      23 /
C     DATA IMACH( 9) / 8388607 /
C     DATA IMACH(10) /       2 /
C     DATA IMACH(11) /      23 /
C     DATA IMACH(12) /    -127 /
C     DATA IMACH(13) /     127 /
C     DATA IMACH(14) /      38 /
C     DATA IMACH(15) /    -127 /
C     DATA IMACH(16) /     127 /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES.
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /   43 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   36 /
C     DATA IMACH( 6) /    6 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   35 /
C     DATA IMACH( 9) / O377777777777 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   27 /
C     DATA IMACH(12) / -127 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   63 /
C     DATA IMACH(15) / -127 /
C     DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     3 WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA IMACH(1) /      5/
C     DATA IMACH(2) /      6 /
C     DATA IMACH(3) /      4 /
C     DATA IMACH(4) /      1 /
C     DATA IMACH(5) /     16 /
C     DATA IMACH(6) /      2 /
C     DATA IMACH(7) /      2 /
C     DATA IMACH(8) /     15 /
C     DATA IMACH(9) /  32767 /
C     DATA IMACH(10)/      2 /
C     DATA IMACH(11)/     23 /
C     DATA IMACH(12)/   -128 /
C     DATA IMACH(13)/    127 /
C     DATA IMACH(14)/     39 /
C     DATA IMACH(15)/   -128 /
C     DATA IMACH(16)/    127 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     4 WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA IMACH(1) /      5 /
C     DATA IMACH(2) /      6 /
C     DATA IMACH(3) /      4 /
C     DATA IMACH(4) /      1 /
C     DATA IMACH(5) /     16 /
C     DATA IMACH(6) /      2 /
C     DATA IMACH(7) /      2 /
C     DATA IMACH(8) /     15 /
C     DATA IMACH(9) /  32767 /
C     DATA IMACH(10)/      2 /
C     DATA IMACH(11)/     23 /
C     DATA IMACH(12)/   -128 /
C     DATA IMACH(13)/    127 /
C     DATA IMACH(14)/     55 /
C     DATA IMACH(15)/   -128 /
C     DATA IMACH(16)/    127 /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
C     THE PERKIN ELMER (INTERDATA) 7/32.
C
C     DATA IMACH( 1) /   5 /
C     DATA IMACH( 2) /   6 /
C     DATA IMACH( 3) /   7 /
C     DATA IMACH( 4) /   6 /
C     DATA IMACH( 5) /  32 /
C     DATA IMACH( 6) /   4 /
C     DATA IMACH( 7) /   2 /
C     DATA IMACH( 8) /  31 /
C     DATA IMACH( 9) / Z7FFFFFFF /
C     DATA IMACH(10) /  16 /
C     DATA IMACH(11) /   6 /
C     DATA IMACH(12) / -64 /
C     DATA IMACH(13) /  63 /
C     DATA IMACH(14) /  14 /
C     DATA IMACH(15) / -64 /
C     DATA IMACH(16) /  63 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    5 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   36 /
C     DATA IMACH( 6) /    5 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   35 /
C     DATA IMACH( 9) / "377777777777 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   27 /
C     DATA IMACH(12) / -128 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   54 /
C     DATA IMACH(15) / -101 /
C     DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    5 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   36 /
C     DATA IMACH( 6) /    5 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   35 /
C     DATA IMACH( 9) / "377777777777 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   27 /
C     DATA IMACH(12) / -128 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   62 /
C     DATA IMACH(15) / -128 /
C     DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN S SUPPORTING
C     32-BIT INTEGER ARITHMETIC.
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    5 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   32 /
C     DATA IMACH( 6) /    4 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   31 /
C     DATA IMACH( 9) / 2147483647 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   24 /
C     DATA IMACH(12) / -127 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   56 /
C     DATA IMACH(15) / -127 /
C     DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN S SUPPORTING
C     16-BIT INTEGER ARITHMETIC.
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    5 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   16 /
C     DATA IMACH( 6) /    2 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   15 /
C     DATA IMACH( 9) / 32767 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   24 /
C     DATA IMACH(12) / -127 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   56 /
C     DATA IMACH(15) / -127 /
C     DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES. FTN COMPILER
C
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    1 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   36 /
C     DATA IMACH( 6) /    4 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   35 /
C     DATA IMACH( 9) / O377777777777 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   27 /
C     DATA IMACH(12) / -128 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   60 /
C     DATA IMACH(15) /-1024 /
C     DATA IMACH(16) / 1023 /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES. FOR COMPILER
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    7 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   36 /
C     DATA IMACH( 6) /    6 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   35 /
C     DATA IMACH( 9) / O377777777777 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   27 /
C     DATA IMACH(12) / -128 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   60 /
C     DATA IMACH(15) /-1024/
C     DATA IMACH(16) / 1023 /
C
C
C     MACHINE CONSTANTS FOR THE VAX 11/780
C
C     DATA IMACH(1) /    5 /
C     DATA IMACH(2) /    6 /
C     DATA IMACH(3) /    5 /
C     DATA IMACH(4) /    6 /
C     DATA IMACH(5) /   32 /
C     DATA IMACH(6) /    4 /
C     DATA IMACH(7) /    2 /
C     DATA IMACH(8) /   31 /
C     DATA IMACH(9) /2147483647 /
C     DATA IMACH(10)/    2 /
C     DATA IMACH(11)/   24 /
C     DATA IMACH(12)/ -127 /
C     DATA IMACH(13)/  127 /
C     DATA IMACH(14)/   56 /
C     DATA IMACH(15)/ -127 /
C     DATA IMACH(16)/  127 /
C
C
C     MACHINE CONSTANTS FOR THE APOLLO DOMAIN SYSTEM
C
      DATA IMACH(1) /    5 /
      DATA IMACH(2) /    6 /
      DATA IMACH(3) /    5 /
      DATA IMACH(4) /    6 /
      DATA IMACH(5) /   32 /
      DATA IMACH(6) /    4 /
      DATA IMACH(7) /    2 /
      DATA IMACH(8) /   31 /
      DATA IMACH(9) /2147483647 /
      DATA IMACH(10)/    2 /
      DATA IMACH(11)/   23 /
      DATA IMACH(12)/ -127 /
      DATA IMACH(13)/  127 /
      DATA IMACH(14)/   52 /
      DATA IMACH(15)/-1023 /
      DATA IMACH(16)/ 1023 /
C
C***FIRST EXECUTABLE STATEMENT  I1MACH
C
      IF (I .LT. 1  .OR.  I .GT. 16) GO TO 10
C
      I1MACH=IMACH(I)
      RETURN
C
   10 CONTINUE
      WRITE(OUTPUT,9000)
9000  FORMAT(39H1ERROR    1 IN I1MACH - I OUT OF BOUNDS )
C
C     CALL FDUMP
C
C
      STOP
      END
c
c
c
*DECK J4SAVE
      FUNCTION J4SAVE(IWHICH,IVALUE,ISET)
C***BEGIN PROLOGUE  J4SAVE
C***REFER TO  XERROR
C     ABSTRACT
C        J4SAVE SAVES AND RECALLS SEVERAL GLOBAL VARIABLES NEEDED
C        BY THE LIBRARY ERROR HANDLING ROUTINES.
C
C     DESCRIPTION OF PARAMETERS
C      --INPUT--
C        IWHICH - INDEX OF ITEM DESIRED.
C                 = 1 REFERS TO CURRENT ERROR NUMBER.
C                 = 2 REFERS TO CURRENT ERROR CONTROL FLAG.
C                 = 3 REFERS TO CURRENT UNIT NUMBER TO WHICH ERROR
C                     MESSAGES ARE TO BE SENT.  (0 MEANS USE STANDARD.)
C                 = 4 REFERS TO THE MAXIMUM NUMBER OF TIMES ANY
C                     MESSAGE IS TO BE PRINTED (AS SET BY XERMAX).
C                 = 5 REFERS TO THE TOTAL NUMBER OF UNITS TO WHICH
C                     EACH ERROR MESSAGE IS TO BE WRITTEN.
C                 = 6 REFERS TO THE 2ND UNIT FOR ERROR MESSAGES
C                 = 7 REFERS TO THE 3RD UNIT FOR ERROR MESSAGES
C                 = 8 REFERS TO THE 4TH UNIT FOR ERROR MESSAGES
C                 = 9 REFERS TO THE 5TH UNIT FOR ERROR MESSAGES
C        IVALUE - THE VALUE TO BE SET FOR THE IWHICH-TH PARAMETER,
C                 IF ISET IS .TRUE. .
C        ISET   - IF ISET=.TRUE., THE IWHICH-TH PARAMETER WILL BE
C                 GIVEN THE VALUE, IVALUE.  IF ISET=.FALSE., THE
C                 IWHICH-TH PARAMETER WILL BE UNCHANGED, AND IVALUE
C                 IS A DUMMY PARAMETER.
C      --OUTPUT--
C        THE (OLD) VALUE OF THE IWHICH-TH PARAMETER WILL BE RETURNED
C        IN THE FUNCTION VALUE, J4SAVE.
C
C     WRITTEN BY RON JONES, WITH SLATEC COMMON MATH LIBRARY SUBCOMMITTEE
C     ADAPTED FROM BELL LABORATORIES PORT LIBRARY ERROR HANDLER
C     LATEST REVISION ---  23 MAY 1979
C
C***REFERENCES
C   JONES R.E., *SLATEC COMMON MATHEMATICAL LIBRARY ERROR HANDLING
C    PACKAGE*, SAND78-1189, SANDIA LABORATORIES, 1978.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  J4SAVE
C
      LOGICAL ISET
      INTEGER IPARAM(9)
      DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,1,0,10/
      DATA IPARAM(5)/1/
      DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/
C***FIRST EXECUTABLE STATEMENT  J4SAVE
      J4SAVE = IPARAM(IWHICH)
      IF (ISET) IPARAM(IWHICH) = IVALUE
      RETURN
      END
c
c
c
*DECK XERABT
      SUBROUTINE XERABT(MESSG,NMESSG)
C***BEGIN PROLOGUE  XERABT
C***REVISION DATE  811015   (YYMMDD)
C***CATEGORY NO.  Q
C***KEYWORDS  XERROR,ERROR,ABORT
C***DATE WRITTEN  AUGUST 1979
C***AUTHOR  JONES R.E. (SLA)
C***PURPOSE
C ABORTS PROGRAM EXECUTION AND PRINTS ERROR MESSAGE
C***DESCRIPTION
C     ABSTRACT
C        ***NOTE*** MACHINE DEPENDENT ROUTINE
C        XERABT ABORTS THE EXECUTION OF THE PROGRAM.
C        THE ERROR MESSAGE CAUSING THE ABORT IS GIVEN IN THE CALLING
C        SEQUENCE IN CASE ONE NEEDS IT FOR PRINTING ON A DAYFILE,
C        FOR EXAMPLE.
C
C     DESCRIPTION OF PARAMETERS
C        MESSG AND NMESSG ARE AS IN XERROR, EXCEPT THAT NMESSG MAY
C        BE ZERO, IN WHICH CASE NO MESSAGE IS BEING SUPPLIED.
C
C     WRITTEN BY RON JONES, WITH SLATEC COMMON MATH LIBRARY SUBCOMMITTEE
C     LATEST REVISION ---  7 JUNE 1978
C
C***REFERENCES
C   JONES R.E., *SLATEC COMMON MATHEMATICAL LIBRARY ERROR HANDLING
C    PACKAGE*, SAND78-1189, SANDIA LABORATORIES, 1978.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  XERABT
C
      implicit double precision (a-h,o-z)
      DIMENSION MESSG(NMESSG)
C***FIRST EXECUTABLE STATEMENT  XERABT
      STOP
C     RETURN
      END
c
c
c
*DECK XERCTL
      SUBROUTINE XERCTL(MESSG1,NMESSG,NERR,LEVEL,KONTRL)
C***BEGIN PROLOGUE  XERCTL
C***REVISION DATE  811015   (YYMMDD)
C***CATEGORY NO.  Q
C***KEYWORDS  XERROR,ERROR,CONTROL
C***DATE WRITTEN  AUGUST 1979
C***AUTHOR  JONES R.E. (SLA)
C***PURPOSE
C ALLOWS USER CONTROL OVER HANDLING OF INDIVIDUAL ERRORS
C***DESCRIPTION
C     ABSTRACT
C        ALLOWS USER CONTROL OVER HANDLING OF INDIVIDUAL ERRORS.
C        JUST AFTER EACH MESSAGE IS RECORDED, BUT BEFORE IT IS
C        PROCESSED ANY FURTHER (I.E., BEFORE IT IS PRINTED OR
C        A DECISION TO ABORT IS MADE) A CALL IS MADE TO XERCTL.
C        IF THE USER HAS PROVIDED HIS OWN VERSION OF XERCTL, HE
C        CAN THEN OVERRIDE THE VALUE OF KONTROL USED IN PROCESSING
C        THIS MESSAGE BY REDEFINING ITS VALUE.
C        KONTRL MAY BE SET TO ANY VALUE FROM -2 TO 2.
C        THE MEANINGS FOR KONTRL ARE THE SAME AS IN XSETF, EXCEPT
C        THAT THE VALUE OF KONTRL CHANGES ONLY FOR THIS MESSAGE.
C        IF KONTRL IS SET TO A VALUE OUTSIDE THE RANGE FROM -2 TO 2,
C        IT WILL BE MOVED BACK INTO THAT RANGE.
C
C     DESCRIPTION OF PARAMETERS
C
C      --INPUT--
C        MESSG1 - THE FIRST WORD (ONLY) OF THE ERROR MESSAGE.
C        NMESSG - SAME AS IN THE CALL TO XERROR OR XERRWV.
C        NERR   - SAME AS IN THE CALL TO XERROR OR XERRWV.
C        LEVEL  - SAME AS IN THE CALL TO XERROR OR XERRWV.
C        KONTRL - THE CURRENT VALUE OF THE CONTROL FLAG AS SET
C                 BY A CALL TO XSETF.
C
C      --OUTPUT--
C        KONTRL - THE NEW VALUE OF KONTRL.  IF KONTRL IS NOT
C                 DEFINED, IT WILL REMAIN AT ITS ORIGINAL VALUE.
C                 THIS CHANGED VALUE OF CONTROL AFFECTS ONLY
C                 THE CURRENT OCCURRENCE OF THE CURRENT MESSAGE.
C
C***REFERENCES
C   JONES R.E., *SLATEC COMMON MATHEMATICAL LIBRARY ERROR HANDLING
C    PACKAGE*, SAND78-1189, SANDIA LABORATORIES, 1978.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  XERCTL
C
C***FIRST EXECUTABLE STATEMENT  XERCTL
      RETURN
      END
c
c
c
*DECK XERPRT
      SUBROUTINE XERPRT(MESSG,NMESSG)
C***BEGIN PROLOGUE  XERPRT
C***REVISION DATE  811015   (YYMMDD)
C***CATEGORY NO.  Q
C***KEYWORDS  XERROR,ERROR,PRINT
C***DATE WRITTEN  AUGUST 1979
C***AUTHOR  JONES R.E. (SLA)
C***PURPOSE
C PRINTS ERROR MESSAGES
C***DESCRIPTION
C     ABSTRACT
C        PRINT THE HOLLERITH MESSAGE IN MESSG, OF LENGTH NMESSG,
C        ON EACH FILE INDICATED BY XGETUA.
C
C     PREPARE FORMAT
C***REFERENCES
C   JONES R.E., *SLATEC COMMON MATHEMATICAL LIBRARY ERROR HANDLING
C    PACKAGE*, SAND78-1189, SANDIA LABORATORIES, 1978.
C***ROUTINES CALLED  XGETUA,S88FMT,I1MACH
C***END PROLOGUE  XERPRT
C
      implicit double precision (a-h,o-z)
      INTEGER F(10),G(14),LUN(5)
      DIMENSION MESSG(NMESSG)
      DATA F(1),F(2),F(3),F(4),F(5),F(6),F(7),F(8),F(9),F(10)
     1   / 1H( ,1H1 ,1HX ,1H, ,1H  ,1H  ,1HA ,1H  ,1H  ,1H) /
      DATA G(1),G(2),G(3),G(4),G(5),G(6),G(7),G(8),G(9),G(10)
     1   / 1H( ,1H1 ,1HX ,1H  ,1H  ,1H  ,1H  ,1H  ,1H  ,1H  /
      DATA G(11),G(12),G(13),G(14)
     1   / 1H   ,1H   ,1H   ,1H)  /
      DATA LA/1HA/,LCOM/1H,/,LBLANK/1H /
C    PREPARE FORMAT FOR WHOLE LINES
C***FIRST EXECUTABLE STATEMENT  XERPRT
      NCHAR = I1MACH(6)
      NFIELD = 72/NCHAR
      CALL S88FMT(2,NFIELD,F(5))
      CALL S88FMT(2,NCHAR,F(8))
C     PREPARE FORMAT FOR LAST, PARTIAL LINE, IF NEEDED
      NCHARL = NFIELD*NCHAR
      NLINES = NMESSG/NCHARL
      NWORD  = NLINES*NFIELD
      NCHREM = NMESSG - NLINES*NCHARL
      IF (NCHREM.LE.0) GO TO 40
         DO 10 I=4,13
10          G(I) = LBLANK
         NFIELD = NCHREM/NCHAR
         IF (NFIELD.LE.0) GO TO 20
C        PREPARE WHOLE WORD FIELDS
            G(4) = LCOM
            CALL S88FMT(2,NFIELD,G(5))
            G(7) = LA
            CALL S88FMT(2,NCHAR,G(8))
20       CONTINUE
         NCHLST = MOD(NCHREM,NCHAR)
         IF (NCHLST.LE.0) GO TO 30
C        PREPARE PARTIAL WORD FIELD
            G(10) = LCOM
            G(11) = LA
            CALL S88FMT(2,NCHLST,G(12))
30       CONTINUE
40    CONTINUE
C     PRINT THE MESSAGE
      NWORD1 = NWORD+1
      NWORD2 = (NMESSG+NCHAR-1)/NCHAR
      CALL XGETUA(LUN,NUNIT)
      DO 50 KUNIT = 1,NUNIT
         IUNIT = LUN(KUNIT)
         IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
         IF (NWORD.GT.0) WRITE (IUNIT,F) (MESSG(I),I=1,NWORD)
         IF (NCHREM.GT.0) WRITE (IUNIT,G) (MESSG(I),I=NWORD1,NWORD2)
50    CONTINUE
      RETURN
      END
c
c
c
*DECK XGETUA
      SUBROUTINE XGETUA(IUNIT,N)
C***BEGIN PROLOGUE  XGETUA
C***REVISION DATE  811015   (YYMMDD)
C***CATEGORY NO.  Q
C***KEYWORDS  XERROR,ERROR,UNIT NUMBER
C***DATE WRITTEN  AUGUST 1979
C***AUTHOR  JONES R.E. (SLA)
C***PURPOSE
C RETURNS UNIT NUMBER(S) TO WHICH ERROR MESSAGES ARE BEING SENT
C***DESCRIPTION
C     ABSTRACT
C        XGETUA MAY BE CALLED TO DETERMINE THE UNIT NUMBER OR NUMBERS
C        TO WHICH ERROR MESSAGES ARE BEING SENT.
C        THESE UNIT NUMBERS MAY HAVE BEEN SET BY A CALL TO XSETUN,
C        OR A CALL TO XSETUA, OR MAY BE A DEFAULT VALUE.
C
C     DESCRIPTION OF PARAMETERS
C      --OUTPUT--
C        IUNIT - AN ARRAY OF ONE TO FIVE UNIT NUMBERS, DEPENDING
C                ON THE VALUE OF N.  A VALUE OF ZERO REFERS TO THE
C                DEFAULT UNIT, AS DEFINED BY THE I1MACH MACHINE
C                CONSTANT ROUTINE.  ONLY IUNIT(1),...,IUNIT(N) ARE
C                DEFINED BY XGETUA.  THE VALUES OF IUNIT(N+1),...,
C                IUNIT(5) ARE NOT DEFINED (FOR N.LT.5) OR ALTERED
C                IN ANY WAY BY XGETUA.
C        N     - THE NUMBER OF UNITS TO WHICH COPIES OF THE
C                ERROR MESSAGES ARE BEING SENT.  N WILL BE IN THE
C                RANGE FROM 1 TO 5.
C
C     WRITTEN BY RON JONES, WITH SLATEC COMMON MATH LIBRARY SUBCOMMITTEE
C
C***REFERENCES
C   JONES R.E., *SLATEC COMMON MATHEMATICAL LIBRARY ERROR HANDLING
C    PACKAGE*, SAND78-1189, SANDIA LABORATORIES, 1978.
C***ROUTINES CALLED   J4SAVE
C***END PROLOGUE  XGETUA
C
      implicit double precision (a-h,o-z)
      DIMENSION IUNIT(5)
C***FIRST EXECUTABLE STATEMENT  XGETUA
      N = J4SAVE(5,0,.FALSE.)
      DO 30 I=1,N
         INDEX = I+4
         IF (I.EQ.1) INDEX = 3
         IUNIT(I) = J4SAVE(INDEX,0,.FALSE.)
   30 CONTINUE
      RETURN
      END
c
c
c
*DECK S88FMT
      SUBROUTINE S88FMT(N,IVALUE,IFMT)
C***BEGIN PROLOGUE  S88FMT
C***REFER TO  XERROR
C     ABSTRACT
C        S88FMT REPLACES IFMT(1), ... ,IFMT(N) WITH THE
C        CHARACTERS CORRESPONDING TO THE N LEAST SIGNIFICANT
C        DIGITS OF IVALUE.
C
C     TAKEN FROM THE BELL LABORATORIES PORT LIBRARY ERROR HANDLER
C     LATEST REVISION ---  7 JUNE 1978
C
C***REFERENCES
C   JONES R.E., *SLATEC COMMON MATHEMATICAL LIBRARY ERROR HANDLING
C    PACKAGE*, SAND78-1189, SANDIA LABORATORIES, 1978.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  S88FMT
C
      implicit double precision (a-h,o-z)
      DIMENSION IFMT(N),IDIGIT(10)
      DATA IDIGIT(1),IDIGIT(2),IDIGIT(3),IDIGIT(4),IDIGIT(5),
     1     IDIGIT(6),IDIGIT(7),IDIGIT(8),IDIGIT(9),IDIGIT(10)
     2     /1H0,1H1,1H2,1H3,1H4,1H5,1H6,1H7,1H8,1H9/
C***FIRST EXECUTABLE STATEMENT  S88FMT
      NT = N
      IT = IVALUE
   10    IF (NT .EQ. 0) RETURN
         INDEX = MOD(IT,10)
         IFMT(NT) = IDIGIT(INDEX+1)
         IT = IT/10
         NT = NT - 1
         GO TO 10
      END
c
c
c
*DECK XERSAV
      SUBROUTINE XERSAV(MESSG,NMESSG,NERR,LEVEL,ICOUNT)
C***BEGIN PROLOGUE  XERSAV
C***REVISION DATE  811015   (YYMMDD)
C***CATEGORY NO.  Q
C***KEYWORDS  XERROR,ERROR,RECORD
C***DATE WRITTEN  MARCH 19, 1980
C***AUTHOR  JONES R.E. (SLA)
C***PURPOSE
C RECORDS THAT AN ERROR OCCURRED.
C***DESCRIPTION
C     ABSTRACT
C        RECORD THAT THIS ERROR OCCURRED.
C
C     DESCRIPTION OF PARAMETERS
C     --INPUT--
C       MESSG, NMESSG, NERR, LEVEL ARE AS IN XERROR,
C       EXCEPT THAT WHEN NMESSG=0 THE TABLES WILL BE
C       DUMPED AND CLEARED, AND WHEN NMESSG IS LESS THAN ZERO THE
C       TABLES WILL BE DUMPED AND NOT CLEARED.
C     --OUTPUT--
C       ICOUNT WILL BE THE NUMBER OF TIMES THIS MESSAGE HAS
C       BEEN SEEN, OR ZERO IF THE TABLE HAS OVERFLOWED AND
C       DOES NOT CONTAIN THIS MESSAGE SPECIFICALLY.
C       WHEN NMESSG=0, ICOUNT WILL NOT BE ALTERED.
C
C     WRITTEN BY RON JONES, WITH SLATEC COMMON MATH LIBRARY SUBCOMMITTEE
C     LATEST REVISION ---  19 MAR 1980
C
C***REFERENCES
C   JONES R.E., *SLATEC COMMON MATHEMATICAL LIBRARY ERROR HANDLING
C    PACKAGE*, SAND78-1189, SANDIA LABORATORIES, 1978.
C***ROUTINES CALLED  XGETUA,S88FMT,I1MACH
C***END PROLOGUE  XERSAV
C
      implicit double precision (a-h,o-z)
      INTEGER F(17),LUN(5)
      DIMENSION MESSG(1)
      DIMENSION MESTAB(10),NERTAB(10),LEVTAB(10),KOUNT(10)
      DATA MESTAB(1),MESTAB(2),MESTAB(3),MESTAB(4),MESTAB(5),
     1     MESTAB(6),MESTAB(7),MESTAB(8),MESTAB(9),MESTAB(10)
     2     /0,0,0,0,0,0,0,0,0,0/
      DATA NERTAB(1),NERTAB(2),NERTAB(3),NERTAB(4),NERTAB(5),
     1     NERTAB(6),NERTAB(7),NERTAB(8),NERTAB(9),NERTAB(10)
     2     /0,0,0,0,0,0,0,0,0,0/
      DATA LEVTAB(1),LEVTAB(2),LEVTAB(3),LEVTAB(4),LEVTAB(5),
     1     LEVTAB(6),LEVTAB(7),LEVTAB(8),LEVTAB(9),LEVTAB(10)
     2     /0,0,0,0,0,0,0,0,0,0/
      DATA KOUNT(1),KOUNT(2),KOUNT(3),KOUNT(4),KOUNT(5),
     1     KOUNT(6),KOUNT(7),KOUNT(8),KOUNT(9),KOUNT(10)
     2     /0,0,0,0,0,0,0,0,0,0/
      DATA KOUNTX/0/
      DATA F(1),F(2),F(3),F(4),F(5),F(6),F(7),F(8),F(9),F(10),
     1     F(11),F(12),F(13),F(14),F(15),F(16),F(17)
     2     /1H( ,1H1 ,1HX ,1H, ,1HA ,1H  ,1H  ,1H, ,1HI ,1H  ,
     3      1H  ,1H, ,1H2 ,1HI ,1H1 ,1H0 ,1H) /
C***FIRST EXECUTABLE STATEMENT  XERSAV
      IF (NMESSG.GT.0) GO TO 80
C     DUMP THE TABLE
         IF (KOUNT(1).EQ.0) RETURN
C        PREPARE FORMAT
         NCHAR = I1MACH(6)
         CALL S88FMT(2,NCHAR,F(6))
         NCOL = 20 - NCHAR
         CALL S88FMT(2,NCOL,F(10))
C        PRINT TO EACH UNIT
         CALL XGETUA(LUN,NUNIT)
         DO 60 KUNIT=1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
C           PRINT TABLE HEADER
            WRITE (IUNIT,10)
   10       FORMAT (32H0          ERROR MESSAGE SUMMARY/
     1              41H FIRST WORD      NERR     LEVEL     COUNT)
C           PRINT BODY OF TABLE
            DO 20 I=1,10
               IF (KOUNT(I).EQ.0) GO TO 30
               WRITE (IUNIT,F) MESTAB(I),NERTAB(I),LEVTAB(I),KOUNT(I)
   20       CONTINUE
   30       CONTINUE
C           PRINT NUMBER OF OTHER ERRORS
            IF (KOUNTX.NE.0) WRITE (IUNIT,40) KOUNTX
   40       FORMAT (41H0OTHER ERRORS NOT INDIVIDUALLY TABULATED=,I10)
            WRITE (IUNIT,50)
   50       FORMAT (1X)
   60    CONTINUE
         IF (NMESSG.LT.0) RETURN
C        CLEAR THE ERROR TABLES
         DO 70 I=1,10
   70       KOUNT(I) = 0
         KOUNTX = 0
         RETURN
   80 CONTINUE
C     PROCESS A MESSAGE...
C     SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
C     OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
      DO 90 I=1,10
         II = I
         IF (KOUNT(I).EQ.0) GO TO 110
         IF (MESSG(1).NE.MESTAB(I)) GO TO 90
         IF (NERR.NE.NERTAB(I)) GO TO 90
         IF (LEVEL.NE.LEVTAB(I)) GO TO 90
         GO TO 100
   90 CONTINUE
C     THREE POSSIBLE CASES...
C     TABLE IS FULL
         KOUNTX = KOUNTX+1
         ICOUNT = 1
         RETURN
C     MESSAGE FOUND IN TABLE
  100    KOUNT(II) = KOUNT(II) + 1
         ICOUNT = KOUNT(II)
         RETURN
C     EMPTY SLOT FOUND FOR NEW MESSAGE
  110    MESTAB(II) = MESSG(1)
         NERTAB(II) = NERR
         LEVTAB(II) = LEVEL
         KOUNT(II)  = 1
         ICOUNT = 1
         RETURN
      END
