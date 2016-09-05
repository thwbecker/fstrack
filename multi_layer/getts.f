c
c aniosotropic reflectivity package from Vera Schulte-Pelkum in March 2005
c based on code by Peter Shearer and various others, see individual READMEs 
c and source code comments for details and description
c  
c minor modifications by Thorsten Becker
c
c $Id: getts.f,v 1.1 2005/03/23 20:29:35 becker Exp $
c
c
cc GETTS gets the time series from a complex spectrum using a cosine
c squared taper.
c
c  Inputs:   spec   =  complex spectrum
c                      spec(1) for zero frequency
c                      spec(nfreq) is Nyquist
c            nfreq  =  number of frequency points (must be power of 2 + 1)
c            df     =  frequency spacing
c            ntap   =  0  don't use taper
c                   =  1  apply cosine squared taper
c  Returns:  ts     =  real time series
c            ntpts  =  number of time points          !***added Oct. 94
c            dt     =  time series spacing
c
c  Requires:  FFS, FFAFFSSUBS
c
      subroutine GETTS(specin,nfreq,df,ntap,tsout,ntpts,dt)
c
      real filcof(16385),tsout(32768),ts(32768)
      complex specin(16385),spec(16385)
      equivalence (spec(1),ts(1))

      if (nfreq.gt.16385) then
         print *,'***ERROR in GETTS: too many freq points'
         print *,'nfreq = ',nfreq
         stop
      end if

      do 10 i=1,nfreq
         spec(i)=specin(i)
         filcof(i)=1.
10    continue
c
      fmax=float(nfreq-1)*df   !this is nyquist frequency
      ntpts=(nfreq-1)*2
      dt=1./(df*float(ntpts))
c
      fl1=0.
      fl2=0.
      fh1=0.
      fh2=fmax
C
C  THE SUBROUTINE TAPER RETURNS A REAL ARRAY NF POINTS LONG
C  CONTAINING THE COSINE SQUARED TAPER REQUESTED.THE LAST
C  LOOP CONTAINS THE FT NORMALIZATION.
C
      NF=nfreq
      if (ntap.eq.1) then
           CALL TAPER(FH1,FH2,FL1,FL2,FILCOF,DT,NF)
      end if
           DO 7362 N3=1,NF
7362       FILCOF(N3)=DF*FILCOF(N3)
C
C  THE ZERO FREQUENCY COMPONENT OF SPCTRA IS SET UP
C  WITH THE DC COMPONENT IN THE REAL PART AND THE
C  NYQUIST FREQUENCY COMPONENT (IDENTICALLY ZERO) IN THE
C  IMAGINARY PART. NOTE THAT WE HAVE ASSUMED THE FIRST
C  FREQUENCY AT WHICH THE RESPONSE WAS CALCULATED IN
C  PROSE WAS ZERO HZ.
C
      SPEC(1)=CMPLX(REAL(SPEC(1)),0.)
C
C  MULTIPLY BY TAPER
C
      DO 301 L=1,nfreq
         SPEC(L)=SPEC(L)*FILCOF(L)
301   CONTINUE
c
      spec(1)=cmplx(real(spec(1)),real(spec(nfreq)))
C
C  DO INVERSE FOURIER TRANSFORM USING FFS. Normalization done in
c  taper.
C
      CALL FFS(SPEC(1),NTPTS,.false.,.false.)
c
      do 500 i=1,ntpts
         tsout(i)=ts(i)
500   continue
C
      return
      end
C
C
      SUBROUTINE TAPER(FH1,FH2,FL1,FL2,FC,DT,NF)
C  THIS SUBROUTINE COMPUTES A COSINE SQUARED TAPER BETWEEN
C  FH1 AND FH2 AND ONE BETWEEN FL1 AND FL2. IT PLACES THE
C  RESULT IN FC.
      DIMENSION FC(1)
      PI=3.1415927
      TIME=(NF-1)*DT*2.
      DF=1./TIME
      ILL=FL1*TIME+1
      IRL=FL2*TIME+1
      ILH=FH1*TIME+1
      IRH=FH2*TIME+1
      IF (IRH.NE.ILH) DFH=PI/FLOAT(IRH-ILH)
      IF (IRL.NE.ILL) DFL=PI/FLOAT(IRL-ILL)
      DO 100 I=1,NF
         IF (I.LE.ILH.AND.I.GE.IRL) THEN
            FC(I)=1.0
            GOTO 100
         ENDIF
         IF (I.GT.IRH.OR.I.LT.ILL) THEN
            FC(I)=0.0
            GOTO 100
         ENDIF
         IF (I.GT.ILH) THEN
            FC(I)=(COS(DFH*(I-ILH))/2.+0.5)**2.
         ELSE
            FC(I)=(COS(DFL*(IRL-I))/2.+0.5)**2.
         ENDIF
100   CONTINUE
      RETURN
      END
C
C
