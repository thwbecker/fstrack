c
c aniosotropic reflectivity package from Vera Schulte-Pelkum in March 2005
c based on code by Peter Shearer and various others, see individual READMEs 
c and source code comments for details and description
c  
c minor modifications by Thorsten Becker
c
c $Id: getspec.f,v 1.1 2005/03/23 20:29:35 becker Exp $
c
c
cc GETSPEC gets the spectrum from a time series.
c
c  Inputs:   ts    =  real time series
c            ntime =  number of time points (must be power of 2)
c            dt    =  time series spacing
c            tap   =  % cosine taper
c                  =  1. for Hanning taper
c                  =  0. for no taper
c  Returns:  spec  =  complex spectrum (nfreq = 1+ntime/2)
c                     spec(1) for zero frequency
c                     spec(nfreq) for Nyquist
c            nfreq =  number of frequency points !*****added Oct. 94
c            df    =  frequency spacing
c
c  Requires:  FFA, FFAFFSSUBS
c
      subroutine GETSPEC(ts_in,ntime,dt,tap,spec_out,nfreq,df)
      real ts(32768),ts_in(32768)
      complex spec(16385),spec_out(16385)
      equivalence (spec(1),ts(1))

      if (ntime.gt.32768) then
         print *,'***ERROR in GETSPEC: too many time points'
         print *,'ntime = ',ntime
         stop
      end if

c
      call TAPER1(ts_in,ntime,tap,ts)
c
      call FFA(ts,ntime,.false.,.false.)
c   
      df=1./(dt*float(ntime))
      nfreq=1+ntime/2
c
      do 100 i=1,nfreq-1
c         spec_out(i)=spec(i)/(float(ntime)*df)
         spec_out(i)=spec(i)*dt
100   continue
      spec_out(nfreq)=cmplx(aimag(spec_out(1)),0.)
      spec_out(1)=real(spec_out(1))
c
      return
      end
c
c
c TAPER1 multiplies a time series by a cosine squared taper
c
c Inputs:      x  =  input time series
c              n  =  number of time points
c           frac  =  % cosine squared taper
c                 =  1. for Hanning taper
c                 =  0. for no taper
c Returns:  xtap  =  tapered time series
c
      subroutine TAPER1(x,n,frac,xtap)
      real x(*),xtap(*)
      pihalf=3.1415927/2.
      do 10 i=1,n
10       xtap(i)=x(i)
      if (frac.eq.0.) go to 999
      xlen=float(n-1)
      taplen=frac*xlen/2.
      i2=nint(taplen)
      do 20 i=1,i2
         angle=pihalf*float(i-1)/(i2-1)
         fact=sin(angle)**2
         xtap(i)=x(i)*fact
         xtap(n+1-i)=x(n+1-i)*fact
20    continue
999   return
      end
