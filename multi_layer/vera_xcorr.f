c
c aniosotropic reflectivity package from Vera Schulte-Pelkum in March 2005
c based on code by Peter Shearer and various others, see individual READMEs 
c and source code comments for details and description
c  
c minor modifications by Thorsten Becker
c
c $Id: vera_xcorr.f,v 1.1 2005/11/19 01:12:53 becker Exp becker $
c
c
cc Date: Tue, 4 Jan 2005 10:38:46 -0800
c From: Peter Shearer <pshearer@ucsd.edu>
c Here is a routine.  It assumes windows within two longer time series. It shifts
c both windows in a symmetric fashion when doing the xcor.  You can make it run 
c faster by setting ndeci to numbers larger than one, at the cost of a less
c accurate result.

c GETXCOR computes cross-correlations
c
c Input: a1  =  time series 1
c        a2  =  time series 2
c        n1  =  number of points in a1
c        n2  =  number of points in a2
c        iref1 = reference window starting point for 1
c        iref2 = reference window starting point for 2
c        nxcor = number of points to cross-correlate (i.e., window length)
c        nshift = maximum number of points to shift
c        ndeci =  decimation factor for faster approximate solution
c                 (1=complete, 2=every other point, etc.)
c
c Returns: xcor = cross-correlation function (-1000:1000)
c          rmax = maximum correlation coef.
c          ir1 = shift (number of points) at rmax
c          d1  = best-fitting scale factor a2/a1 at rmax
c          rmin = minimum correlation coef.
c          ir2 = shift (number of points) at rmin
c          d2  = best-fitting scale factor a2/a1 at rmin
c
c Note:  ir1,ir2 > 0 when a2 pulse is later than a1 pulse
c        d1,d2 > 1 when a2 pulse ampl. is greater than a1 pulse ampl.  

      subroutine vera_getxcor(a1,a2,n1,n2,iref1,iref2,nxcor,nshift,
     &        ndeci,xcor,rmax,ir1,d1,rmin,ir2,d2)
      integer nsmax
      parameter (nsmax=1000)
      real a1(n1),a2(n2),xcor(-nsmax:nsmax)

      if(nshift.gt.nsmax)then
         print *,'vera_getxcor: out of bounds, nshift',nshift,nsmax
         stop
      endif

      ioffmax=nshift/2+1

      if(iref1-ioffmax.lt.1)then
c         print *,ioffmax
         ioffmax = iref1-1
c         print *,ioffmax
      endif

      if(iref2-ioffmax.lt.1)then
         ioffmax = iref2-1
      endif

      if(ioffmax.lt.1)then
         print *,'error ioffmax',ioffmax
         stop
      endif

c automatically reduce offset
      if(iref1+nxcor+ioffmax.gt.n1)then
         nxcor = n1-iref1-ioffmax
      endif
      if(iref2+nxcor+ioffmax.gt.n2)then
         nxcor = n2-iref2-ioffmax
      endif
      

c automatically reduce number of points to shift


      if ((iref1-ioffmax.lt.1).or.(iref2-ioffmax.lt.1).or.
     &     (iref1+nxcor+ioffmax.gt.n1).or.
     &     (iref2+nxcor+ioffmax.gt.n2)) then
         print *,'getxcor: WARNING: outside data array in GETXCOR'
         print *,iref1-ioffmax
         print *,iref2-ioffmax
         print *,iref1+nxcor+ioffmax,n1
         print *,iref2+nxcor+ioffmax,n2
         rmax=0.
         rmin=0.
         ir1=0
         ir2=0
         d1=1.
         d2=1.
         return
      end if

      rmax=-1.
      rmin=1.
      do 100 ishift=-nshift,nshift
         ishifta=ishift/2
         ishiftb=ishift-ishifta
         suma=0.
         sumb=0.
         sumab=0.
         do ixcor=0,nxcor,ndeci
            ita=iref1+ixcor-ishifta
            itb=iref2+ixcor+ishiftb
            suma=suma+a1(ita)**2
            sumb=sumb+a2(itb)**2
            sumab=sumab+a1(ita)*a2(itb)
         enddo
         denom=suma*sumb
         sdenom=sqrt(denom)
         if (sdenom.ne.0.) then
            test=sumab/sdenom
         else
            xcor(ishift)=0.
            d1=1.
            d2=1.
            go to 100
         end if
         xcor(ishift)=test
         if (test.gt.rmax) then
            rmax=test
            ir1=ishift
            d1=sumb/suma
         end if
         if (test.lt.rmin) then
            rmin=test
            ir2=ishift
            d2=sumb/suma
         end if
 100  continue
      d1=sqrt(d1)
      d2=sqrt(d2)

      return
      end
