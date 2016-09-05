c
c
c     cartesian version
c
c     Runge Kutta driver modified from Numerical Recipes
c 
c     $Id: odeintell_cart.F,v 1.3 2011/04/12 06:19:09 becker Exp $
c
      SUBROUTINE odeintell_cart(ystart,nvar,x1,x2,eps,h1,hmin,
     &     vx,vy,vz,nx,ny,nz,dx,dy,xymax,zlevels,vtimes,
     &     nvtimes,iflyap,renorm,rlyap,rlbailout,tbailout,
     &     ifremovestrain,vhdivmax) 
#include "fortran_precision.h"
       dimension vz(*),vy(*),vx(*), zlevels(nz), rlyap(3),xymax(2),
     &      alpha(3) ! original lengths of orthonormalized vectors
       dimension vtimes(*),vhdivmax(*)
       INTEGER nbad,nok,nvar,KMAXX,MAXSTP,NMAX,nvtimes
       COMP_PRECISION TINY
       PARAMETER (MAXSTP=80000,NMAX=12,KMAXX=200,TINY=1.e-30)
       INTEGER i,kmax,kount,nstp
       COMP_PRECISION dxsav,h,hdid,hnext,x,xsav,dwork,xp,
     &      work,yp,yscal,dx,dy,eps,h1,hmin,x1,x2, ystart,
     &      rlbailout,tlalpha,tlalphaold(3),vtimes,xymax
       dimension dwork(NMAX),xp(KMAXX),work(NMAX),yp(NMAX,KMAXX),
     &      yscal(NMAX),ystart(nvar)
       COMMON /path/ kmax,kount,dxsav,xp,yp
       
       logical ifrenorm

       icntbailout=0
       xsav=0.
       
#ifdef FSTRACK_DEBUG
       if(nvar.gt.NMAX)then
          print *,'odeintell_cart: nvar > NMAX, nvar:',nvar
          stop
       endif
       if((iflyap.eq.1).and.(rlbailout.le.0.0))then
          print *,'odeintell_cart: error, rlbailout=',rlbailout
          stop
       endif
#endif
       if(iflyap.ne.0)then
c     init for lyapunov bailout crit
          do i=1,3
             tlalphaold(i) = -1.0d0
          enddo
       endif
c
       x=x1
       h=sign(h1,x2-x1)
       nok=0
       nbad=0
       kount=0
       do i=1,nvar
          work(i)=ystart(i)
       enddo
       if (kmax.gt.0) xsav=x-2.*dxsav
       do nstp=1, MAXSTP
c     get derivatives
          call ellderivs_cart(x,work,dwork,nvar,
     &         vx,vy,vz,nx,ny,nz,dx,dy,xymax,zlevels,
     &         vtimes,nvtimes,ifremovestrain,vhdivmax)
          yscalmax=TINY
          do i=1,nvar
             yscal(i)=abs(work(i))+abs(h*dwork(i))+TINY
             yscalmax=max(yscalmax,yscal(i))
          enddo
          if (nvar.gt.3) then
c     rescale yscal such that all components of F have the
c     same accuracy
             do i=3,nvar
                yscal(i)=yscalmax
             end do
          end if
c          
          if(kmax.gt.0)then
             if(abs(x-xsav).gt.abs(dxsav)) then
                if(kount.lt.kmax-1)then
                   kount=kount+1
                   xp(kount)=x
                   do i=1,nvar
                      yp(i,kount)=work(i)
                   enddo
                   xsav=x
                endif
             endif
          endif
          if((x+h-x2)*(x+h-x1).gt.0.) h=x2-x
c
c     move by one quality controlled RK step
          call rkqsell_cart(work,dwork,nvar,x,h,eps,yscal,hdid,hnext,
     &         vx,vy,vz,nx,ny,nz,dx,dy,xymax,zlevels,
     &         vtimes,nvtimes,ifremovestrain,vhdivmax)
c     
c     if necessary, check to see if we need to renormalize
          if (iflyap.ne.0) then
             ifrenorm=.false.
             do i=4, nvar
                if(abs(work(i)).gt.renorm) ifrenorm=.true.
             end do
             if (ifrenorm) then
                icntbailout=0
c     renormalize
                call gramschmidt(work,alpha)
c     and increment lyapunov exponents
                do i=1,3
                   slalpha = log(alpha(i))
c     check for bailout
                   if(rlyap(i).ne.0.0)then
                      tlalpha = abs(slalpha/rlyap(i))
c     need at least one renorm step (tlalphaold() init as -1)
c     and a decrease in the Delta addition
                      if((tlalphaold(i).gt.tlalpha).and. 
c     will also not bailout if time < 50 Myr, use abs(x) since we might go backward
     &                     (abs(x).gt.50.0).and.
     &                     (tlalpha.lt.rlbailout))then
                         icntbailout = icntbailout + 1
                      endif
c     save increment for next check
                      tlalphaold(i) = tlalpha
                   endif
                   rlyap(i)=rlyap(i) + slalpha
                end do
             end if
             if(icntbailout.eq.3)then 
c     only bailout if all exponents < threshold
c     set time of bailout
                tbailout = x
                return 
             endif
          end if
          if(hdid.eq.h)then
             nok=nok+1
          else
             nbad=nbad+1
          endif
          if((x-x2)*(x2-x1).ge.0.)then
             do i=1,nvar
                ystart(i)=work(i)
             enddo
             if(kmax.ne.0)then
                kount=kount+1
                xp(kount)=x
                do i=1,nvar
                   yp(i,kount)=work(i)
                enddo
             endif
             tbailout = x
             return
          endif
          if(abs(hnext).lt.hmin) then
             print *,'odeintell_cart: steps <  minimum in odeint'
             stop
          endif
          h=hnext
       enddo
       print *,'odeintell_cart: too many steps in odeint'
       stop
       return
       END
c
c adaptive RK routine, same as in odeintell except for routine names
c

      SUBROUTINE rkqsell_cart(work,dwork,n,x,htry,eps,yscal,hdid,hnext,
     &     vx,vy,vz,nx,ny,nz,dx,dy,xymax,zlevels,
     &     vtimes,nvtimes,ifremovestrain,vhdivmax)
#include "fortran_precision.h"
      dimension vx(*), vy(*), vz(*), zlevels(*),xymax(2),
     &     vtimes(*),vhdivmax(*)
      INTEGER n,NMAX,nvtimes
      COMP_PRECISION eps,hdid,hnext,htry,x,dwork(n),work(n),yscal(n),
     &     dy, dx,vtimes,xymax
      PARAMETER (NMAX=12)
      INTEGER i
      COMP_PRECISION errmax,h,xnew,yerr(NMAX),ytemp(NMAX),SAFETY,
     &     PGROW,PSHRNK,ERRCON
      PARAMETER (SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.89e-4)

      h=htry
c
c     try forward advection
c
 1    call rkckell_cart(work,dwork,n,x,h,ytemp,yerr,
     &     vx,vy,vz,nx,ny,nz,dx,dy,xymax,zlevels,
     &     vtimes,nvtimes,ifremovestrain,vhdivmax)
c
c     init error conditions
c
      errmax=0.d0
      do i=1,n
c     
c determine error condition
c     
c     errmax=max(errmax,abs(yerr(i)/yscal(i)))
         tmp=abs(yerr(i)/yscal(i))
         if (tmp.gt.errmax) then
            errmax=tmp
            imax=i
         end if
      enddo
      errmax = errmax/eps
      if(errmax.gt.1.0)then ! error too large
         htemp=SAFETY*h*(errmax**PSHRNK)
         h=sign(max(abs(htemp),0.1*abs(h)),h)
         xnew=x+h
         if(xnew.eq.x)then
            print *,'rkqsell_cart: stepsize underflow: time:',x
            print *,'rkqsell_cart: stepsize underflow: x:',
     &           work(1),work(2),work(3)
            stop
         endif
         goto 1
      else
         if(errmax.gt.ERRCON)then
            hnext=SAFETY*h*(errmax**PGROW)
         else
            hnext=5.*h
         endif
         hdid=h
         x=x+h
         do i=1,n
            work(i)=ytemp(i)
         enddo
         return
      endif
      END
c
c     take RK step
c
      SUBROUTINE rkckell_cart(y,dydx,n,x,h,yout,yerr,
     &     vx,vy,vz,nx,ny,nz,dx,dy,xymax,zlevels,
     &     vtimes,nvtimes,ifremovestrain,vhdivmax)
#include "fortran_precision.h"
      dimension vz(*),vy(*), vx(*), zlevels(nz),xymax(2)
      INTEGER n,NMAX,nvtimes
      COMP_PRECISION h,x,dydx(n),y(n),yerr(n),yout(n),dx,dy,
     *     vtimes(*),vhdivmax(*),xymax,vx,vy,vz
      PARAMETER (NMAX=12)
C     USES derivs
      INTEGER i
      COMP_PRECISION ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),
     &     ak6(NMAX),
     &     ytemp(NMAX),A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,
     &     B53,B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
      PARAMETER (A2=.2d0,A3=.3d0,A4=.6d0,A5=1.d0,A6=.875d0,B21=.2d0,
     &     B31=3.d0/40.d0,B32=9.d0/40.d0,B41=.3d0,B42=-.9d0,B43=1.2d0,
     &     B51=-11.d0/54.d0,B52=2.5d0,B53=-70.d0/27.d0,B54=35.d0/27.d0,
     &     B61=1631.d0/55296.d0,B62=175.d0/512.d0,B63=575.d0/13824.d0,
     &     B64=44275.d0/110592.d0,B65=253.d0/4096.d0,C1=37.d0/378.d0,
     &     C3=250.d0/621.d0,C4=125.d0/594.d0,C6=512.d0/1771.d0,
     &     DC1=C1-2825.d0/27648.d0,DC3=C3-18575.d0/48384.d0,
     &     DC4=C4-13525.d0/55296.d0,DC5=-277.d0/14336.d0,DC6=C6-.25d0)
      do  i=1,n
         ytemp(i)=y(i)+B21*h*dydx(i)
      enddo
      call ellderivs_cart(x+A2*h,ytemp,ak2,n,
     &     vx,vy,vz,nx,ny,nz,dx,dy,xymax,zlevels,
     &     vtimes,nvtimes,ifremovestrain,vhdivmax) 
      do i=1,n
         ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
      enddo
      call ellderivs_cart(x+A3*h,ytemp,ak3,n,
     &     vx,vy,vz,nx,ny,nz,dx,dy,xymax,zlevels,
     &     vtimes,nvtimes,ifremovestrain,vhdivmax) 
      do i=1,n
         ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
      enddo
      call ellderivs_cart(x+A4*h,ytemp,ak4,n,
     &     vx,vy,vz,nx,ny,nz,dx,dy,xymax,zlevels,
     &     vtimes,nvtimes,ifremovestrain,vhdivmax) 
      do i=1,n
         ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
      enddo
      call ellderivs_cart(x+A5*h,ytemp,ak5,n,
     &     vx,vy,vz,nx,ny,nz,dx,dy,xymax,zlevels,
     &     vtimes,nvtimes,ifremovestrain,vhdivmax) 
      do i=1,n
         ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+
     &        B65*ak5(i))
      enddo
      call ellderivs_cart(x+A6*h,ytemp,ak6,n,
     &     vx,vy,vz,nx,ny,nz,dx,dy,xymax,zlevels,
     &     vtimes,nvtimes,ifremovestrain,vhdivmax) 
      do i=1,n
         yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
      enddo
      do i=1,n
         yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*
     &        ak6(i))
      enddo
      return
      END




