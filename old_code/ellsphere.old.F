c
c $Id: ellsphere.F,v 1.37 2004/02/12 18:58:25 becker Exp becker $
c
c
c     runge kutta integration driver, calls odeintell
c     this routine also has ellderivs which determines velocities and their
c     derivatives. ellderivs uses routines in vel_interpol.F
c
c     this routine assumes that velocities are given as 
c     phi = 0. ... (nphi-1)*dphi
c     theta = dtheta/2 ... pi-dtheta/2
c
c
      subroutine ellsphere(vr,vtheta,vphi,nr,ntheta,nphi,dtheta,dphi,
     &     rlevels,vtimes,nvtimes,
     &     work,ti,tf,dt,ifellipse,eps,iflyap,renorm,rlyap,
     &     rlbailout,tbailout,ifremovestrain,vhdivmax)
c_______________________________________________________________________
c
c     ti  -> initial time 
c     tf  -> final time   both of these have to be in the right units
c                         for the time interpolation scheme to work
c     dt -> suggested first time step
c     work -> array with r,theta,phi coordinates and F components, if 
c     ifellipse is set, ie. work is work(3) or work(12)
c     vtimes(nvtimes*2) times at which velocity arrays are given
c     iflyap -> =0 - do not calculate lyapunov exponents
c               =1 calculate lyapunov exponents
c     renorm -> renormalization condition for deformation matrix
c            (if def(i,j)>renorm, renormalize)
c     rlyap(i) -> ith lyapunov exponent of tracer (should be 
c            dimensionalized to 3)
c     rlbailout -> bailout criterion for changes in lyapunov exponents
c     tbailout -> output, time of bailout on exit
c     ifremovestrain: remove the symmtetric part of the velocity gradient
c                    matrix depending on the rules laid out in ellderiv
c     vhdivmax(nvtimes): 
c               max horizontal divergence at surface, needed in ellderiv
c               for the strain removal to decide at what depth the critical 
c               lies, given at each time like vtimes
c
c_______________________________________________________________________

#include "fortran_precision.h"
      dimension vr(*), vtheta(*), vphi(*), rlevels(*),rlyap(3),
     &     alpha(3),vtimes(*),vhdivmax(*)
      parameter (smaxstep = 1.0)
      dimension work(12)

      !set parameters and call odeint
      if (ifellipse.eq.0) then
c no strain tracing and vel gradient matrix needed
         nwork=3 
         if (iflyap.ne.0) then
            write (6,*) 'You cannot calculate lyapunov exponents'
            write (6,*) 'without calculating strain (ifellipse)!'
            write (6,*) 'ifellipse,iflyap',ifellipse,iflyap
            stop
         end if
      else
         nwork=12
      end if
      
      t1 = ti   ! start time
      t2 = tf   ! final time
      step = dt ! time step
c
c     max step is internally limited to smaxstep
c
      if(abs(step) .gt. smaxstep)then
         step=(abs(step)/step) * smaxstep
      endif
      stepmin=0.d0
c
c     integrate to obtain position and deformation matrix at t2
c
      call odeintell(work,nwork,t1,t2,eps,step,stepmin,
     &     vr,vtheta,vphi,nr,ntheta,nphi,dtheta,dphi,rlevels,vtimes,
     &     nvtimes,iflyap,renorm,rlyap,rlbailout,tbailout,
     &     ifremovestrain,vhdivmax)

      if (iflyap.ne.0) then
      !increment lyapunov exponents
         call gramschmidt(work,alpha)
         !and increment lyapunov exponents
         do i=1,3
            rlyap(i)=rlyap(i)+log(alpha(i))
         end do
      end if

c     tracer should be within radial bounds since this is corrected for in  rkqsell
      if(XTRACER(1).gt.1.0)then
         print *,'ellsphere: tracer x_r > 1: ',XTRACER(1)
         stop
      endif

!     adjust the location vector, work(1...3)
      call check_phys_lim_tracer(XTRACER(1))
      return
      end

