c
c $Id: ellsphere.F,v 1.16 2001/10/31 18:34:29 kellogg Exp $
c
      subroutine ellsphere(vr,vtheta,vphi,nr,ntheta,nphi,dtheta,dphi,
     &     rlevels,work,dt,ifellipse,eps,iflyap,renorm,rlyap)
c_______________________________________________________________________
c
c     iflyap -> =0 - do not calculate lyapunov exponents
c              =1 calculate lyapunov exponents
c     renorm -> renormalization condition for deformation matrix
c            (if def(i,j)>renorm, renormalize)
c     rlyap(i) -> ith lyapunov exponent of tracer (should be 
c            dimensionalized to 3)
c_______________________________________________________________________

#include "precision.h"
      dimension vr(*), vtheta(*), vphi(*), rlevels(*),rlyap(3)
     &     , alpha(3)
      parameter (smaxstep = 0.5)
      dimension work(12)

      !set parameters and call odeint
      if (ifellipse.eq.0) then
c no strain tracing and vel gradient matrix needed
         nwork=3 
         if (iflyap.ne.0) then
            write (6,*) 'You cannot calculate lyapunov exponents'
            write (6,*) 'without calculating ellipses!!'
            write (6,*) 'ifellipse,iflyap',ifellipse,iflyap
            stop
         end if
      else
         nwork=12
      end if

      t1=0.d0
      t2=dt
c initial step size
      step=dt/30.d0
      if(abs(step) .gt. smaxstep)then
         step=(abs(step)/step) * smaxstep
      endif
      stepmin=0.d0
c     integrate position and deformation matrix
      call odeintell(work,nwork,t1,t2,eps,step,stepmin,
     *     vr,vtheta,vphi,nr,ntheta,nphi,dtheta,dphi,rlevels,
     &     iflyap,renorm,rlyap)

      if (iflyap.ne.0) then
      !increment lyapunov exponents
         call gramschmidt(work,alpha)
         !and increment lyapunov exponents
         do i=1,3
            rlyap(i)=rlyap(i)+log(alpha(i))
         end do
      end if
         
#ifdef DEBUG
      detf=DEF1(1)*DEF2(2)*DEF3(3)+
     &     DEF1(2)*DEF2(3)*DEF3(1)+
     &     DEF1(3)*DEF2(1)*DEF3(2)-
     &     DEF1(1)*DEF2(3)*DEF3(2)-
     &     DEF1(2)*DEF2(1)*DEF3(3)-
     &     DEF1(3)*DEF2(2)*DEF3(1)
      if(abs(detf-1.0).gt.100.*eps)then
         print *,'detf=',detf
      endif
#endif
      !adjust position values
      !0<=theta<=pi
      do while (XTRACER(2).gt.twopi)
         XTRACER(2)=XTRACER(2)-twopi
      end do
      do while (XTRACER(2).lt.0.d0)
         XTRACER(2)=XTRACER(2)+twopi
      end do
      if (XTRACER(2).gt.pi) then
         XTRACER(2)=twopi-XTRACER(2)
         XTRACER(3)=XTRACER(3)+pi
      end if

      !0<=phi<=2*pi
      do while (XTRACER(3).gt.twopi)
         XTRACER(3)=XTRACER(3)-twopi
      end do
      do while (XTRACER(3).lt.0.d0)
         XTRACER(3)=XTRACER(3)+twopi
      end do

      return
      end



      subroutine ellderivs(t,work,workder,nwork,
     &     vr,vtheta,vphi,nr,ntheta,nphi,dtheta,dphi,rlevels) 
c_______________________________________________________________________
c
c     Note:  passed muti-dimensional arrays are in c-order
c
c_______________________________________________________________________

#include "precision.h"
      dimension vr(*), vtheta(*), vphi(*), rlevels(nr)      
      dimension work(12),workder(12)
      dimension grid(0:3,3),weight(0:3,0:1,3)
      dimension igrid(0:3,3),ider(0:3), dvr(0:3),dvtheta(0:3),
     &     dvphi(0:3),vgm(3,3)
c     debugging
      dimension tmat(3,3),eval(3),evec(3,3)
      
      logical icalc_vgm

c     decide if we need to calculate the velocity gradient matrix
      if(nwork.gt.3)then 
                                ! for vel diff loop
         icalc_vgm=.true. 
         idindex=3 
      else 
                                ! only velocities
         icalc_vgm=.false.
         idindex=0 
      endif

      !pick cell for tracer
      !first, radial component
      ixtracer=0
      i=0
      do while ((ixtracer.eq.0).and.(i.lt.nr))
         i=i+1
         if (XTRACER(1).le.rlevels(i)) ixtracer=i
      end do
      if(ixtracer .eq. 0)then 
         ixtracer=nr
         print *,'out of box: t:',t,' r:',work(1)
         stop
      endif
      

      !pick grid points
      do i=0,3
         igrid(i,1)=ixtracer-1+i
      end do

      !make sure all grid points exist
      do while (igrid(0,1).le.0)
         do i=0,3
            igrid(i,1)=igrid(i,1)+1
         end do
      end do
      do while (igrid(3,1).gt.nr)
         do i=0,3
            igrid(i,1)=igrid(i,1)-1
         end do
      end do

      !find values of r for each grid point
      do i=0,3
         grid(i,1)=rlevels(igrid(i,1))
      end do


      !calculate 1-d weights for each grid point
      ngrid=3
      iorder=1
      call weights(XTRACER(1),grid(0,1),ngrid,iorder,weight(0,0,1))

      !now, do theta component
      do while (XTRACER(2).gt.twopi)
         XTRACER(2)=XTRACER(2)-twopi
      end do
      do while (XTRACER(2).lt.0.d0)
         XTRACER(2)=XTRACER(2)+twopi
      end do
      if (XTRACER(2).gt.pi) then
         XTRACER(2)=twopi-XTRACER(2)
         XTRACER(3)=XTRACER(3)+pi
      end if

      ixtracer=nint(XTRACER(2)/dtheta)
      
      !pick grid points
      do i=0,3
         igrid(i,2)=ixtracer-1+i
      end do

      !adust grid points to avoid wrap-around
      do while (igrid(0,2).le.0)
         do i=0,3
            igrid(i,2)=igrid(i,2)+1
         end do
      end do
      do while (igrid(3,2).gt.ntheta)
         do i=0,3
            igrid(i,2)=igrid(i,2)-1
          end do
      end do

      !find values of theta 
      do i=0,3
         grid(i,2)=(igrid(i,2)-.5d0)*dtheta
      end do

      !calculate 1-d weights for each grid point
      call weights(XTRACER(2),grid(0,2),ngrid,iorder,weight(0,0,2))
      
      !and, finally, the phi component
      do while (XTRACER(3).gt.twopi)
         XTRACER(3)=XTRACER(3)-twopi
      end do
      do while (XTRACER(3).lt.0.d0)
         XTRACER(3)=XTRACER(3)+twopi
      end do

      ixtracer=int(XTRACER(3)/dphi)+1

      !pick grid points
      do i=0,3
         igrid(i,3)=ixtracer-1+i
      end do

      !find values of phi
      do i=0,3
         grid(i,3)=(igrid(i,3)-1)*dphi
      end do
      !calculate 1-d weights for each grid point
      call weights(XTRACER(3),grid(0,3),ngrid,iorder,weight(0,0,3))
      
      !now make sure grid points exist
      do i=0,3
         if (igrid(i,3).lt.1) igrid(i,3)=igrid(i,3)+nphi
         if (igrid(i,3).gt.nphi) igrid(i,3)=igrid(i,3)-nphi
      end do

      ntnp = ntheta * nphi

      !check if we are in polar circle
      if ((XTRACER(2).gt.dtheta).or.((pi-XTRACER(2)).gt.dtheta)) then
c
c     OK, ie. not in polar regions
c
c     first calculate velocities only (idindex=0) or vel and  derivatives 
c     of velocity (e.g. v_(r,r)) if needed (idindex=3)
         do m=0, idindex        !m=0 -> no derivative
            dvr(m)=0.d0         !m=1 -> derivative wrt r
            dvtheta(m)=0.d0     !m=2 -> derivative wrt theta
            dvphi(m)=0.d0       !m=3 -> derivative wrt phi
            do msub=0,3
               ider(msub)=0
            end do
            ider(m)=1 
            do i=0,3
               do j=0,3
                  do k=0,3
                     rnet=weight(i,ider(1),1)*
     &                    weight(j,ider(2),2)*
     &                    weight(k,ider(3),3)
                     index=(igrid(i,1)-1)*ntnp + (igrid(j,2)-1)*nphi +
     &                    igrid(k,3)
                     dvr(m)=    dvr(m)     + rnet*vr(index)
                     dvtheta(m)=dvtheta(m) + rnet*vtheta(index)
                     dvphi(m)=  dvphi(m)   + rnet*vphi(index)
                  end do
               end do
            end do         
         end do
         
      else    
c
c     polar regions
c     
         print *,'polar! theta(for example):',work(2)
c
         do m=0,idindex
            dvr(m)=0.d0      
            dvtheta(m)=0.d0  
            dvphi(m)=0.d0    
            do msub=0,3
               ider(msub)=0
            end do
            ider(m)=1 
            do i=0,3            ! radial node
               do k=0,3         !phi node
                  if (m.ne.2) then ! just use values for theta=dtheta
                     j=1
                     !check if n pole or s pole
                     if (pi-XTRACER(2).lt.dtheta) j=ntheta
                     index=(igrid(i,1)-1)*ntnp+(j-1)*nphi + 
     &                    igrid(k,3)
                     rnet=0.5*weight(i,ider(1),1)*weight(i,ider(3),3)
                     dvr(m)=dvr(m)+rnet*vr(index)
                     dvtheta(m)=dvtheta(m)+rnet*vtheta(index)
                     dvphi(m)=dvphi(m)+rnet*vphi(index)
                  else
                  !for theta derivative, project derivative in from 
                  ! two nodes closest to pole
                     do j=0,1
                        jindex=j
                        ishift=1 !see rnet calculation.  this gets the
                                 !sign of the weight right for the deriv
                        if (pi-XTRACER(2).lt.dtheta) then !s pole
                           jindex=ntheta-j
                           shift=0
                        end if
                        index=(igrid(i,1)-1)*ntnp+(jindex-1)*nphi + 
     &                       igrid(k,3)
                        rnet=((-1)**(j+ishift)/dtheta)*
     &                       weight(i,ider(1),1)*
     &                       weight(i,ider(3),3)
                        dvr(m)=dvr(m)+rnet*vr(index)
                        dvtheta(m)=dvtheta(m)+rnet*vtheta(index)
                        dvphi(m)=dvphi(m)+rnet*vphi(index)
                     end do
                  end if
               end do
            end do
         end do
      end if
         
c     factors
      finvr=1.d0/XTRACER(1)
      sintheta=sin(XTRACER(2))

c     repack velocities
      workder(1)=dvr(0)
      workder(2)=finvr*dvtheta(0)
      workder(3)=finvr*dvphi(0)/sintheta


      if(icalc_vgm)then
         costheta=cos(XTRACER(2))
         cottheta=1.d0/tan(XTRACER(2))
c
c     build velocity gradient matrix (transpose of grad v)
c     (see, e.g., Dahlen & Tromp, p. 836)
         
      
         vgm(1,1)=dvr(1)
         vgm(2,1)=dvtheta(1)
         vgm(3,1)=dvphi(1)   
         vgm(1,2)=finvr*(dvr(2)-dvtheta(0))
         vgm(2,2)=finvr*(dvtheta(2)+dvr(0))
         vgm(3,2)=finvr*dvphi(2)
         vgm(1,3)=finvr*(dvr(3)/sintheta-dvphi(0)) 
         vgm(2,3)=finvr*(dvtheta(3)/sintheta-cottheta*dvphi(0)) 
         vgm(3,3)=finvr*(dvphi(3)/sintheta+dvr(0)+cottheta*dvtheta(0))
c     calculate divergence and adjust vgm matrix so that the divergence
c     is zero (not ideal, but functional)
         div=(vgm(1,1)+vgm(2,2)+vgm(3,3))/3.0d0
         do i=1,3
            vgm(i,i)=vgm(i,i)-div
         end do
c     write divergence
c        write (6,*) 'ellderivs: trace vgm=', div
c
c     the time derivative of F is L dot F, repack into workder array
c

         index=3
         do i=1,3
            do j=1,3
               index=index+1
               workder(index)=vgm(i,1)*DEF1(j)+vgm(i,2)*DEF2(j)+
     &              vgm(i,3)*DEF3(j)
            end do
         end do
      endif

      return
      end





      subroutine weights(xi,x,n,m,c)
c*********************************************
c copied from Fornberg(1996),p.168
c calculates weights for 1-d interpolations
c INPUT PARAMETERS:
c xi: point at which approximations are to be accurate
c x : xcoords for grid points, array dimensioned to x(0:n)
c n : # of grid points minus 1
c m : highest order of derivative to be approximated
c
c OUTPUT PARAMETER:
c c : weights, array dimensioned  c(0:n, 0:m)
c     the element c(j,k) contains the weight to be applied
c     at x(j) when the kth derivative is approximated by
c     a stencil extending over x(0), x(1),...,x(n)
c*********************************************
#include "precision.h"
      COMP_PRECISION x(0:n),c(0:n,0:m)
      c1=1.0d0
      c4=x(0)-xi
      do k=0,m
         do j=0,n
            c(j,k)=0.0d0
         end do
      end do
      c(0,0)=1.0d0
      do i=1,n
         mn=min(i,m)
         c2=1.0d0
         c5=c4
         c4=x(i)-xi
         do j=0,i-1
            c3=x(i)-x(j)
            c2=c2*c3
            do k=mn,1,-1
               c(i,k)=c1*(k*c(i-1,k-1)-c5*c(i-1,k))/c2
            end do
            c(i,0)=-c1*c5*c(i-1,0)/c2
            do k=mn,1,-1
               c(j,k)=(c4*c(j,k)-k*c(j,k-1))/c3
            end do
            c(j,0)=c4*c(j,0)/c3
         end do
         c1=c2
      end do
      return
      end

