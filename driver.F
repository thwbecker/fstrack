      program driver
#include "fortran_precision.h"

c
c     fortran program to test RK part
c
c     $Id: driver.F,v 1.7 2002/03/18 20:16:44 becker Exp $
c
      parameter (nr=10,ntheta=10,nphi=10,nv=nr*ntheta*nphi)

      VPRECF vr,vtheta,vphi

      dimension vr(nv),vtheta(nv),vphi(nv),rlevels(nr),
     &     xpart(3), def(9),strain(9)

      dimension tprime(200),wwork(12,200),work(12)

      common /path/ kmax,kount,dxsav,tprime,wwork

      kmax=0

      dtheta=pi/(ntheta-1)
      dphi=2.d0*pi/nphi
      dr=1.d0/(nr-1)
      do i=1,nr
         rlevels(i)=(i-1)*dr
      end do
      index=0
      psi=0.d0
      romega1=twopi*sin(psi)
      romega3=twopi*cos(psi)
      do ir=1,nr
         r=rlevels(ir)
         twopir=twopi*r
         do itheta=1,ntheta
            theta=(itheta-1)*dtheta
            temp=twopi*sin(theta)
            sintheta=sin(theta)
            do iphi=1,nphi
               sinphi=sin(phi)
               cosphi=cos(phi)
               phi=(iphi-1)*dphi
               index=index+1
               vr(index)=0.d0
               vtheta(index)=-r*romega1*sinphi
               vphi(index)=-r*
     &              (-romega3*sintheta+romega1*costheta*cosphi)
            end do
         end do
      end do

      do i=1,nr
         rlevels(i)=(i-1)*dr
      end do
      
      XTRACER(1)=.5d0
      XTRACER(2)=pi/2.d0
      XTRACER(3)=0.d0

      index=0
      do i=1,3
         DEF1(i)=0.d0
         DEF2(i)=0.d0
         DEF3(i)=0.d0
      end do
      DEF1(1)=1.d0
      DEF2(2)=1.d0
      DEF3(3)=1.d0

      nt=10
      dt=1.d0/nt
      write (6,*) 'XTRACER'
      write (6,*) 0.d0,(XTRACER(i),i=1,3)
      !write (6,*) (def(j),j=1,9)
      index=0
      do i=1,3
         do j=1,3
            index=index+1
            strain(index)=0.d0
            strain(index)=DEF1(i)*DEF1(j)+DEF2(i)*DEF2(j)
     $           +DEF3(i)*DEF3(j)
            if (i.eq.j) strain(index)=strain(index)-1.d0
            strain(index)=strain(index)/2.d0
         end do
      end do
c      write (6,*) (strain(i),i=1,9)
      eps=1.e-5
      ifellipse=1
      do it=1,nt
c     repack
         call ellsphere(vr,vtheta,vphi,nr,ntheta,nphi,dtheta,dphi,
     &        rlevels,work,dt,ifellipse,eps)
         write (6,*) 'time=',it*dt
         write (6,*) 'XTRACER'
         write (6,*) (XTRACER(j),j=1,3)
         !write (6,*) (def(j),j=1,9)
         !calculate lagrangian strain tensor e=.5*(F^T dot F - I) 
         index=0
         do i=1,3
            do j=1,3
               index=index+1
               strain(index)=0.d0
               strain(index)=DEF1(i)*DEF1(j)+DEF2(i)*DEF2(j)
     $              +DEF3(i)*DEF3(j)
               if (i.eq.j) strain(index)=strain(index)-1.d0
               strain(index)=strain(index)/2.d0
            end do
         end do
         write (6,*) 'F'
         write (6,*) (DEF1(j),j=1,3)
         write (6,*) (DEF2(j),j=1,3)
         write (6,*) (DEF3(j),j=1,3)
         write (6,*) 'strain'
         do i=0,2
            write (6,*) (strain(i*3+j),j=1,3)
         end do
         detf=DEF1(1)*DEF2(2)*DEF3(3)+
     &        DEF1(2)*DEF2(3)*DEF3(1)+
     &        DEF1(3)*DEF2(1)*DEF3(2)-
     &        DEF1(1)*DEF2(3)*DEF3(2)-
     &        DEF1(2)*DEF2(1)*DEF3(3)-
     &        DEF1(3)*DEF2(2)*DEF3(3)
         write (6,*) 'determinant of def=',detf
      end do
      end

      
