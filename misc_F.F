c
c
c invert a 3 x 3 Matrix that 
c
c     input: a
c     output: ainv
c     
c     this routine is now defunct, use the faster C routine instead
c
c
c $Id: misc_F.F,v 1.3 2004/04/23 02:57:45 becker Exp $
c
      subroutine invert3x3(a,ainv)
#include "fortran_precision.h"
      dimension a(3,3),ainv(3,3),am(2,2)
      integer i,j,m,n,p,q
c     det(a)
      deta = a(1,1)*a(2,2)*a(3,3)+a(1,2)*a(2,3)*a(3,1)
     &     +a(1,3)*a(2,1)*a(3,2)-a(1,1)*a(2,3)*a(3,2)
     &     -a(1,2)*a(2,1)*a(3,3)-a(1,3)*a(2,2)*a(3,1)

      if(abs(deta).lt.eps_prec)then
         print *,'invert3x3: runtime error, determinate is:',deta
         print *,'invert3x3: output of (almost) singular input matrix:'
         print *,a
         stop
      endif
      do n=1,3
         do m=1,3
            q=0
            do j=1,2
               q=q+1
               if (q.eq.m) q=q+1
               p=0
               do i=1,2
                  p=p+1
                  if (p.eq.n) p=p+1
                  am(i,j)=a(p,q)
               end do
            end do
            ainv(m,n)=((-1)**(m+n))*
     &           (am(1,1)*am(2,2)-am(1,2)*am(2,1))/deta
         end do
      end do
      return
      end
c
c
c     rotate a fortran 3x3 matrix by euler angles alpha,beta,gamma
c
      subroutine rotate_cart_mat_ftrn(a,b,alpha,beta,gamma)
      COMP_PRECISION a,b,alpha,beta,gamma,rot
      dimension a(3,3),b(3,3),rot(3,3)
c     get the rotation matrix
      call calc_rotmat_cart_ftrn(rot,alpha,beta,gamma)
c     rotate the tensor
      call rotate_mat_ftrn(a,b,rot)
      end
      

c
c     get a fortran style rotation matrix given euler angles a,b,g
c
c
      subroutine calc_rotmat_cart_ftrn(rot,alpha,beta,gamma)
      COMP_PRECISION alpha,beta,gamma,rot,sin_alpha,cos_alpha,
     &     cos_beta,sin_beta,sin_gamma,cos_gamma
      dimension rot(3,3)
      
      sin_alpha=sin(alpha)
      cos_alpha=cos(alpha)
      sin_beta=sin(beta)
      cos_beta=cos(beta)
      sin_gamma=sin(gamma)
      cos_gamma=cos(gamma)

      rot(1,1) = cos_alpha*cos_beta*cos_gamma - sin_alpha*sin_gamma
      rot(1,2) = sin_alpha*cos_beta*cos_gamma + cos_alpha*sin_gamma
      rot(1,3) = -sin_beta*cos_gamma
      rot(2,1) = -cos_alpha*cos_beta*sin_gamma -sin_alpha*cos_gamma
      rot(2,2) = -sin_alpha*cos_beta*sin_gamma +cos_alpha*cos_gamma
      rot(2,3) = sin_beta*sin_gamma
      rot(3,1) = cos_alpha*sin_beta
      rot(3,2) = sin_alpha*sin_beta
      rot(3,3) = cos_beta
      end
c
c     rotate a matrix given a rotation matrix rot
c
      subroutine rotate_mat_ftrn(xin, xout,rot)
      COMP_PRECISION xin,xout,rot,tmp
      dimension tmp(3,3),rot(3,3),xin(3,3),xout(3,3)
      integer i,j,k
c     calculate xin . r^T
      do i=1,3
         do j=1,3
            tmp(i,j)=0.0d0
            do k=1,3
               tmp(i,j) = tmp(i,j) + xin(i,k) * rot(j,k)
            enddo
         enddo
      enddo
c     calculate r . tmp
      do i=1,3
         do j=1,3
            xout(i,j) = 0.0d0
            do k=1,3
               xout(i,j) = xout(i,j) + rot(i,k) * tmp(k,j)
            enddo
         enddo
      enddo

      end
