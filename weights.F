c
c     
c $Id: weights.F,v 1.1 2002/10/12 20:37:13 tbecker Exp $
c     
c
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
#include "fortran_precision.h"
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

