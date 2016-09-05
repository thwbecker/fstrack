c $Id: gramschmidt.F,v 1.3 2002/03/18 20:16:44 becker Exp $
      subroutine gramschmidt(work,alpha)
#include "fortran_precision.h"
c___________________________________
c     given the work array containing the components of 
c     position (r,theta,phi), then the deformation matrix
c     def(i,j) in C order (i.e. j varies most quickly)
c     orthonormalize def(i,j) with Gram-Schmidt
c     e.g. p. 148 of Ott, _Chaos in Dynamical Systems_
c     and output the length (alpha(i)) of each of the 
c     column vectors
c__________________________________
      dimension work(*),alpha(3)

c     We'll take our originally orthogonal vectors to be the unit
c     vectors in the r, theta, and phi directions.  That means that
c     the columns of def(i,j) are the originally orthogonal vectors
c     that we want to orthonormalize

c     for each column vector
      do ivec=1,3
         !if this isn't the first column vector...
         if (ivec.gt.1) then
c     for each previously normalized vector...
            do ipnv=1,ivec-1
c     take dot product of vector with pnv
               dotprod=DEF1(ivec)*DEF1(ipnv)+
     &              DEF2(ivec)*DEF2(ipnv)+
     &              DEF3(ivec)*DEF3(ipnv)
c     and subtract off projection along basis vector
               DEF1(ivec)=DEF1(ivec)-
     &              dotprod*DEF1(ipnv)
               DEF2(ivec)=DEF2(ivec)-
     &              dotprod*DEF2(ipnv)
               DEF3(ivec)=DEF3(ivec)-
     &              dotprod*DEF3(ipnv)
            end do
         end if
c     calculate length of orthogonalized vector
         alpha(ivec)=sqrt(
     &        DEF1(ivec)**2+
     &        DEF2(ivec)**2+
     &        DEF3(ivec)**2)
c     and normalize it
         DEF1(ivec)=DEF1(ivec)/alpha(ivec)
         DEF2(ivec)=DEF2(ivec)/alpha(ivec)
         DEF3(ivec)=DEF3(ivec)/alpha(ivec)
      end do

      return
      end
