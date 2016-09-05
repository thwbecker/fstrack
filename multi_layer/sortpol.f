c
c aniosotropic reflectivity package from Vera Schulte-Pelkum in March 2005
c based on code by Peter Shearer and various others, see individual READMEs 
c and source code comments for details and description
c  
c minor modifications by Thorsten Becker
c
c $Id: sortpol.f,v 1.2 2005/11/28 01:39:23 becker Exp $
c
c
cc SORTPOL helps sort polarization vectors by returning an index
c
c Inputs:    q(3,6)    =  (complex) slowness vectors for 6 solutions
c            cpol(3,6) =  (complex) polarization vectors for 6 solutions
c            sym(3)    =  (real) hexagonal symmetry axis orientation
c                      =  (0,0,0) for simple sorting by size of q
c Returns    index(6)  =  index for sorting 6 solutions
c
c For hexagonal symmetry, solutions are sorted as:
c     index(1) = downgoing qP position   (positive vertical slowness
c     index(2) = downgoing qSP position         or larger imag. part)
c     index(3) = downgoing qSR position
c     index(4) = upgoing qP position     (negative vertical slowness
c     index(5) = upgoing qSP position          or smaller imag. part)
c     index(6) = upgoing qSR position
c
c If no symmetry axis is specified, i.e. sym(i)=0., then the qSP and qSR
c indices correspond to the faster and slower solutions, respectively.
c This sorting scheme has not been entirely debugged, and may give
c unpredictable results, especially for extremely anisotropic models
c and/or evanescent waves.
c
      SUBROUTINE SORTPOL(q,cpol,sym,index)
      implicit double precision (a-h,o-z)
      double complex cpol(3,6),q(3,6)
      real*8 pol(6,3),sym(3),f_real(6),f_imag(6)
      real*8 ab_real(6),ab_imag(6),dot(6)
      integer index(6),list(6)
      logical simple
c      print *,'SORTPOL sym = ',sym
      if (sym(1).eq.0.0d0.and.sym(2).eq.0.0d0.and.sym(3).eq.0.0d0) then
         simple=.true.
      else
         simple=.false.
      end if
c
c Set up unit normalized real polarization vector (pol), etc.
c      print *,'SORTPOL sym=',sym
      do 20 i=1,6
            sum=0.0d0
            do 10 j=1,3
               pol(i,j)=dreal(cpol(j,i))
10             sum=sum+pol(i,j)**2
            if (sum.ne.0.0d0) then
               sum=dsqrt(sum)
               do 15 j=1,3
15                pol(i,j)=pol(i,j)/sum
            end if
            list(i)=i
            f_real(i)=dreal(q(3,i))
            f_imag(i)=dimag(q(3,i))
            ab_real(i)=dabs(f_real(i))
            ab_imag(i)=dabs(f_imag(i))
c            print *,' '
c            print *,'q=',real(q(1,i)),real(q(2,i)),real(q(3,i))
c            print *,'pol=',pol(i,1),pol(i,2),pol(i,3)
c            scr=real(q(1,i))*pol(i,1)+real(q(2,i))*pol(i,2)+
c     &          real(q(3,i))*pol(i,3)
c            print *,'q dot pol=',scr
            if (.not.simple) then
              qsr1=dreal(q(2,i))*sym(3)-dreal(q(3,i))*sym(2) !qsr polarization is
              qsr2=dreal(q(3,i))*sym(1)-dreal(q(1,i))*sym(3) !normal to both sym
              qsr3=dreal(q(1,i))*sym(2)-dreal(q(2,i))*sym(1) !and slowness vectors
              scr=dsqrt(qsr1**2+qsr2**2+qsr3**2)
              if (scr.ne.0.0d0) then
                 qsr1=qsr1/scr              !normalize qsr to unit length
                 qsr2=qsr2/scr
                 qsr3=qsr3/scr
              end if
c              print *,'qsr pol=',qsr1,qsr2,qsr3
              dot(i)=abs(pol(i,1)*qsr1+pol(i,2)*qsr2+pol(i,3)*qsr3)
c              print *,'dot=',dot(i)
            end if
20    continue
c
c sort into list by increasing ab_real
      do 50 i=1,6
      do 40 j=1,5
         if (ab_real(list(j)).gt.ab_real(list(j+1))) then
            ii=list(j)
            list(j)=list(j+1)
            list(j+1)=ii
         end if
40    continue
50    continue
c
c sort solutions with near-zero ab_real by decreasing ab_imag
      do 70 i=1,6
      do 60 j=1,5
         if (ab_real(list(j)).lt.0.001d0.and.
     &        ab_real(list(j+1)).lt.0.001d0
     &       .and.ab_imag(list(j)).lt.ab_imag(list(j+1))) then
            ii=list(j)
            list(j)=list(j+1)
            list(j+1)=ii
         end if
60    continue
70    continue
c
c qP is now in list(1) and list(2), qS is in list(3) to list(6)
c Now sort qS into qSP and qSR
      if (simple) then
         do 90 i=1,4  !sort qS with non-zero ab_real by decreasing f_real
         do 80 j=3,5
            if (ab_real(list(j)).ge..001d0.and.
     &           ab_real(list(j+1)).ge..001d0
     &          .and.f_real(list(j)).lt.f_real(list(j+1))) then
               ii=list(j)
               list(j)=list(j+1)
               list(j+1)=ii
            end if
80       continue
90       continue
      else
         do 110 i=1,4   !sort qS solutions by increasing dot(i)
         do 100 j=3,5
            if (dot(list(j)).gt.dot(list(j+1))) then
               ii=list(j)
               list(j)=list(j+1)
               list(j+1)=ii
            end if
100      continue
110      continue
      end if
c
c qP  is now at list(1) and list(2)
c qSP is now at list(3) and list(4)
c qsR is now at list(5) and list(6)
c order each by decreasing f_real or by decreasing f_imag
      do 120 i=1,5,2
         if (f_real(list(i)).lt.f_real(list(i+1))) then
            ii=list(i)
            list(i)=list(i+1)
            list(i+1)=ii
         end if
         if (ab_real(list(i)).lt.0.001d0.and.
     &        ab_real(list(i+1)).lt.0.001d0
     &       .and.f_imag(list(i)).lt.f_imag(list(i+1))) then
            ii=list(i)
            list(i)=list(i+1)
            list(i+1)=ii
         end if
120   continue
c
      index(1)=list(1)
      index(4)=list(2)
      index(2)=list(3)
      index(5)=list(4)
      index(3)=list(5)
      index(6)=list(6)
c
      if (.not.simple) then
         if (f_real(list(2)).gt.0.0d0.or.f_real(list(4)).gt.0.0d0.or.
     &       f_real(list(6)).gt.0.0d0) then
            print *,'***Possible problem in SORTPOL'
            print *,'list=',list
            print *,'f_real=',f_real
            print *,'dot=',dot
            print *,'sym=',sym
            print *,'pol=',pol
            print *,'q=',q
         end if
      end if
c
      return
      end
