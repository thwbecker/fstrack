c Subroutine SLOWNESS solves for the slowness, polarization,
c and group velocity vectors for a specified phase direction.
c For the case of general anisotropy, this is solved as an
c eigenvalue problem.
c  Given:
c    c(3,3,3,3) = elastic tensor
c    qhat(3)    = phase (slowness) direction.  This is a unit vector.
c  Returned:
c    qmag(3) = slowness magnitudes (eigenvalues)
c    q(3,3)  = slowness vectors
c    a(3,3)  = polarization vectors (eigenvectors)
c    umag(3) = group velocity magnitudes
c    u(3,3)  = group velocity vectors
c All eigenvalues and eigenvectors are sorted as
c      1) qP
c      2) qSV
c      3) qSH
c Requires: JACOBI
c
      subroutine SLOWNESS(c,qhat,qmag,q,a,umag,u)
      implicit none 
      double precision c(3,3,3,3),qhat(3),qmag(3),q(3,3),a(3,3),
     &     umag(3),u(3,3),err,dotcor,scr,emax,amax
      double precision m(3,3),h(3,3),eigen(3,3)
      integer index(3),i,j,k,l,n
c
      print *,'ok 1'
      do 40 i=1,3
         do 30 k=1,3
            m(i,k)=0.d0
            do 20 j=1,3
               do 10 l=1,3
10                m(i,k)=m(i,k)+c(i,j,k,l)*qhat(j)*qhat(l)
20          continue
            h(i,k)=m(i,k)
30       continue
40    continue
c      call DISP33(h)
c
      err=.00001d0
      call vera_jacobi(h,eigen,3,3,err)
c
      print *,'ok 2'
c Now sort by polarization
c qP is largest eigenvalue (smallest slowness)
      emax=0.d0
      do 50 i=1,3
         if (h(i,i).gt.emax) then
            emax=h(i,i)
            index(1)=i
         end if
50    continue
c qSV has most vertical (001) polarization
      amax=0.d0
      do 60 i=1,3
         if (i.eq.index(1)) go to 60
         if (abs(eigen(3,i)).ge.amax) then
            amax=abs(eigen(3,i))
            index(2)=i
         end if
60    continue
c qSH is the remaining polarization
      do 70 i=1,3
         if (i.ne.index(1).and.i.ne.index(2)) then
            index(3)=i
         end if
70    continue
c
      do 90 i=1,3
         scr=h(index(i),index(i))
         if (scr.lt.0.) then
            print *,'We are going to bomb!!!!!'
            print *,'i,index= ',i,index
            print *,'h= ',h
         end if
         if (scr.ne.0.) then
            qmag(i)=1./sqrt(h(index(i),index(i)))
         else
            qmag(i)=0.
         end if
         do 80 j=1,3
            q(j,i)=qhat(j)*qmag(i)
            a(j,i)=eigen(j,index(i))
80       continue
90    continue
c
      print *,'ok 3'
c Now find group velocity vectors
      do 150 i=1,3
         dotcor=0.d0
         do 130 j=1,3
            u(j,i)=0.d0
            do 120 k=1,3
            do 110 l=1,3
            do 100 n=1,3
               u(j,i)=u(j,i)+c(j,k,l,n)*a(k,i)*q(l,i)*a(n,i)
100         continue
110         continue
120         continue
            dotcor=dotcor+u(j,i)*q(j,i)
130      continue
         umag(i)=0.
         do 140 j=1,3
            if (dotcor.ne.0.) u(j,i)=u(j,i)/dotcor
            umag(i)=umag(i)+u(j,i)**2
140      continue
         umag(i)=sqrt(umag(i))
150   continue
c
      return
      end
c
