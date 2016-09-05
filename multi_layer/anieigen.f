c
c aniosotropic reflectivity package from Vera Schulte-Pelkum in March 2005
c based on code by Peter Shearer and various others, see individual READMEs 
c and source code comments for details and description
c  
c minor modifications by Thorsten Becker
c
c $Id: anieigen.f,v 1.4 2010/03/04 21:51:31 becker Exp $
c
c
cc subroutine ANIEIGEN calculates the eigenvalues and eigenvectors
c for a general anisotropic solid for a given horizontal slowness.
c
c Six solutions are obtained,and are sorted as follows:
c     (1) = downgoing qP    (positive or larger vertical slowness
c     (2) = downgoing qSP         or larger imag. part)
c     (3) = downgoing qSR
c     (4) = upgoing qP      (negative or smaller vertical slowness
c     (5) = upgoing qSP           or smaller imag. part)
c     (6) = upgoing qSR
c
c If no symmetry axis is specified, i.e. sym(i)=0., then the qSP and qSR
c indices correspond to the faster and slower solutions, respectively.
c This sorting scheme has not been entirely debugged, and may give
c unpredictable results, especially for extremely anisotropic models
c and/or evanescent waves.
c
c Program is based on equations which can be found in Keith & Crampin
c (GJRAS 49, 209-223, 1977), Crampin (Wave Motion 3, 343-391,1981), and
c Garmany (GJRAS 75, 565-569, 1983).
c
c Inputs:   rho           =  density
c           c(3,3,3,3)    =  elastic tensor (normalized by density)
c           hslow         =  given horizontal slowness (p1 where p2=0)
c           sym(3)        =  symmetry axis direction,
c                            used to determine qS wave ordering.
c                         =  (0.,0.,0.) for simple ordering by eigenvalue size.
c Returns:  q(3,6)        =  (complex) slowness vectors for the
c                            six solutions.
c           pol(3,6)      =  (complex) polarization vectors,
c                            normalized to unit length.
c           e(6,6)        =  (complex) displacement-stress eigenvector matrix,
c                            assuming unit length polarizations, stress
c                            normalized by i/hslow factor as in Crampin papers.
c           d(6,6)        =  (complex) same as e(6,6), but without the
c                            i/hslow factor normalizing the stresses.
c           dnorm(6,6)    =  (complex) same as d(6,6), but with eigenvectors
c                            normalized by Garmany Eqn. (20).
c           dnorminv(6,6) =  (complex) inverse of ee matrix.
c           energy(3,6)   =  (complex) energy flux vector for 6 waves.
c                            Note that this must be normalized to
c                            get group velocity vector.
c
c Requires:  SGEEV   : calculates eigenvalues, eigenvectors
c            INVERT  : finds inverse of matrix
c            SORTPOL : finds sorting index for solutions
c
      SUBROUTINE vera_anieigen(rho,c,hslow,sym,q,pol,
     &     e,d,dnorm,dnorminv,energy)
      implicit double precision (a-h,o-z)
      real*8 c(3,3,3,3),sym(3)
      double complex q(3,6),q_ran(3,6),pol(3,6),pol_ran(3,6)
      double complex eval(6),eval_ran(6),evec(6,6),evec_ran(6,6)
      double complex e(6,6),d(6,6),dnorm(6,6),dnorminv(6,6),energy(3,6)
      integer index(6)
      real*8 r(3,3),v(3,3),t(3,3),s(3,3),that(3,3),rinv(3,3)
      real*8 scrm1(3,3),scrm2(3,3),scrm3(3,3)
      real*8 big(6,6),work(12)
      double complex p(3,3),ptil(3,3),a(3,3),atil(3,3),atr(3,3),
     &     atiltr(3,3)
      double complex rcom(3,3),thatcom(3,3),vcom(3,3)
      double complex cscr1(3,3),cscr2(3,3),cscr3(3,3),cscr4(3,3),
     &     cscr5(3,3)
      double complex cbigscr1(6,6),cbigscr2(6,6)
      double complex cs1,cs2,cs3,asum,cfact,cxprod,ci,csum,czero

      ci=dcmplx(0.0d0,1.0d0)
      czero=dcmplx(0.0d0,0.0d0)

c      print *,'ANIEIGEN sym = ',sym
c
c Set up all the little matrices in Crampin 1981
      hvel=1.0d0/hslow
      do i=1,3
         do j=1,3
            r(i,j)=c(i,3,j,3)*rho
            rinv(i,j)=r(i,j)
            v(i,j)=c(i,3,j,1)*rho !note error in Crampin 1981 definition
            t(i,j)=c(i,1,j,1)*rho
            that(i,j)=t(i,j)
            if (i.eq.j) that(i,j)=that(i,j)-rho*hvel**2
         enddo
      enddo
      do  i=1,3
         do  j=1,3
            s(i,j)=v(i,j)+v(j,i)
         enddo
      enddo
      call INVERT(rinv,3,determ)
      if (determ.eq.0.0d0) then
         print *,'anieigen: error:  R is singular!'
         stop
      end if
      do  i=1,3                 !scrm1=rinv*s
         do  j=1,3              !scrm2=rinv*that
            scrm1(i,j)=0.0d0
            scrm2(i,j)=0.0d0
            do  k=1,3
               scrm1(i,j)=scrm1(i,j)+rinv(i,k)*s(k,j)
               scrm2(i,j)=scrm2(i,j)+rinv(i,k)*that(k,j)
            enddo
         enddo
      enddo

c
c Set up big matrix for finding eigenvectors
      call ZERO66(big)
      do  i=1,3
         do  j=1,3
            ii=i+3
            jj=j+3
            big(i,j)=-scrm1(i,j)
            big(i,jj)=-scrm2(i,j)
            big(ii,j)=0.0d0
            if (i.eq.j) big(ii,j)=1.0d0
            big(ii,jj)=0.0d0
         enddo
      enddo
c     
c     Get eigenvalues and eigenvectors of big matrix
      call SGEEV(big,6,6,eval_ran,evec_ran,6,work,1,info)
c
c _ran suffix indicates unsorted solutions
c
c Get polarizations and vertical slownesses
c Normalize polarizations to unit length and larger real part
      q1=hslow
      do i=1,6
         q_ran(1,i)=dcmplx(q1,0.0d0)
         q_ran(2,i)=dcmplx(0.0d0,0.0d0)
         q_ran(3,i)=eval_ran(i)/dcmplx(hvel,0.0d0)
      enddo
      do  j=1,6
         cs1=dcmplx(0.0d0,0.0d0)
         sreal=0.0d0
         simag=0.0d0
         do i=1,3
            pol_ran(i,j)=evec_ran(i+3,j)
            cs1=cs1+pol_ran(i,j)*conjg(pol_ran(i,j))
            sreal=sreal+(dreal(pol_ran(i,j)))**2
            simag=simag+(dimag(pol_ran(i,j)))**2
         enddo
         cs1=cdsqrt(cs1)
         do  i=1,3
            pol_ran(i,j)=pol_ran(i,j)/cs1
            evec_ran(i,j)=evec_ran(i,j)/cs1
            evec_ran(i+3,j)=evec_ran(i+3,j)/cs1
            if (simag.gt.sreal) then
               pol_ran(i,j)=pol_ran(i,j)*ci
               evec_ran(i,j)=evec_ran(i,j)*ci
               evec_ran(i+3,j)=evec_ran(i+3,j)*ci
            end if
         enddo
      enddo
c
c Get sorting index and sort everything
      call SORTPOL(q_ran,pol_ran,sym,index)
      do  j=1,6
         eval(j)=eval_ran(index(j))
         do i=1,3
            pol(i,j)=pol_ran(i,index(j))
            q(i,j)=q_ran(i,index(j))
            evec(i,j)=evec_ran(i,index(j))
            evec(i+3,j)=evec_ran(i+3,index(j))
         enddo
      enddo
c
c     Now compute e using Crampin 1981 Eqn. 3.11
c     or Keith & Crampin 1977 Eqn. 11
      do  j=1,3
         do  n=1,6
            e(j,n)=evec(j+3,n)
            e(j+3,n)=dcmplx(0.0d0,0.0d0)
            do m=1,3
               e(j+3,n)=e(j+3,n)+r(j,m)*evec(m,n)+
     &              v(j,m)*evec(m+3,n)
            enddo
         enddo
      enddo
c
c remove i*c factor from stresses in e (Crampin tau factor)
      cs1=ci/dcmplx(hslow,0.0d0)
      do  j=1,6
         do i=1,3
            d(i,j)=e(i,j)
            d(i+3,j)=e(i+3,j)/cs1
         enddo
      enddo
c     
c     Normalize d to dnorm as in Garmany 1983
      do  j=1,6
         csum=czero
         do i=1,3
            csum=csum+d(i,j)*d(i+3,j)
         enddo
         csum=cdsqrt(dcmplx(0.5d0,0.0d0)/csum)
         do  i=1,6
            dnorm(i,j)=d(i,j)*csum
         enddo
      enddo
c
c Inverse of dnorm is easily computed
      do  i=1,3
         do  j=1,3
            dnorminv(i,j)=dnorm(j+3,i)
            dnorminv(i,j+3)=dnorm(j,i)
            dnorminv(i+3,j)=dnorm(j+3,i+3)
            dnorminv(i+3,j+3)=dnorm(j,i+3)
         enddo
      enddo
c
c Now calculate energy flux vector (Keith & Crampin Eqn. 15)
      do  i=1,6
         do j=1,3
            energy(j,i)=dcmplx(0.0d0,0.0d0)
            do  k=1,3
               do  m=1,3
                  do  n=1,3
                     cs1=conjg(pol(k,i))*q(m,i)
     &                    *pol(n,i)
                     cs2=pol(k,i)*conjg(q(m,i))
     &                    *conjg(pol(n,i))
                     cs3=dcmplx(0.25d0*rho*c(j,k,m,n),0.0d0)
                     energy(j,i)=energy(j,i)+cs3*(cs1+cs2)
                  enddo
               enddo
            enddo
         enddo
      enddo
c     
      return
      end
