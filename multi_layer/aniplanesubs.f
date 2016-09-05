c
c aniosotropic reflectivity package from Vera Schulte-Pelkum in March 2005
c based on code by Peter Shearer and various others, see individual READMEs 
c and source code comments for details and description
c  
c minor modifications by Thorsten Becker
c
c $Id: aniplanesubs.f,v 1.5 2005/11/28 01:39:23 becker Exp $
c
c
cc ANIPLANESUBS contains subroutines necessary for ANIPLANE
c
c
c ANILAYER calculates the diagonal phase matrix necessary for propagation
c within an anisotropic layer.
c
c   Given:     vslow(6)   = vertical slownesses (1st 3 are negative)
c              h          = layer thickness
c              omega      = frequency
c   Returns:   lamda      = 6x6 diagonal phase matrix
c              lam_d      = 3x3 "downgoing" part of lamda
c              lam_u      = 3x3 "upgoing" part of lamda
c                           Note:  For some anisotropies lamda will have
c                           2 or 4 downgoing parts, not 3!  lam_d and lam_d
c                           are calculated with positive phase terms and
c                           such that any exponential part is always decaying,
c                           and are thus not simply subsets of lamda.
c
c NOTE:  All variables are complex
c
c Requires:  CZEROOUT
c
      subroutine vera_anilayer(vslow,h,omega,lamda,lam_d,lam_u)
      implicit double complex (a-h,o-z)
      double complex vslow(6),lamda(6,6),lam_d(3,3),lam_u(3,3)
      ci=dcmplx(0.0d0,1.0d0)

c      print *,'inside ANILAYER, vslow = ',vslow
c
      call CZEROOUT(lamda,6)
      do i=1,6
         lamda(i,i)=exp(ci*omega*vslow(i)*h)
      enddo
      call CZEROOUT(lam_d,3)
      call CZEROOUT(lam_u,3)
      do  i=1,3
         if (aimag(vslow(i)).le.0.0d0) then
            lam_d(i,i)=exp(-ci*omega*vslow(i)*h)
         else
            lam_d(i,i)=exp(ci*omega*vslow(i)*h)
         end if
         if (aimag(vslow(i+3)).gt.0.0d0) then !*****changed from .ge. 6/17/08
            lam_u(i,i)=exp(ci*omega*vslow(i+3)*h)
         else
            lam_u(i,i)=exp(-ci*omega*vslow(i+3)*h)
         end if
      enddo
c
      return
      end
c
c
c
c
c ANIFACE calculates reflection/transmission coefficients
c for an interface between two anisotropic layers.
c
c   Given:    e1(6,6)      =  eigenvector matrix for layer 1
c             e2(6,6)      =  eigenvector matrix for layer 2
c             p            =  horizontal slowness
c   Returns:  rt(6,6)      =  reflection/tranmission matrix
c             tran_d(3,3)  =  submatrix of rt
c             tran_u(3,3)  =        "
c             ref_d(3,3)   =        "
c             ref_u(3,3)   =        "
c
c NOTE:  All variables are complex!
c
c Requires:  CINVERT,UTILCMAT
c
      subroutine vera_aniface(e1,e2,p,rt,tran_d,tran_u,ref_d,ref_u)
      implicit double complex (a-h,o-z)
      double complex e1(6,6),e2(6,6),einc(6,6),esca(6,6),escainv(6,6)
      double complex rt(6,6),tran_d(3,3),tran_u(3,3),ref_d(3,3),
     &     ref_u(3,3)
      double complex z(6),det(2)
      integer ipvt(6)
      real*8 cond
      
      cone=dcmplx(1.0d0,0.0d0)
      hvel=cone/p
      
c       print*, 'aniface'
c       print*, ((e1(i,j),i=1,3),j=1,3)
c       print*, ' '
c       print*, ((e2(i,j),i=1,3),j=1,3)

c Rearrange e matrices to separate incident and scattered waves
      do  i=1,6
         do  j=1,3
            einc(i,j)=e1(i,j)/hvel !hvel scales to match A. & R.
            einc(i,j+3)=-e2(i,j+3)/hvel
            esca(i,j)=e2(i,j)/hvel
            esca(i,j+3)=-e1(i,j+3)/hvel
            escainv(i,j)=esca(i,j)
            escainv(i,j+3)=esca(i,j+3)
         enddo
      enddo
c
      call CGECO(escainv,6,6,ipvt,cond,z)
c      print *,'e condition #= ',cond
      job=1
      call CGEDI(escainv,6,6,ipvt,det,z,job)
      call CMMULT(escainv,einc,rt,6)

c Set up submatrices
      do  i=1,3
         do  j=1,3
            tran_d(i,j)=rt(i,j)
            ref_d(i,j)=rt(i+3,j)
            ref_u(i,j)=rt(i,j+3)
            tran_u(i,j)=rt(i+3,j+3)
         enddo
      enddo
      

c
      return
      end
