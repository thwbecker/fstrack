!
!
!  subroutines for D-REX by E. Kaminski as of 01/11/2003
!
!  those are all minor subroutines that have been stripped of 
!  a common block 
!
!  modifications:
!  - Voigt Reuss Hill averaging
!  - Enstatite is fixed
!
!
!
!  $Id: drex_util.f90,v 1.20 2005/07/25 17:33:01 becker Exp becker $
!
!
!
! read in parameters for model simulation
!
! if 
! iread_file == 0: use size3,Xol,tau(1),tau(2),tau(3),tau(4),
!                      tau_ens,Mob,chi
!                  as input
! iread_file == 1: reads from file
! iread_file == 2: use default values
!
!
#include "drex_fconst.h"

subroutine drex_init_para(size,size3,Xol,&
     tau,tau_ens,Mob,chi,lambda,stressexp,&
     alt,del,ijkl,l1,l2,S0,S0_ens,iread_file)

  IMPLICIT NONE
  integer iread_file ! read file or use standard parameters
  INTEGER :: size3, size 
  DOUBLE PRECISION :: lambda, Mob, chi,Xol
  ! lambda = nucleation parameter
  ! Mob = grain mobility
  ! threshold volume fraction for activation of grain boundary sliding
  ! Xol = fraction of olivine in the aggregate

  DOUBLE PRECISION, DIMENSION(3,3,3) :: alt ! \epsilon_{ijk} tensor
  DOUBLE PRECISION, DIMENSION(3,3) :: del ! \delta_{ij} tensor
  INTEGER, DIMENSION(3,3) :: ijkl ! tensor of indices to form Cijkl from Sij
  INTEGER, DIMENSION(6) :: l1,l2 ! tensot of indices to form Sij from Cijkl

  DOUBLE PRECISION, DIMENSION(4) :: tau
  ! RSS for the 4 slip systems of olivine (only 3 are activated)
  DOUBLE PRECISION :: tau_ens ! RSS of enstatite slip system
  DOUBLE PRECISION :: stressexp ! stress exponent for olivine and enstatite

!!! Elastic tensors

  DOUBLE PRECISION, DIMENSION(6,6) :: S0,S0_ens


  if(iread_file .eq. 1)then
     !
     ! read parameters from file
     !
     OPEN(15,file="drex_input.dat")
     ! general parameters
     READ(15,*) size3,Xol ! read in as %
     READ(15,*) tau(1),tau(2),tau(3),tau(4),tau_ens
     READ(15,*) Mob
     READ(15,*) chi
!!! Nucleation parameter set to 5 as in Kaminski and Ribe, 2001
     read(15,*)lambda
     CLOSE(15)
  else if(iread_file .eq. 2)then
     !
     ! assign standard values
     ! 
     size3 = 12 ! total number is size3^3
     Xol = 70d0 !  fraction of olivine (%)
     tau(1)  = 1d0  !! tau1, tau2, tau3, tau4 (olivine), tau_ens (enstatite)
     tau(2)  = 2d0 
     tau(3)  = 3d0 
     tau(4)  = 1d60 
     tau_ens = 1d0
     Mob = 125d0 ! grain boundary mobility
     chi = 0.30d0 !! chi (threshold volume fraction for GBS)
     
!!! Nucleation parameter set to 5 as in Kaminski and Ribe, 2001
     lambda = 5d0

  endif

  ! go to fractional
  Xol = Xol/1d2


!!! initial size
  size = size3**3

!!! stress exponent for olivine and enstatite
  stressexp = 3.5d0

!!! tensor \epsilon_{ijk}

  alt=0d0
  alt(1,2,3) = 1d0 ; alt(2,3,1) = 1d0 ; alt(3,1,2) = 1d0
  alt(1,3,2) = -1d0 ; alt(2,1,3) = -1d0 ; alt(3,2,1) = -1d0

!!! tensor \delta_{ij}

  del=0d0
  del(1,1) = 1d0 ; del(2,2) = 1d0 ; del(3,3) = 1d0

!!! tensors of indices

  ijkl(1,1) = 1 ; ijkl(1,2) = 6 ; ijkl(1,3) = 5
  ijkl(2,1) = 6 ; ijkl(2,2) = 2 ; ijkl(2,3) = 4
  ijkl(3,1) = 5 ; ijkl(3,2) = 4 ; ijkl(3,3) = 3

  l1(1) = 1 ; l1(2) = 2 ; l1(3) = 3
  l1(4) = 2 ; l1(5) = 3 ; l1(6) = 1
  l2(1) = 1 ; l2(2) = 2 ; l2(3) = 3
  l2(4) = 3 ; l2(5) = 1 ; l2(6) = 2
!
!
!!! get reference tensors 
! (those will only get overridden if varying modules are allowed)
!
!
  call drex_get_olivine0_sav_ftrn(s0)
  call drex_get_enstatite0_sav_ftrn(s0_ens)


end subroutine drex_init_para

!
!
!
!!! initialization of orientations - uniformally random distribution
!!! Rmq cos(theta) used to sample the metric Eulerian space
!
! irmode = 0: random, built in random number
! irmode = 1: random, ran1 from numerical recipes
! irmode = 2: oriented, only for testing purposes
!
subroutine drex_init_acs_random_ftrn(size,size3,acs,irmode)

  IMPLICIT NONE

  INTEGER :: i,j1,j2,j3! loop counters
  INTEGER :: size,size3,sizet3,irmode
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ran0,xe1,xe2,xe3
  ! matrix of random numbers used to generate initial random LPO
  DOUBLE PRECISION, DIMENSION(size,3,3) :: acs
  DOUBLE PRECISION :: sin_phi1,sin_phi2,cos_phi1,cos_phi2
  DOUBLE PRECISION :: cos_theta,sin_theta,pi,pif,pih
  parameter (pi = DREX_PI )
  logical :: random_init
  ! matrixes of initial random eulerian angles

  sizet3 = size * 3
  ALLOCATE(ran0(sizet3))
  ALLOCATE(xe1(size),xe2(size),xe3(size))
  pif = pi / DBLE(size3)
  pih = pi / 2.0d0
  if(irmode.eq.0)then
     !
     ! built in random
     !
     call random_number(ran0)
     random_init = .true.
  else if(irmode.eq.1)then
     !
     ! numerical recipes random
     !
     call drex_random_number(ran0,sizet3)
     random_init = .true.
  else if(irmode .eq. 2)then
     random_init = .false. 
  else
     print *,'drex_init_acs_random_ftrn: error: irmode ',irmode,&
          ' undefined'
     stop
  endif
  if(random_init)then
     i = 1
     
     DO j1 =1, size3 
        DO j2 =1, size3 
           DO j3 =1, size3
              ! phi_1 
              xe1(i) = (dble(j1)-ran0(i)) * pif
              ! theta 
              xe2(i) = acos(-1d0 + (dble(j2)-ran0(size+i))/dble(size3)*2d0)
              ! phi_2 angle
              xe3(i) = (dble(j3)-ran0(i+2*size)) * pif
              i = i + 1
           END DO
        END DO
     END DO
  else                          !non-random branch
     print *,'drex_init_acs_random_ftrn: WARNING: using oriented axes initialization'
     j1=1
     DO i =1, size
        xe1(i) = (-1.0d0+2*ran0(j1)) * 0.05d0
        j1=j1+1
        xe2(i) = pif + (-1.0d0+2*ran0(j1)) * 0.05d0
        j1=j1+1
        xe3(i) = (-1.0d0+2*ran0(j1)) * 0.05d0
        j1=j1+1
     enddo
#ifndef DREX_DEBUG
     print *,'drex_init_acs_random_ftrn: exiting, we are not in debug mode'
     stop
#endif
  endif
  
  DO i = 1 , size

     cos_phi1 = cos(xe1(i))
     cos_phi2 = cos(xe3(i))
     sin_phi1 = sin(xe1(i))
     sin_phi2 = sin(xe3(i))
     sin_theta = sin(xe2(i))
     cos_theta = cos(xe2(i))
!!! Direction cosine matrix

     acs(i,1,1)=cos_phi2*cos_phi1-cos_theta*sin_phi1*sin_phi2
     acs(i,1,2)=cos_phi2*sin_phi1+cos_theta*cos_phi1*sin_phi2
     acs(i,1,3)=sin_phi2*sin_theta

     acs(i,2,1)=-sin_phi2*cos_phi1-cos_theta*sin_phi1*cos_phi2
     acs(i,2,2)=-sin_phi2*sin_phi1+cos_theta*cos_phi1*cos_phi2
     acs(i,2,3)=cos_phi2*sin_theta

     acs(i,3,1)=sin_theta*sin_phi1
     acs(i,3,2)=-sin_theta*cos_phi1
     acs(i,3,3)=cos_theta
  END DO

  DEALLOCATE(xe1,xe2,xe3,ran0)

end subroutine drex_init_acs_random_ftrn

!
! init the parameters anew for each pathline
!
subroutine drex_init_pathline_ftrn(rt,odf,acs,acs0,isize)
  implicit none
  integer :: isize
  DOUBLE PRECISION, DIMENSION(isize,3,3) :: acs,acs0
  DOUBLE PRECISION, DIMENSION(isize) :: odf,rt

  rt = 0d0 
  odf = 1d0/dble(isize)         ! random ODF
  acs = acs0                    ! same random initial distribution
end subroutine drex_init_pathline_ftrn
!
! compute the strain rate tensor e and characteristic strain rate from 
! a velocity gradient matrix l
!
subroutine drex_strain_rate_ftrn(l,e,epsnot)

  use drex_nrmod
  implicit none 
  double precision, dimension(3,3) :: l,e,ee
  double precision :: epsnot
  INTEGER :: nrot
  ! number of rotations for the Jacobi
  DOUBLE PRECISION, DIMENSION(3) :: evals
  DOUBLE PRECISION, DIMENSION(3,3) :: evects
  ! eigen values and vectors in jacobi
  e = 0.0d0
  ! get symmetric strain rate tensor from l
  e(1,1) = l(1,1) ; e(2,2) = l(2,2); e(3,3) = l(3,3) 
  e(3,1) = (l(1,3)+l(3,1))/2d0 ; e(1,3) = e(3,1)
  e(1,2) = (l(2,1)+l(1,2))/2d0 ; e(2,1) = e(1,2)
  e(2,3) = (l(3,2)+l(2,3))/2d0 ; e(3,2) = e(2,3)

  !! reference strain rate, store e in ee, since it gets destroyed
  ee = e
  call drex_nr_jacobi(ee,evals,evects,nrot)
  epsnot = MAXVAL(ABS(evals))
  
end subroutine drex_strain_rate_ftrn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Subroutine VHR - Calculates elastic tensor cav_{ijkl} for an  !!!
! aggregate using Voigt (f=0), Reuss (f=1), or VHR (f=0.5) averaging
!
! output is Sav, the average stiffness tensor
!
! this assumes fortran sorting of sav, which doesn't matter since the 
! sav matrix is symmetric
!
!                                        !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE drex_vhr_avg_ftrn(S0,S0_ens,acs,acs_ens,odf,odf_ens,&
     Xol,isize,ijkl,l1,l2,Sav,favg)

  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(6,6) :: S0,Sav,S0_ens,cc6,sc6
  ! stiffness matrix
  INTEGER, DIMENSION(3,3) :: ijkl ! tensor of indices to form Cijkl from Sij
  INTEGER, DIMENSION(6) :: l1,l2 ! tensor of indices to form Sij from Cijkl
  integer :: isize
  ! averaging type,  Voigt (f=0), Reuss (f=1), or VHR (f=0.5)
  double precision :: favg,favg1,weight
  DOUBLE PRECISION, DIMENSION(isize,3,3) :: acs,acs_ens
  DOUBLE PRECISION, DIMENSION(isize) :: odf,odf_ens
  ! volume fraction of the olivine and enstatite grains
  ! fraction of olivine in the aggregate
  DOUBLE PRECISION :: Xol,one_minus_ol
  !
  logical voigt_avg

  !local variable
  INTEGER :: i,j,k,ll,nu,p,q,r,ss

  DOUBLE PRECISION, DIMENSION(3,3,3,3) :: Cavv,Cavr,C0,Cav2,cc4
  
  !  open(99,file='cavg.dat')

  if(abs(favg).lt.1d-8)then
     voigt_avg = .true.        ! voigt avg
  else
     voigt_avg = .false.
  endif
  
  C0 = 0d0 ; Cavv = 0d0 ; Cavr = 0d0;  Sav = 0d0

  one_minus_ol = 1.0d0 - Xol    !en fraction
  favg1 = 1.0d0 - favg          !voigt fraction

  !
  ! Single-xl elastic tensors c0_{ijkl}
  ! go from Voigt to 3,3,3,3 tensor
  !
  ! olivine
  !
  call drex_tens4_ijkl_ftrn(s0,ijkl,c0)
  !
  DO nu = 1 , isize             ! grain loop
     Cav2 = 0d0
     weight=odf(nu)*Xol
     !
     ! get the grain's elasticity tensor from all axes
     !
     DO i = 1 , 3 ; DO j = 1 , 3 ; DO k = 1 , 3 ; DO ll = 1 , 3
        DO p = 1 , 3 ; DO q = 1 , 3 ; DO r = 1 , 3 ; DO ss = 1 , 3
           Cav2(i,j,k,ll) = Cav2(i,j,k,ll) + &
                acs(nu,p,i)*acs(nu,q,j)*acs(nu,r,k)*acs(nu,ss,ll)*C0(p,q,r,ss)
        END DO; END DO ; END DO ; END DO
        Cavv(i,j,k,ll) = Cavv(i,j,k,ll) + Cav2(i,j,k,ll)*weight ! voigt
!        print *,Cav2(i,j,k,ll),weight 
     END DO; END DO ; END DO ; END DO
     if(.not.voigt_avg)then
        ! not only the voigt average needs to be computed
        call drex_c42mat6(Cav2,cc6)  ! convert to 6x6   
        call drex_stiffconv(cc6,sc6)! compute compliance
        call drex_tens4_ftrn(sc6,cc4)! convert back to 3,3,3,3
        DO i = 1 , 3 ; DO j = 1 , 3 ; DO k = 1 , 3 ; DO ll = 1 , 3
           Cavr(i,j,k,ll) = Cavr(i,j,k,ll) + cc4(i,j,k,ll)*weight ! reuss     
        END DO; END DO ; END DO ; END DO
     endif
  END DO
  !
  ! Enstatite-xl elastic tensors c0_{ijkl}
  !
  call drex_tens4_ijkl_ftrn(s0_ens,ijkl,c0)
  !
  ! enstatite
  !
  DO nu = 1 , isize
     Cav2 = 0d0
     weight = odf_ens(nu)*one_minus_ol
     DO i = 1 , 3 ; DO j = 1 , 3 ; DO k = 1 , 3 ; DO ll = 1 , 3
        DO p = 1 , 3 ; DO q = 1 , 3 ; DO r = 1 , 3 ; DO ss = 1 , 3
           Cav2(i,j,k,ll) = Cav2(i,j,k,ll) + &
                acs_ens(nu,p,i)*acs_ens(nu,q,j)*acs_ens(nu,r,k)*acs_ens(nu,ss,ll)*&
                c0(p,q,r,ss)
        END DO; END DO ; END DO ; END DO
        Cavv(i,j,k,ll) = Cavv(i,j,k,ll) + Cav2(i,j,k,ll)*weight !voigt 
!        print *,Cav2(i,j,k,ll),weight
     END DO; END DO ; END DO ; END DO
     if(.not.voigt_avg)then
        call drex_c42mat6(Cav2,cc6)  ! convert to 6x6   
        call drex_stiffconv(cc6,sc6)! compute compliance
        call drex_tens4_ftrn(sc6,cc4)! convert back to 3,3,3,3
        DO i = 1 , 3 ; DO j = 1 , 3 ; DO k = 1 , 3 ; DO ll = 1 , 3
           Cavr(i,j,k,ll) = Cavr(i,j,k,ll) + cc4(i,j,k,ll)*weight ! reuss     
        END DO; END DO ; END DO ; END DO
     endif
  END DO
  if(.not.voigt_avg)then
     !
     ! invert the average compliance, some favg != 0
     !
     ! Average stiffness matrix, convert to 6,6
     ! Reuss avg
     !
     call drex_c42mat6(Cavr,cc6)  ! convert to 6x6   
     call drex_stiffconv(cc6,sc6)       ! invert
     DO i = 1 , 6 ; DO j = 1 , 6
        Sav(i,j) = favg1 * Cavv(l1(i),l2(i),l1(j),l2(j)) + &
             favg * sc6(i,j) 
     END DO; END DO
  else
     ! voigt average
     DO i = 1 , 6 ; DO j = 1 , 6
        Sav(i,j) = Cavv(l1(i),l2(i),l1(j),l2(j)) 
     END DO; END DO
  endif
  

  RETURN

END SUBROUTINE drex_vhr_avg_ftrn
!
! compute the VHR average of two Voigt matrices for volume fraction fvol1 for the first
! materia, fvhr factor 0: Voigt 0.5: VHR 1 : Reuss
!
subroutine drex_st_vhr_avg_ftrn(c1,c2,fvol1,fvhr,cavg)
  implicit none
  double precision, intent(in), dimension(6,6) :: c1,c2 ! input tensors to be averaged
  double precision, intent(in) :: fvol1 ! volume fraction (0...1) of first tensor
  double precision, intent(in) :: fvhr  ! VHR factor 0: Voigt 0.5: VHR 1 : Reuss
  double precision, intent(out), dimension(6,6) :: cavg ! average output tensor
  double precision, dimension(6,6) :: s1,s2,savg ! compliances
  double precision :: fm1
  call drex_stiffconv(c1,s1)! compute compliance
  call drex_stiffconv(c2,s2)! compute compliance
  fm1 = 1.0d0 - fvol1
  cavg = fvol1 * c1 + fm1 * c2 ! average tensor
  savg = fvol1 * s1 + fm1 * s2 ! average compliance
  call drex_stiffconv(savg,s2)! inverse
  fm1 = 1.0d0-fvhr
  cavg = fm1 * cavg + fvhr * s2
end subroutine drex_st_vhr_avg_ftrn

!
! make sure cosine arrays and orientation density functions
! lie within physical bounds
!
subroutine drex_check_phys_limits_ftrn(acs,acs_ens,&
     odf,odf_ens,isize)
  implicit none
  ! input/output
  integer :: isize
  DOUBLE PRECISION, DIMENSION(isize,3,3) :: acs,acs_ens
  DOUBLE PRECISION, DIMENSION(isize) :: odf,odf_ens

  ! local
  integer :: j,j1,j2
  !
  ! limit asc arrays to lie between -1 and 1
  !
  DO j = 1 , isize 
     DO j1 = 1 , 3  
        DO j2 = 1, 3
           IF (acs(j,j1,j2) .GT. 1d0) acs(j,j1,j2) = 1d0
           IF (acs_ens(j,j1,j2) .GT. 1d0) acs_ens(j,j1,j2) = 1d0
           IF (acs(j,j1,j2) .LT. -1d0) acs(j,j1,j2) = -1d0
           IF (acs_ens(j,j1,j2) .LT. -1d0) acs_ens(j,j1,j2) = -1d0
        END DO; 
     END DO; 
  END DO
  !
  ! odf shouldn't be smaller than 0
  !
  DO j = 1 , isize
     IF (odf(j) .LE. 0 ) odf(j) = 0d0
     IF (odf_ens(j) .LE. 0 ) odf_ens(j) = 0d0
  END DO
  !
  ! normalize odf
  !
  odf = odf/SUM(odf)
  odf_ens = odf_ens/SUM(odf_ens)
end subroutine drex_check_phys_limits_ftrn


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Subroutine  eigen - find 3 eigen values of velocity gradient tensor    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE drex_eigen(l,F)
  !
  ! input: l(3,3) velocity gradient tensor
  !

  IMPLICIT NONE

  DOUBLE PRECISION :: a2,a3,Q,R,theta,xx,lambda1,lambda2,lambda3

  DOUBLE PRECISION, DIMENSION(3,3) :: F,Id,l

  Id = 0d0 ; Id(1,1) = 1d0 ; Id(2,2) = 1d0 ; Id(3,3) = 1d0

!!! looking for the eigen-values of L (using tr(l)=0)

  a2 = l(1,1)*l(2,2)+l(2,2)*l(3,3)+l(3,3)*l(1,1) &
       -l(1,2)*l(2,1)-l(2,3)*l(3,2)-l(3,1)*l(1,3)
  a3 = l(1,1)*l(2,3)*l(3,2)+l(1,2)*l(2,1)*l(3,3)+l(1,3)*l(2,2)*l(3,1) &
       -l(1,1)*l(2,2)*l(3,3)-l(1,2)*l(2,3)*l(3,1)-l(1,3)*l(2,1)*l(3,2)

  Q = -a2/3d0
  R = a3/2d0

  IF (ABS(Q) .LT. 1d-9) THEN
     !! simple shear, isa=veloc
     F = -1d0
  ELSE IF (Q**3-R**2 .GE. 0) THEN
     theta = ACOS(R/Q**1.5)
     lambda1 = -2d0*SQRT(Q)*COS(theta/3d0)
     lambda2 = -2d0*SQRT(Q)*COS((theta+2d0*ACOS(-1d0))/3d0)
     lambda3 = -2d0*SQRT(Q)*COS((theta+4d0*ACOS(-1d0))/3d0)

     IF (ABS(lambda1-lambda2) .LT. 1e-13) lambda1=lambda2
     IF (ABS(lambda2-lambda3) .LT. 1e-13) lambda2=lambda3
     IF (ABS(lambda3-lambda1) .LT. 1e-13) lambda3=lambda1

     IF (lambda1 .GT. lambda2 .AND. lambda1 .GT. lambda3) F=matmul(l-lambda2*Id,l-lambda3*Id)
     IF (lambda2 .GT. lambda3 .AND. lambda2 .GT. lambda1) F=matmul(l-lambda3*Id,l-lambda1*Id)
     IF (lambda3 .GT. lambda1 .AND. lambda3 .GT. lambda2) F=matmul(l-lambda1*Id,l-lambda2*Id)

     IF (lambda1 .EQ. lambda2 .AND. lambda3 .GT. lambda1) F=matmul(l-lambda1*Id,l-lambda2*Id)
     IF (lambda2 .EQ. lambda3 .AND. lambda1 .GT. lambda2) F=matmul(l-lambda2*Id,l-lambda3*Id)
     IF (lambda3 .EQ. lambda1 .AND. lambda2 .GT. lambda3) F=matmul(l-lambda3*Id,l-lambda1*Id)

     IF (lambda1 .EQ. lambda2 .AND. lambda3 .LT. lambda1) F=0d0
     IF (lambda2 .EQ. lambda3 .AND. lambda1 .LT. lambda2) F=0d0
     IF (lambda3 .EQ. lambda1 .AND. lambda2 .LT. lambda3) F=0d0

  ELSE

     xx = (SQRT(R**2-Q**3)+ABS(R))**(1d0/3d0)
     lambda1 = -SIGN(1d0,R)*(xx+Q/xx)
     lambda2 = -lambda1/2d0
     lambda3 = -lambda1/2d0
     IF (lambda1 .GT. 1d-9) THEN
        F=matmul(l-lambda2*Id,l-lambda3*Id)
     ELSE
        F = 0d0
     END IF

  END IF

  RETURN

END SUBROUTINE drex_eigen


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Subroutine drex_isacalc_ftrn - 
!!! Calculates ISA orienation at grid point           !!!
!!!
!!! gol only gets modified if F == 0
!!!
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE drex_isacalc_ftrn(isa,gol,l)
  !
  ! input: l: *normalized* velocity gradient tensor, in FORTRAN sorting
  !
  ! output: isa(3): infite strain axis, 
  ! gol: grain orientation lag
  !
  use drex_nrmod
  IMPLICIT NONE

  DOUBLE PRECISION, DIMENSION(3) :: isa
  ! orientation local of infinite strain axis
  double precision :: gol,ltrace

  DOUBLE PRECISION, DIMENSION(3,3) :: F,U,l,lloc
  ! deformation gradient tensor and left strech tensor
   ! velocity gradient tensor
  INTEGER :: nrot
  ! number of rotations for the Jacobi

  DOUBLE PRECISION, DIMENSION(3) :: evals
  DOUBLE PRECISION, DIMENSION(3,3) :: evects
  ! eigen values and vectors in jacobi

!!! Rmq: if GOL is not defined, then isa=0 and GOL=-1d0

!!! calculation of the ISE orientation using Sylvester's formula

!!!
  
  ! make sure l does not have a trace
  ltrace=(l(1,1)+l(2,2)+l(3,3))/3.0d0
  lloc = l
  lloc(1,1)=lloc(1,1)-ltrace; lloc(2,2)=lloc(2,2)-ltrace; lloc(3,3)=lloc(3,3)-ltrace;
!!! Limit deformation gradient tensor for infinite time
  CALL drex_eigen(lloc,F)

  IF (SUM(F) .EQ. -9d0) THEN
!!! isa is flow
     isa = -1d0 
!!! 2. formation of the left-stretch tensor U = FFt
  ELSE IF (SUM(ABS(F)) .EQ. 0) THEN
     isa = 0d0 
     gol = -1d0
  ELSE
     U = MATMUL(F,TRANSPOSE(F)) 

!!! 3. eigen-values and eigen-vectors of U
     !
     ! this destroys U!
     !
     call drex_nr_jacobi(u,evals,evects,nrot)

     IF (evals(1) .EQ. MAXVAL(evals)) THEN
        isa = evects(:,1)
     ELSE IF(evals(2) .EQ. MAXVAL(evals)) THEN
        isa = evects(:,2)
     ELSE
        isa = evects(:,3)
     END IF

  END IF

  RETURN

END SUBROUTINE drex_isacalc_ftrn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! 
! get random numbers using numerical recipes
!
subroutine drex_random_number(xran,n)
  implicit none
  DOUBLE PRECISION :: ran1
  integer :: n
  DOUBLE PRECISION, dimension(n):: xran
  integer :: seed,i
  logical :: init
  save seed,init
  data seed /-1/, init /.false./
  ! use ran1 from numerical recipes
  if(.not.init)then
     ! init ran1
     xran(1) = ran1(seed)
     init = .true.
  endif
  do i=1,n
     xran(i) = ran1(seed)
  enddo
end subroutine drex_random_number


!
! ran1 random number generator from numerical recipes
!
FUNCTION ran1(idum)
  implicit none
  INTEGER :: idum,IA,IM,IQ,IR,NTAB,NDIV
  DOUBLE PRECISION :: ran1,AM,EPS,RNMX
  PARAMETER (IA=16807,IM=2147483647,AM=1.d0/IM,IQ=127773,IR=2836,&
       NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-14,RNMX=1.d0-EPS)
  INTEGER :: j,k,iv(NTAB),iy
  SAVE iv,iy
  DATA iv /NTAB*0/, iy /0/
  if (idum.le.0.or.iy.eq.0) then
     idum=max(-idum,1)
     do j=NTAB+8,1,-1
        k=idum/IQ
        idum=IA*(idum-k*IQ)-IR*k
        if (idum.lt.0) idum=idum+IM
        if (j.le.NTAB) iv(j)=idum
     enddo
     iy=iv(1)
  endif
  k=idum/IQ
  idum=IA*(idum-k*IQ)-IR*k
  if (idum.lt.0) idum=idum+IM
  j=1+iy/NDIV
  iy=iv(j)
  iv(j)=idum
  ran1 = min(AM*iy,RNMX)
  return
END FUNCTION ran1

!
! compute pole figures for three axes and olivine and enstative components
!
! this routine projects the Cartesian system orientations
! of acs into a local spherical system at locp
!
! acs(k,i,j) and acs_ens hold the orientations of the crystallographic axis i
! with respect to the reference frame axis j for grain k
!
! xlocp(r,theta,phi) is the location at which we are to perform the 
!
! if irotate == 0, will leave ODF in Cartesian frame, 
! if irotate == 1, will rotate into N-E-up frame
!
! o1xyz_out: 1: print the direction cosines for [100] olivine and ODF for
!               each grain to a o1xyz.time.dat file
!            2: print the direction cosines for [100], [010], [001] and ODF for olivine
!
! icall: number for o1xyz output file
! time and icall is only passed for information
!
!
SUBROUTINE drex_compute_pdens_ftrn(acs,acs_ens,odf,odf_ens,&
     isize,pdens,pdens_ens,npx,npy,xp,time,irotate,o1xyz_out,&
     icall,iuse_old_pdens)

  IMPLICIT NONE
  integer, intent(in) :: isize,npx,npy,irotate,o1xyz_out,icall,iuse_old_pdens
  double precision :: time

  DOUBLE PRECISION, DIMENSION(isize,3,3) :: acs,acs_ens ! direction cosines
  DOUBLE PRECISION, DIMENSION(isize) :: odf,odf_ens ! orientation density function
  double precision, dimension(3) :: xp ! location vector
  double precision, dimension(3) :: rvec, rvec2, rvec3 ! rotated vectors
  DOUBLE PRECISION, dimension (npx*npy*3) :: pdens,pdens_ens ! output array
                                                     ! with summed projections
  double precision, dimension(3,3) :: rotmat ! for conversion
  double precision :: theta,phi,odf_fac
  integer :: i,j,index,nxny,nxny3,offset
  character*10  timestring,callstring
  character*200 filename
  logical rotate,warned
  save warned
  data warned /.false./

  nxny = npx * npy
  nxny3 = nxny * 3

  if(irotate.ne.0)then
     rotate = .true.
     !
     ! compute a local cartesian system for 
     ! which the summation will be performed
     !
     call drex_calc_rot_cart_basis_ftrn(rotmat,xp)
     !write(*,"(9(f5.2,1x))")rotmat
  else
     if(.not.warned)then
        write(0,*)'drex_pdens: WARNING: not rotating pole figures into geographic system'
        warned=.true.
     endif
     rotate = .false.
  endif

  !
  ! normalize the pdens array by the uniform distribution
  !
  odf_fac =     1.0d0/dble(isize)
  
  if(iuse_old_pdens.ne.0)then
     !
     ! begin summation, using my old routines
     !
     pdens=0.0d0                   !reset density arrays
     pdens_ens=0.0d0
     
     do i=1,isize                  !loop through grains
        offset = 1
        do j=1,3                   !loop through axes
           ! olivine component, this will return phi as azimuth
           call drex_compute_pdens_index(acs(i,j,:),theta,phi,npx,npy,index,&
                rotmat,rotate)
           pdens(offset+index) =   pdens(offset+index) + odf(i)
           ! enstatite component
           call drex_compute_pdens_index(acs_ens(i,j,:),theta,phi,npx,npy,index,&
                rotmat,rotate)
           pdens_ens(offset+index) =   pdens_ens(offset+index) + odf_ens(i)
           offset = offset + nxny 
        enddo
     enddo
     
     pdens = pdens/odf_fac
     pdens_ens = pdens_ens/odf_fac
  else
     
     !
     ! use Jules' routines
     !
     call drex_lambert_proj(isize,acs,odf,rotmat,rotate,npy,pdens)
     call drex_lambert_proj(isize,acs_ens,odf_ens,rotmat,rotate,npy,pdens_ens)
  endif
  if(o1xyz_out.ne.0)then
     write(callstring,'(i2)') icall
     write(timestring,'(f6.2)') time
     if(o1xyz_out.eq.1)then
        filename = 'o1xyz.' // trim(adjustl(callstring)) // '.' // trim(adjustl(timestring)) // '.dat'
        write(0,"('compute_pdens_ftrn: printing ODF for [100] olivine to ',a30)")filename
     else
        filename = 'o3xyz.' // trim(adjustl(callstring)) // '.' // trim(adjustl(timestring)) // '.dat'
        write(0,"('compute_pdens_ftrn: printing ODF for [***] olivine to ',a30)")filename
     endif

     !
     ! write the orientations of the first axis of olivine to file
     !
     OPEN(99,FILE=filename)
     if(rotate)then
        ! header
        write(99,"(i10,1x,i2)")isize,1
        if(o1xyz_out.eq.1)then
           do i=1,isize         !only [100]
              ! rotate into local geographic frame
              call drex_abase_vec2bbase_vec_ftrn(acs(i,1,:),rvec,rotmat)
              write(99,"(3(f9.5  ,1x),4x,f9.5)")rvec,odf(i)/odf_fac
           enddo
        else                    ![100] [010] [001]
           do i=1,isize
              call drex_abase_vec2bbase_vec_ftrn(acs(i,1,:),rvec, rotmat)
              call drex_abase_vec2bbase_vec_ftrn(acs(i,2,:),rvec2,rotmat)
              call drex_abase_vec2bbase_vec_ftrn(acs(i,3,:),rvec3,rotmat)
              write(99,"(9(f9.5  ,1x),4x,f9.5)")rvec,rvec2,rvec3,odf(i)/odf_fac
           enddo
        endif
     else                       ! no rotation, use Cartesian system
        write(99,"(i10,1x,i2)")isize,0
        if(o1xyz_out.eq.1)then
           do i=1,isize
              write(99,"(3(f9.5  ,1x),4x,f9.5)")acs(i,1,:),odf(i)/odf_fac
           enddo
        else
           do i=1,isize
              write(99,"(9(f9.5  ,1x),4x,f9.5)")acs(i,1,:),acs(i,2,:),acs(i,3,:),odf(i)/odf_fac
           enddo
        endif
     endif
     close(99)
  endif
end SUBROUTINE drex_compute_pdens_ftrn

!
! compute the index and theta,phi coordinate for a 
! pole figure given direction cosines a1,a2,a3; 
! those are the x,y,z components. the summation is performed 
! in a local,rotated Cartesian system, based on geographic orientation
! if rotate is true. else, leave in global cartesian
!
! if rotate is false, rotmat will not get referenced
!
! WARNING: this routine assumes that nx and ny don't change 
!
! this will sort longitudes CW from N 0 .. 90 E ... 180 S
!
! uses the upper hemisphere only
!
subroutine drex_compute_pdens_index(cart_vec,theta,phi,&
     nx,ny,index,rotmat,rotate)
  IMPLICIT NONE
  ! input
  double precision, intent(in), dimension(3) :: cart_vec
  double precision, intent(in), dimension(3,3) :: rotmat
  integer, intent(in) :: nx,ny
  logical, intent(in) :: rotate
  ! output
  integer, intent(out):: index
  double precision, intent(out) :: theta,phi
  !
  ! local 
  !
  double precision, dimension(3):: vec
  double precision :: dtheta, dphi
  integer :: o1,o2
  logical init
  save init
  data init /.false./
  save dphi,dtheta
  if(.not.init)then
     dphi =   DREX_TWO_PI/nx     ! phi goes from 0 ... 2pi
     dtheta = DREX_HALF_PI/ny   ! theta from 0 ...pi/2
     init = .true.
  endif
  !
  ! convert the cartesian vector axis, x, into a local, rotated system
  !
  if(rotate)then
     call drex_abase_vec2bbase_vec_ftrn(cart_vec,vec,rotmat)
  else
     vec = cart_vec
  endif
  !
  ! convert unity x,y,z vector to r, theta, phi coordinates 
  !
  theta = acos(vec(3)) ! theta
  !
  ! phi, CW from north
  !
  phi = atan2(vec(2),vec(1)) 
  !
  if(theta.gt.DREX_HALF_PI)then    ! flip theta, if in lower hemisphere
     theta = DREX_PI - theta
     phi = phi + DREX_PI
  endif
  if(phi < 0)then ! adjust phi
     phi = phi + DREX_TWO_PI
  endif
  if(phi > DREX_TWO_PI)then
     phi = phi - DREX_TWO_PI
  endif
  ! lon index from phi, runs from dx/2 ... 360-dx/2
  o1 = int(phi/dphi)        
  ! lat index runs from  90 ... down 
  o2 = int(theta/dtheta)    ! lat index from theta
  if(o1.ge.nx)then
     o1 = nx-1
  endif
  if(o2.ge.ny)then
     o2 = ny-1
  endif  
  index = o2 * nx + o1
end subroutine drex_compute_pdens_index

!
! compute euler angles from orientation cosine matrix(3,3) 
! (see KR, EPSL, 189, 253, 2001, eq.1), euler is phi_1, theta, phi_2
!
subroutine drex_euler(acs,euler)
  implicit none 
  double precision, intent(in), dimension(3,3) :: acs
  double precision, intent(out), dimension(3) :: euler
  double precision :: twopi,pi,aunity
  parameter(pi = DREX_PI, twopi = DREX_TWO_PI,aunity = 1.0d0 - 1d-13)
  !
  ! euler angles: phi_1, theta, phi_2 (or Phi_1, phi, Phi_2)
  !
  ! theta
  euler(2) = acos(acs(3,3))
  if(abs(acs(3,3)).lt.aunity)then
     ! phi_1
     euler(1) = atan2(acs(3,1),-acs(3,2))
     ! phi_2
     euler(3) = atan2(acs(1,3),acs(2,3))
  else
     ! can only determine phi_1 + or - phi_2
     euler(1) = atan2(-acs(2,1),acs(2,2))
     euler(3) = 0.0d0
  endif
  if(euler(1).lt.0.0d0)then
     euler(1)=euler(1)+twopi
  endif
  if(euler(3).lt.0.0d0)then
     euler(3)=euler(3)+twopi
  endif
     
end subroutine drex_euler

subroutine drex_abase_vec2bbase_vec_ftrn(a,b,basis)
  IMPLICIT NONE
  double precision, dimension(3) :: a,b
  double precision, dimension(3,3) :: basis

  b(1) = basis(1,1) * a(1) + basis(1,2) * a(2) + basis(1,3) * a(3);
  b(2) = basis(2,1) * a(1) + basis(2,2) * a(2) + basis(2,3) * a(3);
  b(3) = basis(3,1) * a(1) + basis(3,2) * a(2) + basis(3,3) * a(3);

end subroutine drex_abase_vec2bbase_vec_ftrn

!
!
subroutine drex_calc_rot_cart_basis_ftrn(rot,rvec)
  IMPLICIT NONE
  double precision, intent(in),dimension(3) :: rvec ! r,theta,phi vector
  double precision, intent(out),dimension(3,3) :: rot
  double precision :: alpha,beta,gamma
  integer :: coord_convention
  parameter(coord_convention = DREX_REG_CART_CONVENTION)
  !
  ! get the euler angles for this particular location
  !
  call drex_calc_euler_rad_from_xp(rvec,alpha,beta,gamma, &
       coord_convention)
  call drex_calc_rotmat_cart(rot,alpha,beta,gamma) !get the rotation matrix
end subroutine drex_calc_rot_cart_basis_ftrn

!
!
! compute some statistical properties of the a axes grain orientations for 
! weighted and unweighted means
!
!
SUBROUTINE drex_compute_avg_axes_ftrn(acs,acs_ens,odf,odf_ens,&
     isize,xp,irot,avg_axes,avg_axes_ens,dev_axes,dev_axes_ens,&
     jindex)
  
  IMPLICIT NONE
  integer, intent(in) :: isize,irot
  DOUBLE PRECISION, intent(in),DIMENSION(isize,3,3) :: acs,acs_ens ! direction cosines
  DOUBLE PRECISION, intent(in),DIMENSION(isize) :: odf,odf_ens ! orientation density function
  ! unweighted and weighted for a
  double precision, intent(out), dimension(2,3) :: avg_axes,avg_axes_ens ! mean axes, unweighted and weighted
  double precision, intent(out), dimension(2) :: dev_axes,dev_axes_ens ! average deviation from mean, based
                                                                       ! on weighted and unweighted properties
  double precision, intent(out) :: jindex ! J index
  double precision, dimension(3) :: xp ! location vector
  ! local
  ! for j index
!  integer :: jsindex(3),iaxis,idim
!  double precisions :: sum1,sum2
  integer :: np,npm1,nphm1,ipe,ite
  double precision :: pi,twopi,dphi,theta,phi
  parameter (np = 72, pi = DREX_PI, twopi = DREX_TWO_PI)
  real :: jsum(0:np-1,0:np/2-1)
  double precision, dimension(4,3) :: rvec
  double precision, dimension(3,3) :: rotmat
  double precision, dimension(2) :: odfsum,diff
!  double precision, dimension(3) :: euler
  double precision :: theta_e(0:np/2-1),phi_e(0:np-1),&
       sxd,xdist,sin_te(0:np/2-1)
       
  logical rotate,init
  integer igrain,itype
  save init,sxd
  data init /.false./
  !
  ! for j index
  !
  dphi = twopi/dble(np)
  nphm1= np/2-1
  npm1 = np - 1
  !
  ! change to logical type for rotation
  !
  if(irot.eq.1)then 
     rotate = .true.
     !
     ! compute a local cartesian system for 
     ! which the summation will be performed
     !
     call drex_calc_rot_cart_basis_ftrn(rotmat,xp)
  else 
     rotate = .false. 
  endif
  !
  ! total orientation density
  !
  odfsum(1) = sum(odf)
  odfsum(2) = sum(odf_ens)
  ! init
  avg_axes = 0.0d0; avg_axes_ens = 0.0d0
  dev_axes = 0.0d0; dev_axes_ens = 0.0d0
  !
  ! we do the summing/mean computations in the original coordinate frame, and rotate later
  !
  do igrain=1,isize             !loop through grains
     !
     ! use the ODF as weighting for the mean, only [100] axes
     !
     ! sum
     avg_axes(1,:) =     avg_axes(1,:) +     acs(igrain,1,:)
     avg_axes_ens(1,:) = avg_axes_ens(1,:) + acs_ens(igrain,1,:) 
     avg_axes(2,:) =     avg_axes(2,:) +     acs(igrain,1,:)     * odf(igrain)
     avg_axes_ens(2,:) = avg_axes_ens(2,:) + acs_ens(igrain,1,:) * odf_ens(igrain)
 
  enddo
  !
  ! normalize averages, hence don't need to divide by odfsum() weights
  !
  do itype = 1,2
     call drex_normalize_vector(avg_axes(itype,:),3) !those are the mean a axis orientations
     call drex_normalize_vector(avg_axes_ens(itype,:),3)
  enddo
  !
  ! olivine [100] J index
  !
  ! 5 deg width of Gaussian for J
  !
  if(.not.init)then
     ! initialize some helper arrays for the Euler angles,
     ! which we parameterize by phi and theta coordinates
     phi = dphi/2.0d0
     do ipe=0,npm1                 
        if(ipe.le.nphm1)then
           theta_e(ipe) = phi
           sin_te(ipe) = sin(phi)
        endif
        phi_e(ipe) = phi
        phi = phi + dphi
     enddo
     ! width of Gaussian
     sxd=1d0/(5D0/180d0*DREX_PI/sqrt(log(2D0)))**2
     init = .true.
  endif
  

  jsum = 0.0d0
  do ipe=0,npm1                 ! euler phi1
     do ite=0,nphm1              ! euler_theta loop
        do igrain=1,isize             !loop through grains
           !
           !
           ! J index summation
           !
           ! colat and longitude of grain
           theta = acos(acs(igrain,1,3))
           phi = atan2(acs(igrain,1,2),acs(igrain,1,1))
           !
           ! compute distance squared between grain and euler angle orientation
           !
           xdist = sin((theta_e(ite)-theta)/2.0d0)**2
           xdist = xdist + sin((phi-phi_e(ipe))/2.0d0)**2 * sin(theta) * sin_te(ite)
           xdist = 4.0d0 * asin(sqrt(xdist))**2
           ! add contribution
           jsum(ipe,ite) = jsum(ipe,ite) + odf(igrain) * exp(-xdist*sxd)
        enddo
     enddo
  enddo

  ! if phi1 (0;2pi), theta (0;pi), phi2(0;2pi)
  ! gets integrates like sin(theta)dtheta dphi1 dphi2, the
  ! results if 8 pi^2
  ! the J index is normalized by pi^2 (and not 8pi^2)???
  
  jindex = 0d0
  do ipe=0,npm1                  !phi 
     do ite=0,nphm1              !theta
        jindex = jindex + sin_te(ite) * jsum(ipe,ite)**2 
     enddo                     !end theta
  enddo                        
  jindex = jindex * dphi**2/pi**2
  !
  ! compute avg dot products to get deviations from mean
  !
  do igrain=1,isize
     diff(1) = abs(dot_product(avg_axes(1,:),acs(igrain,1,:) ))
     diff(2) = abs(dot_product(avg_axes_ens(1,:),acs_ens(igrain,1,:) ))
     ! unweighted
     dev_axes(1)     = dev_axes(1)     + diff(1) 
     dev_axes_ens(1) = dev_axes_ens(1) + diff(2) 
     ! weighted
     dev_axes(2)     = dev_axes(2)     + diff(1) * odf(igrain)
     dev_axes_ens(2) = dev_axes_ens(2) + diff(2) * odf_ens(igrain)
  enddo
  !
  ! define the deviation as (1-<a * <a>>)
  !
  do itype=1,2
     dev_axes(itype) =     (1.0d0-dev_axes(itype)/isize)
     dev_axes_ens(itype) = (1.0d0-dev_axes_ens(itype)/isize)
  enddo
  !
  if(rotate)then             ! rotate mean axes into local system, ie. E-N-U
     do itype = 1,2
        call drex_abase_vec2bbase_vec_ftrn(avg_axes(itype,:),    rvec(1,:),rotmat)
        call drex_abase_vec2bbase_vec_ftrn(avg_axes_ens(itype,:),rvec(2,:),rotmat)
        avg_axes(itype,:)     = rvec(1,:)
        avg_axes_ens(itype,:) = rvec(2,:)
     enddo
  endif
 
end SUBROUTINE drex_compute_avg_axes_ftrn

subroutine drex_normalize_vector(x,n)
  implicit none
  integer n
  double precision, dimension(n) :: x

  x=x/sqrt(sum(x(:)**2))
end subroutine drex_normalize_vector


subroutine drex_stiffconv(c, s)
  !
  !  Subroutine to convert elasticity tensor c(i,j) to 
  !  compliance s(i,j)
  !  or vice-versa, tensor to be converted first
  
  !  From Dan Raymer 3/00; minor mods by dkb to
  !  clean up double precision lingo
  !  
  !  From Donna 06/05
  !
  !  minute changes by TWB
  !
  use drex_nrmod
  implicit none 
  integer :: np,i,j
  double precision, intent(in) :: c(6,6)
  double precision, intent(out) :: s(6,6)
  parameter (np=6)
  double precision, dimension(np,np) :: rmat
  integer, dimension(np) :: indx
  double precision :: d
  !
  rmat = c                      ! make copy
  s = 0.0d0                     ! initialize
  do i=1,np 
     s(i,i)=1.d0 
  enddo
  !
  !  Actual inversion routine
  !
  call ludcmp(rmat,np,indx,d)
  do j=1,np            
     call lubksb(rmat,np,indx,s(1,j))
  enddo
  
  return
END subroutine drex_stiffconv

   

!
! print Cijkl tensor to filep
!
subroutine drex_print_cijkl_ftrn(cij,filep)
  implicit none
  double precision :: cij(3,3,3,3)
  integer :: i,j,k,l,filep
! prints the indices
!  write(filep,'(9(4(i1),1x))') ((((i,j,k,l,&
!       i=1,3),j=1,3),k=1,3),l=1,3)
  write(filep,'(1p,9(e14.6,1x))') ((((cij(i,j,k,l),&
       i=1,3),j=1,3),k=1,3),l=1,3)
  
end subroutine drex_print_cijkl_ftrn
!
!  function reads the upper right half of a Cij elastic
!  voigt matrix (6x6) 
!
! and computes Love anisotropy parameters for an equivalnet VTI medium that has been
! averaged over all azimuths
!
! convention follows Montagner and Nataf (JGR, 1986)
! 
! return the paramaters in table 2, for cos and sin respectively
! also returns l=rho v_sv^2 and n=rho v_sh^2
!
! this routine does depend on the reference frame!
!
!
subroutine drex_compute_swpar_ftrn(cij,gl,cn,ca,ba,hf,l,n,a,c,f)
  implicit none
  double precision, dimension(6,6) :: cij ! input
  double precision, dimension(2) :: gl,cn,ca,ba,hf !output
  double precision,intent(out) :: l,n,a,c,f
  ! local
  double precision bc,bs,gc,gs,hc,hs,cc,cs
  !
  ! constant terms
  a = 0.375d0*(cij(1,1)+cij(2,2)) + 0.25d0 * cij(1,2) + .5d0 * cij(6,6) ! rho v_PH^2

  c = cij(3,3)                                                          ! rho v_PV^2

  f = .5d0 * (cij(1,3)+cij(2,3)) ! for ellipticity
  !
  l = .5d0 * (cij(4,4)+cij(5,5)) ! rho v_SV^2

  n = 0.125d0 * (cij(1,1)+cij(2,2)) - 0.25d0* cij(1,2) + .5d0 * cij(6,6) ! rho v_SH^2
  !
  ! 2 theta 
  bc = .5d0 * (cij(1,1) - cij(2,2)) ! cos 2t
  bs = cij(1,6) + cij(2,6)          ! sin 2t
  gc = .5d0 * (cij(5,5) - cij(4,4)) ! 
  gs = cij(5,4)
  hc = .5d0 * (cij(1,3) - cij(2,3))
  hs = cij(3,6)
  !
  ! 4 theta 
  cc = 0.125d0*(cij(1,1) + cij(2,2)) - 0.25d0*cij(1,2) - .5d0*cij(6,6)
  cs = 0.5d0*(cij(1,6) - cij(2,6))
  !
  ! also compute the relative strength parameters from Table two of M&N(1986)
  !
  gl(1) = gc/l;gl(2) = gs/l     ! Rayleigh wave third 2phi term

  cn(1) = cc/n;cn(2) = cs/n     ! Love wave 4phi term 

  ca(1) = cc/a;ca(2) = cs/a     ! Rayleigh 4phi

  ba(1) = bc/a;ba(2) = bs/a     ! Rayleigh first 2phi term
  hf(1) = hc/f;hf(2) = hs/f     ! Rayleigh second 2phi term

end subroutine drex_compute_swpar_ftrn

