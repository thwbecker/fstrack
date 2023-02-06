!
!
!  subroutines for D-REX by E. Kaminski as of 01/11/2003
!
!  minor changes by TWB
!
!  $Id: drex_deriv.f90,v 1.10 2011/04/12 06:18:38 becker Exp $
!
!  input:
!
!  acs, acs_ens: directions cosines for olivine and enstatite
!  odf, odf_ens: orientation density functions
!
!  output:
!
!  dotacs, dotacs_ens: derivative of direction cosine
!  dotodf, dotodf_ens: derivative of orientation density functions
!
!
! WARNING: this routine has been modified such that the derivatives 
! are already scale by epsnot. before one had to add quantities like
!
! u(t+dt) = u + du * dt * epsnot
!
! where du is the output of this routine. now you can simply perform the 
! operation:
!
! u(t+dt) = u + du * dt
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine DERIV, calculation of the rotation vector and slip rate     !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE drex_deriv_ftrn(l,e,epsnot,size,acs,acs_ens, &
     dotacs,dotacs_ens,odf,odf_ens,dotodf,dotodf_ens,&
     rt,rt_ens,tau,tau_ens,stressexp,alt,Xol,lambda,Mob,chi,&
     isave_gamma,gamma_save)

  IMPLICIT NONE

!!! LPO calculation

  DOUBLE PRECISION, DIMENSION(3,3) :: l,e
  ! velocity gradient tensor and strain rate tensor

  DOUBLE PRECISION :: epsnot
  ! reference strain rate

!!! Dynamic recrystallization

  INTEGER :: size ! size = size3^3
  ! number of points in the (metric) Eulerian space

  DOUBLE PRECISION :: lambda, Mob, chi
  ! lambda = nucleation parameter
  ! Mob = grain mobility
  ! threshold volume fraction for activation of grain boundary sliding

  DOUBLE PRECISION :: Xol
  ! fraction of olivine in the aggregate

  DOUBLE PRECISION, DIMENSION(3,3,3) :: alt ! \epsilon_{ijk} tensor

  DOUBLE PRECISION, DIMENSION(4) :: tau
  ! RSS for the 4 slip systems of olivine (only 3 are activated)
  DOUBLE PRECISION :: tau_ens ! RSS of enstatite slip system
  DOUBLE PRECISION :: stressexp ! stress exponent for olivine and enstatite

  DOUBLE PRECISION, DIMENSION(size):: odf,dotodf

  REAL, dimension(size*3) :: gamma_save
  integer :: isave_gamma,ioff

  ! volume fraction of the olivine grains
  DOUBLE PRECISION, DIMENSION(size):: odf_ens,dotodf_ens
  ! volume fraction of the enstatite grains

  DOUBLE PRECISION, DIMENSION(size) :: rt,rt_ens
  ! dislocations density for olivine and enstatite

  DOUBLE PRECISION, DIMENSION(size,3,3) :: dotacs,acs
  DOUBLE PRECISION, DIMENSION(size,3,3) :: dotacs_ens,acs_ens
  !! matrix of direction cosine

  !
  ! local variables start here
  !

  INTEGER :: i,i1,i2,i3,i4,j,k ! counters
  INTEGER :: imax,iint,imin,iinac ! dummies
  INTEGER , DIMENSION(1) :: ti ! reordering array

  DOUBLE PRECISION :: Emean,rt1,rt2,rt3,Emean_ens,rt0_ens
  !! surface averaged aggregate NRJ
  !! dislocation density for each slip system

  DOUBLE PRECISION :: gam0
  ! slip rate on the softest slip system

  DOUBLE PRECISION :: R1,R2
  DOUBLE PRECISION :: qint,qmin,rat
!!! dummies

  DOUBLE PRECISION, DIMENSION(4) :: bigi,q,qab ! intermediates for G calc

  DOUBLE PRECISION, DIMENSION(4) :: gam
  ! ratios of strain between softest slip system and slip system s for Olivine

  DOUBLE PRECISION, DIMENSION(3) :: rot
  !! rotation rate vector

  DOUBLE PRECISION, DIMENSION(3,3) :: g
  ! slip tensor

  DOUBLE PRECISION, DIMENSION(3,3) :: lx,ex
  ! dimensionless velocity gradient and strain rate tensors

  double precision :: xtmp,xc1,xc2,xc3,xc4,xc5

  !
  !
!!!! start program
  !
  !
  ! define a few constants
  !
  xc1 = 1.5d0-stressexp
  xc2 = 1.5d0/stressexp
  xc3 = stressexp-1d0
  xc4 = chi/dble(size)
  xc5 = 1.0d0/stressexp

!!! Dimensionless strain rate and velocity gradient tensors
  lx = l/epsnot ; ex = e/epsnot

!!! Plastic deformation + dynamic recrystallization

  DO i=1,size

!!! Calculate invariants e_{pr} T_{pr} for the 
!!! four slip systems of olivine

     bigi=0d0 ; gam = 0d0 ; g = 0d0

     DO i1 = 1,3 ; DO i2 = 1,3
        bigi(1) = bigi(1)+ex(i1,i2)*acs(i,1,i1)*acs(i,2,i2)
        bigi(2) = bigi(2)+ex(i1,i2)*acs(i,1,i1)*acs(i,3,i2)
        bigi(3) = bigi(3)+ex(i1,i2)*acs(i,3,i1)*acs(i,2,i2)
        bigi(4) = bigi(4)+ex(i1,i2)*acs(i,3,i1)*acs(i,1,i2)
     ENDDO; ENDDO

!!! Quotients I/tau
     q = bigi/tau

!!! Reorder quotients I/tau according to absolute magnitude
     qab = ABS(q)
     ti = MAXLOC(qab) ; imax = ti(1) ; qab(imax)=-1d0
     ti = MAXLOC(qab) ; iint = ti(1) ; qab(iint)=-1d0
     ti = MAXLOC(qab) ; imin = ti(1) ; qab(imin)=-1d0
     ti = MAXLOC(qab) ; iinac = ti(1)

!!! Calculate weighting factors gam_s relative to value 
!!! gam_i for which
!!! I/tau is largest

     gam(imax)=1d0

     rat = tau(imax)/bigi(imax)
     qint = rat*bigi(iint)/tau(iint)
     qmin = rat*bigi(imin)/tau(imin)


     gam(iint)=qint*(abs(qint))**xc3
     gam(imin)=qmin*(abs(qmin))**xc3
     gam(iinac)=0d0

!!! calculation of G tensor

     DO i1 = 1,3 ; DO i2 = 1,3
        g(i1,i2)=2d0*(gam(1)*acs(i,1,i1)*acs(i,2,i2) + &
             gam(2)*acs(i,1,i1)*acs(i,3,i2) + &
             gam(3)*acs(i,3,i1)*acs(i,2,i2) + &
             gam(4)*acs(i,3,i1)*acs(i,1,i2))
     END DO; END DO

!!! calculation of strain rate on the softest slip system
!!! (this is fixed to correspond to DREX v2)
     R1 = 0d0 ; R2 = 0d0

     DO j= 1 , 3
        i2 = j + 2
        IF (i2 .GT. 3) i2 = i2 - 3

        R1 = R1 - (g(j,i2)-g(i2,j))*(g(j,i2)-g(i2,j))
        R2 = R2 - (g(j,i2)-g(i2,j))*(lx(j,i2)-lx(i2,j))

        DO k = 1 , 3

           R1 = R1 + 2d0*g(j,k)*g(j,k)
           R2 = R2 + 2d0*lx(j,k)*g(j,k)

        END DO
     END DO

     gam0 = R2/R1

     if(isave_gamma.eq.1)then   ! save the gamma factors for ...
        ioff = (i-1)*3
! those should be local shear stresses normalized by the 
! reference shear stress        
        gamma_save(ioff+1) = ABS(gam(1)*gam0)**xc5
        gamma_save(ioff+2) = ABS(gam(2)*gam0)**xc5
        gamma_save(ioff+3) = ABS(gam(3)*gam0)**xc5
     endif

!!! dislocation density calculation

 
     rt1=tau(1)**xc1*ABS(gam(1)*gam0)**xc2
     rt2=tau(2)**xc1*ABS(gam(2)*gam0)**xc2
     rt3=tau(3)**xc1*ABS(gam(3)*gam0)**xc2

     rt(i) = rt1*exp(-lambda*rt1**2)+ &
             rt2*exp(-lambda*rt2**2)+ &
             rt3*exp(-lambda*rt3**2)

!!! calculation of the rotation rate:
     xtmp = gam0/2.0d0

     rot(3) = (lx(2,1)-lx(1,2))/2d0-(g(2,1)-g(1,2))*xtmp

     rot(2) = (lx(1,3)-lx(3,1))/2d0-(g(1,3)-g(3,1))*xtmp

     rot(1) = (lx(3,2)-lx(2,3))/2d0-(g(3,2)-g(2,3))*xtmp

!!! derivative of the matrix of direction cosine

     dotacs(i,:,:) = 0d0

     DO i1 = 1 , 3 ; DO i2 = 1 , 3 ; DO i3 = 1 , 3 ; DO i4 = 1 , 3
        dotacs(i,i1,i2)=dotacs(i,i1,i2)+&
             alt(i2,i3,i4)*acs(i,i1,i4)*rot(i3)
     END DO; END DO ; END DO ; END DO

!!! grain boundary sliding for small grains
     IF (odf(i) .LT. xc4) THEN
        ! switch off this grain's rotation
        dotacs(i,:,:) = 0d0
        rt(i) = 0d0
     END IF
  END DO

!!! Volume averaged energy
  Emean = SUM(odf*rt)

!!! Change of volume fraction by grain boundary migration
  xtmp = Xol * Mob 
  DO i = 1 , size
     dotodf(i) =  xtmp * odf(i) * (Emean-rt(i))
  END DO

!!! ENSTATITE

  DO i=1,size

        
!!! Calculate slip tensor
        
     g = 0d0
     
     DO i1 = 1,3 ; DO i2 = 1,3
        
        g(i1,i2) = 2d0*acs_ens(i,3,i1)*acs_ens(i,1,i2)
        
     ENDDO; ENDDO
     !! calculation of strain rate

     R1 = 0d0 ; R2 = 0d0
     
     DO j= 1 , 3
        i2 = j + 2
        IF (i2 .GT. 3) i2 = i2 - 3
        
        R1 = R1 - (g(j,i2)-g(i2,j))*(g(j,i2)-g(i2,j))
        R2 = R2 - (g(j,i2)-g(i2,j))*(lx(j,i2)-lx(i2,j))
        
        DO k = 1 , 3
           
           R1 = R1 + 2d0*g(j,k)*g(j,k)
           R2 = R2 + 2d0*lx(j,k)*g(j,k)
           
        END DO
     END DO
     
     gam0 = R2/R1
     
     ! weight factor between olivine and enstatite
     
     gam0 = gam0*(1d0/tau_ens)**stressexp
     
!!! dislocation density calculation
     
     rt0_ens=tau_ens**xc1*ABS(gam0)**xc2
     
     rt_ens(i) = rt0_ens*exp(-lambda*rt0_ens**2)
     
!!! calculation of the rotation rate:
     xtmp = gam0/2d0
        
     rot(3) = (lx(2,1)-lx(1,2))/2d0-(g(2,1)-g(1,2))*xtmp
     
     rot(2) = (lx(1,3)-lx(3,1))/2d0-(g(1,3)-g(3,1))*xtmp
     
     rot(1) = (lx(3,2)-lx(2,3))/2d0-(g(3,2)-g(2,3))*xtmp
        
     
     dotacs_ens(i,:,:) = 0d0

     DO i1 = 1 , 3 ; 
        DO i2 = 1 , 3 ; 
           DO i3 = 1 , 3 ; 
              DO i4 = 1 , 3
                 dotacs_ens(i,i1,i2)=dotacs_ens(i,i1,i2)+&
                      alt(i2,i3,i4)*acs_ens(i,i1,i4)*rot(i3)
              END DO; 
           END DO; 
        END DO; 
     END DO
              
!!! grain boundary sliding for small grains
     IF (odf_ens(i) .LT. xc4) THEN
        dotacs_ens(i,:,:) = 0d0
        rt_ens(i) = 0d0
     END IF

  END DO

!!! Volume averaged energy
  Emean_ens = SUM(odf_ens*rt_ens)

!!! Change of volume fraction by grain boundary migration
  xtmp = (1.0d0-Xol) * Mob
  DO i = 1 , size
     dotodf_ens(i) = xtmp * odf_ens(i) * (Emean_ens-rt_ens(i))
  END DO

!!!
!!! rescale all properties by epsnot
!!!

  dotodf = dotodf * epsnot
  dotodf_ens = dotodf_ens * epsnot
  dotacs = dotacs * epsnot
  dotacs_ens = dotacs_ens * epsnot


  RETURN

END SUBROUTINE drex_deriv_ftrn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
