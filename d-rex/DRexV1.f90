
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Paris 01/11/2003 - Calculates the evolutions of texture an olivine     !!!
!!! aggregate in a ridge+shear flow (flow from Neil Ribe).                 !!!
!!! LPO evolves by dynamic recrystallization - Calculated using D-Rex      !!!
!!! Equivalent transverse isotropic media calculated using the theory of   !!!
!!! Jules Browaeys. A couple of bugs fixed thanks to Teresa Lassak         !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
! this version slightly modified by TWB
!
! $Id: DRexV1.f90,v 1.18 2011/04/12 06:18:38 becker Exp becker $
!
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Main program                                                           !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   PROGRAM DR_CF 

   USE comvar
   use drex_nrmod
   IMPLICIT NONE

   DOUBLE PRECISION :: X10,X30 
   !! coordinates of the grid point where calculation is started

   INTEGER :: i1,i3 !! loop counters
   INTEGER :: i1a,i3a !! indices of the grid point at the end of streamline

   INTEGER :: step1 !! number of points on a streamline

   DOUBLE PRECISION, DIMENSION(3,3) :: LSij
   ! left-strech tensor for FSE calculation

   INTEGER :: nrot
   ! number of rotations for the Jacobi

   DOUBLE PRECISION, DIMENSION(3) :: evals
   DOUBLE PRECISION, DIMENSION(3,3) :: evects,scc_rot
   ! eigen values and vectors in jacobi
   
   double precision, dimension (3) :: tiaxis
   ! symmetry axis of approximated hexagonal tensor

   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: X1gr,X3gr
   !! actual coordinates of calculation point when not a grid point

   double precision :: k,g,favg
   double precision, dimension(9) :: vel ! various velocities from tensor
   double precision, dimension(6,6) :: dummyt
   double precision,dimension(6) :: symm_frac

   ! averaging is voigt
   favg = 0.0d0
   
!!! initialization

   CALL init0

   ALLOCATE(X1gr(nx1),X3gr(nx3))

   
   OPEN(18,FILE='Cijkl',form='unformatted')

   X1gr = 0d0 ; X3gr = 0d0
   phi_fse = 0d0 ; ln_fse = 0d0 ; phi_a = 0d0 ; perc_a = 0d0
   grout = 0d0 ; GOL = 0d0

!!! Calculation of GOL parameter

   DO i1 = 1, nx1-1
      DO i3 = 1 , nx3-1

!!! Location of the calculation point 
         X10 = (X1(i1)+X1(i1+1))/2d0 ; X30 = (X3(i3)+X3(i3+1))/2d0
         X1gr(i1) = X10 ; X3gr(i3) = X30

         CALL pipar(i1,i3,X1gr(i1),X3gr(i3))

      END DO
   END DO

!   print *,'pipar ok'
!!! FSE and LPO calculation

! loop through all points
   DO i1 = i1first , i1last-1 , stepx1 
      DO i3 = i3first , i3last-1 , stepx3

!!! Initial location of the grid point 
         X10 = (X1(i1)+X1(i1+1))/2d0 
         X30 = (X3(i3)+X3(i3+1))/2d0
         X1gr(i1) = X10 
         X3gr(i3) = X30

!!! Backward calculation of the pathline for each tracer

         step1 = 0
         CALL pathline(X10,X30,i1,i3,step1,i1a,i3a)

!!! Inward calculation of the LPO
!
!!! Random initial LPO
!
         call drex_init_pathline_ftrn(rt,odf,acs,acs0,isize)
         call drex_init_pathline_ftrn(rt_ens,odf_ens,&
              acs_ens,acs0,isize)


!
! advect material using Runge Kutta
         CALL strain(X10,X30,i1a,i3a,step1,Fij(:,:,i1,i3))

!         print *,'strain ok'
!!! Left-stretch tensor for FSE calculation

         LSij = MATMUL(Fij(:,:,i1,i3),TRANSPOSE(Fij(:,:,i1,i3)))

! find eigensystem of left stretch tensor
         CALL drex_nr_jacobi(LSij,evals,evects,nrot)

!         print *,'Jacobi ok'
!!! Pick up the orientation of the long axis of the FSE

         IF (evals(1) .EQ. MAXVAL(evals)) THEN
            phi_fse(i1,i3) = ATAN2(evects(3,1),evects(1,1))
         ELSE IF(evals(2) .EQ. MAXVAL(evals)) THEN
            phi_fse(i1,i3) = ATAN2(evects(3,2),evects(1,2))
         ELSE
            phi_fse(i1,i3) = ATAN2(evects(3,3),evects(1,3))
         END IF

!!! natural strain = ln(a/c) where a is the 
!!! long axis = maxval(evals)**0.5
         ln_fse(i1,i3) = 0.5*LOG(MAXVAL(evals)/MINVAL(evals))
!!! Activate graphic output
         grout(i1,i3) = 1d0
!!!
!!! Cijkl tensor (using Voigt average)
         CALL drex_vhr_avg_ftrn(S0,S0_ens,acs,acs_ens,odf,odf_ens,&
              Xol,isize,ijkl,l1,l2,Sav,favg)
         WRITE(18) X1gr(i1),X3gr(i3),Sav
!!
!!! Percentage of anisotropy and orientation of 
!!! axis of hexagonal symmetry
         CALL drex_decsym_ftrn(Sav,k,g,symm_frac,tiaxis,dummyt,dummyt,dummyt,&
              dummyt,dummyt,dummyt,vel,dummyt,scc_rot,1)
         perc_a(i1,i3) = symm_frac(2) * 100d0 ! VTI percent
!!! inclination angle with third axis
         phi_a(i1,i3) = asin(tiaxis(3))
      END DO
   END DO

   CLOSE(18)

!!! Graphic output
   
   CALL check_veloc

   OPEN(10,file='fse.m',status='replace')
   OPEN(11,file='a_axis.m',status='replace')
   OPEN(12,file='GOL.m',status='replace')

   WRITE(10,*)
   WRITE(10,*) "clear"
   WRITE(10,*) "figure"

   WRITE(11,*)
   WRITE(11,*) "clear"
   WRITE(11,*) "figure"

   WRITE(12,*)
   WRITE(12,*) "clear"
   WRITE(12,*) "figure"

   WRITE(10,*) "Ux=["
   DO i3 = 1,nx3
      WRITE(10,"(81(1PE12.4))") grout(:,i3)*ln_fse(:,i3)*COS(phi_fse(:,i3))
   END DO
   WRITE(10,*) "];"

   WRITE(10,*) "Uz=["
   DO i3 = 1,nx3
      WRITE(10,"(81(1PE12.4))") grout(:,i3)*ln_fse(:,i3)*SIN(phi_fse(:,i3))
   END DO
   WRITE(10,*) "];"

   WRITE(11,*) "Ux=["
   DO i3 = 1,nx3
      WRITE(11,"(81(1PE12.4))") grout(:,i3)*COS(phi_a(:,i3))*perc_a(:,i3)
   END DO
   WRITE(11,*) "];"

   WRITE(11,*) "Uz=["
   DO i3 = 1,nx3
      WRITE(11,"(81(1PE12.4))") grout(:,i3)*SIN(phi_a(:,i3))*perc_a(:,i3)
   END DO
   WRITE(11,*) "];"

   WRITE(12,*) "Ux=["
   DO i3 = 1,nx3
      WRITE(12,"(81(1PE12.4))") GOL(:,i3)
   END DO
   WRITE(12,*) "];"

   WRITE(10,*) "X=["
   WRITE(10,"(81(1PE12.4))") X1gr
   WRITE(10,*) "];"

   WRITE(10,*) "Z=["
   WRITE(10,"(81(1PE12.4))") X3gr
   WRITE(10,*) "];"

   WRITE(11,*) "X=["
   WRITE(11,"(81(1PE12.4))") X1gr
   WRITE(11,*) "];"

   WRITE(11,*) "Z=["
   WRITE(11,"(81(1PE12.4))") X3gr
   WRITE(11,*) "];"

   WRITE(12,*) "X=["
   WRITE(12,"(81(1PE12.4))") X1gr
   WRITE(12,*) "];"

   WRITE(12,*) "Z=["
   WRITE(12,"(81(1PE12.4))") X3gr
   WRITE(12,*) "];"

   WRITE(10,*) "quiver(X,Z,Ux,Uz,'.')"
   WRITE(10,*) "axis([-4 4 0 3])"
   WRITE(10,*) "xlabel('X1, dimensionless distance from the ridge axis')"
   WRITE(10,*) "ylabel('X3, dimensionless depth')"
   WRITE(10,*) "title('Orientation of the long axis of the FSE')"
   
   WRITE(10,*) "axis square"
   WRITE(10,*) "axis ij"
   
   WRITE(11,*) "quiver(X,Z,Ux,Uz,'.')"
   WRITE(11,*) "axis([-4 4 0 3])"
   WRITE(11,*) "xlabel('X1, dimensionless distance from the ridge axis')"
   WRITE(11,*) "ylabel('X3, dimensionless depth')"
   WRITE(11,*) "title('Orientation of the hexagonal symmetry axis')"
   
   WRITE(11,*) "axis square"
   WRITE(11,*) "axis ij"

   WRITE(12,*) "surf(X,Z,log(abs(Ux)))"
   WRITE(12,*) "axis([-4 4 0 3])"
   WRITE(12,*) "caxis([-6 3])"
   WRITE(12,*) "shading 'interp'"
   WRITE(12,*) "colorbar 'horiz'"
   WRITE(12,*) "xlabel('X1, dimensionless distance from the ridge axis')"
   WRITE(12,*) "ylabel('X3, dimensionless depth')"
   WRITE(12,*) "title('Values of GOL parameter (in log units)')"
   
   WRITE(12,*) "axis square"
   WRITE(12,*) "axis ij"

   CLOSE(10)
   CLOSE(11)
   CLOSE(12)

   STOP

   END PROGRAM DR_CF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine PATHLINE - Calculation of tracers pathlines                 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE pathline(X1s,X3s,i1s,i3s,step2,i1r,i3r)

   USE comvar 

   IMPLICIT NONE

   INTEGER :: step2
   ! number of steps to construct the streamline

   INTEGER :: i1s,i3s,i1,i3,i1r,i3r
   ! initial, intermediate and final indices of the closest upper-left grid point

   DOUBLE PRECISION :: ds,dt
   ! elementary displacement along the streamline and corresponding time step

   DOUBLE PRECISION :: X1s,X3s,X1i,X3i
   ! initial and intermediate position of the tracer

   DOUBLE PRECISION :: U1i,U3i
   ! intermediate velocity components on the streamline

   DOUBLE PRECISION :: kx1,kx2,kx3,kx4
   DOUBLE PRECISION :: kz1,kz2,kz3,kz4
   ! RGK intermediate for streamlines calculation 

!! initial indices and position of the tracers

   i1 = i1s ; i3 = i3s

!!! construction of the streamline

   DO WHILE (X3s .LT. X3(nx3) .AND. X1s .GT. X1(1))

      X1i = X1s ; X3i = X3s
      step2 = step2 + 1

      CALL velocitycalc(X1i,X3i,U1i,U3i,i1,i3)
      CALL gradientcalc(X1i,X3i,i1,i3)

      ds = 0.075d0*MIN(ABS(X1(i1+1)-X1(i1)),ABS(X3(i3+1)-X3(i3)))
      dt = MIN(ds/SQRT(U1i**2+U3i**2),1d-2/epsnot)

      kx1 = -U1i*dt ; kz1 = -U3i*dt

      X1i = X1s + 0.5d0*kx1
      X3i = X3s + 0.5d0*kz1

!!! find the UPPER-LEFT grid point closest to the point of calculation
      DO WHILE (X1(i1) .GT. X1i .AND. i1 .GT. 1) ; i1 = i1-1 ; ENDDO
      DO WHILE (X1(i1+1) .LT. X1i .AND. i1 .LT. nx1-1) ; i1 = i1+1 ; ENDDO

      DO WHILE (X3(i3) .GT. X3i .AND. i3 .GT. 1) ; i3 = i3-1 ; ENDDO
      DO WHILE (X3(i3+1) .LT. X3i .AND. i3 .LT. nx3-1) ; i3 = i3+1 ; ENDDO

      CALL velocitycalc(X1i,X3i,U1i,U3i,i1,i3)

      kx2 = -U1i*dt ; kz2 = -U3i*dt

      X1i = X1s + 0.5d0*kx2
      X3i = X3s + 0.5d0*kz2

!!! find the UPPER-LEFT grid point closest to the point of calculation
      DO WHILE (X1(i1) .GT. X1i .AND. i1 .GT. 1) ; i1 = i1-1 ; ENDDO
      DO WHILE (X1(i1+1) .LT. X1i .AND. i1 .LT. nx1-1) ; i1 = i1+1 ; ENDDO

      DO WHILE (X3(i3) .GT. X3i .AND. i3 .GT. 1) ; i3 = i3-1 ; ENDDO
      DO WHILE (X3(i3+1) .LT. X3i .AND. i3 .LT. nx3-1) ; i3 = i3+1 ; ENDDO

      CALL velocitycalc(X1i,X3i,U1i,U3i,i1,i3)

      kx3 = -U1i*dt ; kz3 = -U3i*dt

      X1i = X1s + kx3
      X3i = X3s + kz3

!!! find the UPPER-LEFT grid point closest to the point of calculation
      DO WHILE (X1(i1) .GT. X1i .AND. i1 .GT. 1) ; i1 = i1-1 ; ENDDO
      DO WHILE (X1(i1+1) .LT. X1i .AND. i1 .LT. nx1-1) ; i1 = i1+1 ; ENDDO

      DO WHILE (X3(i3) .GT. X3i .AND. i3 .GT. 1) ; i3 = i3-1 ; ENDDO
      DO WHILE (X3(i3+1) .LT. X3i .AND. i3 .LT. nx3-1) ; i3 = i3+1 ; ENDDO

      CALL velocitycalc(X1i,X3i,U1i,U3i,i1,i3)

      kx4 = -U1i*dt ; kz4 = -U3i*dt

      X1s = X1s + (kx1/2d0+kx2+kx3+kx4/2d0)/3d0
      X3s = X3s + (kz1/2d0+kz2+kz3+kz4/2d0)/3d0

!!! find the UPPER-LEFT grid point closest to the point of calculation
      DO WHILE (X1(i1) .GT. X1s .AND. i1 .GT. 1) ; i1 = i1-1 ; ENDDO
      DO WHILE (X1(i1+1) .LT. X1s .AND. i1 .LT. nx1-1) ; i1 = i1+1 ; ENDDO

      DO WHILE (X3(i3) .GT. X3s .AND. i3 .GT. 1) ; i3 = i3-1 ; ENDDO
      DO WHILE (X3(i3+1) .LT. X3s .AND. i3 .LT. nx3-1) ; i3 = i3+1 ; ENDDO

      i1r = i1 ; i3r = i3

   END DO

   RETURN

   END SUBROUTINE pathline

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine VELOCITYCALC, calculation of velocity at a given point      !!!
!!! by interpolation method given in Numerical Recipies, Press et al., p96 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE velocitycalc(X1s,X3s,U1s,U3s,i1s,i3s)

   USE comvar

   IMPLICIT NONE

   DOUBLE PRECISION :: X1s,X3s
   ! x and y coordinates of the point on the streamline

   DOUBLE PRECISION :: U1s,U3s
   ! interpolated velocity at the point

   INTEGER :: i1s,i3s
   ! indices of the UP-LEFT grid point closest to the extrapolation point

   DOUBLE PRECISION :: y1,y2,y3,y4
   ! dummies for interpolation

   y1 = Ui(1,i1s,i3s) ; y2 = Ui(1,i1s+1,i3s)
   y3 = Ui(1,i1s+1,i3s+1) ; y4 = Ui(1,i1s,i3s+1)

   CALL interp(X1s,X3s,i1s,i3s,y1,y2,y3,y4,U1s)

   y1 = Ui(3,i1s,i3s) ; y2 = Ui(3,i1s+1,i3s)
   y3 = Ui(3,i1s+1,i3s+1) ; y4 = Ui(3,i1s,i3s+1)

   CALL interp(X1s,X3s,i1s,i3s,y1,y2,y3,y4,U3s)

   RETURN

   END SUBROUTINE velocitycalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine interp, bilinear interpolation of a function                !!!
!!! following the method given in Numerical Recipies, Press et al., p96    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE interp(X1s,X3s,i1s,i3s,y1,y2,y3,y4,res)

   USE comvar

   IMPLICIT NONE

   DOUBLE PRECISION :: X1s,X3s
   ! x and y coordinates of the point on the streamline

   DOUBLE PRECISION :: res
   ! result of the interpolation at the point

   INTEGER :: i1s,i3s
   ! indices of the UP-LEFT grid point closest to the considered point

   DOUBLE PRECISION :: y1,y2,y3,y4,tt,uu
   ! dummies for interpolation (numerical recipies, p96)

   tt = (X1s-X1(i1s))/(X1(i1s+1)-X1(i1s))
   uu = (X3s-X3(i3s))/(X3(i3s+1)-X3(i3s))

   res = (1d0-tt)*(1d0-uu)*y1+tt*(1d0-uu)*y2+tt*uu*y3+(1d0-tt)*uu*y4

   RETURN

   END SUBROUTINE interp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine GRADIENTCALC, interpolation of 
!!! velocity gradient tensor     !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE gradientcalc(X1s,X3s,i1s,i3s)

   USE comvar

   IMPLICIT NONE

   DOUBLE PRECISION :: X1s,X3s
   ! x and y coordinates of the point on the streamline

   INTEGER :: i1s,i3s
   ! indices of the UPPER-LEFT grid point closest to the considered point

   DOUBLE PRECISION :: y1,y2,y3,y4
   ! dummies for interpolation

   
   ! fill in l, the velocity gradient tensor

   l = 0d0

   y1 = Dij(1,1,i1s,i3s) ; y2 = Dij(1,1,i1s+1,i3s)
   y3 = Dij(1,1,i1s+1,i3s+1) ; y4 = Dij(1,1,i1s,i3s+1)

   CALL interp(X1s,X3s,i1s,i3s,y1,y2,y3,y4,l(1,1))

   y1 = Dij(1,3,i1s,i3s) ; y2 = Dij(1,3,i1s+1,i3s)
   y3 = Dij(1,3,i1s+1,i3s+1) ; y4 = Dij(1,3,i1s,i3s+1)

   CALL interp(X1s,X3s,i1s,i3s,y1,y2,y3,y4,l(1,3))

   y1 = Dij(3,1,i1s,i3s) ; y2 = Dij(3,1,i1s+1,i3s)
   y3 = Dij(3,1,i1s+1,i3s+1) ; y4 = Dij(3,1,i1s,i3s+1)

   CALL interp(X1s,X3s,i1s,i3s,y1,y2,y3,y4,l(3,1))

   y1 = Dij(2,3,i1s,i3s) ; y2 = Dij(2,3,i1s+1,i3s)
   y3 = Dij(2,3,i1s+1,i3s+1) ; y4 = Dij(2,3,i1s,i3s+1)

   CALL interp(X1s,X3s,i1s,i3s,y1,y2,y3,y4,l(2,3))
!
! same trick as in our code! (TWB comment)
   l(3,3) = -l(1,1)-l(2,2)

   !
   ! compute strain rate tensor and characteristic strain rate 
   !
   call drex_strain_rate_ftrn(l,e,epsnot)

   RETURN

   END SUBROUTINE gradientcalc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine STRAIN - Calculation of strain along pathlines              !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE strain(X1s,X3s,i1s,i3s,step2,fse)

   USE comvar 

   IMPLICIT NONE


   INTEGER :: step2
   ! number of steps used to construct the streamline

   INTEGER :: i1s,i3s,i1,i3
   ! initial and intermediate indices of the closest upper-left grid point

   DOUBLE PRECISION :: ds,dt
   ! elementary displacement along the streamline and corresponding time step

   DOUBLE PRECISION :: X1s,X3s,X1i,X3i
   ! initial and intermediate position of the tracer

   DOUBLE PRECISION :: U1i,U3i
   ! intermediate velocity components for the tracer

   DOUBLE PRECISION, DIMENSION(3,3) :: fse,fsei
   ! local finite deformation tensor

   DOUBLE PRECISION :: kx1,kx2,kx3,kx4
   DOUBLE PRECISION :: kz1,kz2,kz3,kz4
   ! RGK intermediate for streamlines calculation

   DOUBLE PRECISION, DIMENSION(3,3) :: kfse1,kfse2,kfse3,kfse4
   ! RGK intermediates for Finite Strain Ellipsoid

   DOUBLE PRECISION, DIMENSION(isize) :: kodf1,kodf2,kodf3,kodf4
   DOUBLE PRECISION, DIMENSION(isize) :: kodf1_ens,kodf2_ens,kodf3_ens,kodf4_ens
   ! intermediates for odf

   DOUBLE PRECISION, DIMENSION(isize,3,3) :: kac1,kac2,kac3,kac4
   DOUBLE PRECISION, DIMENSION(isize,3,3) :: kac1_ens,kac2_ens,kac3_ens,kac4_ens
   ! intermediates for matrix of direction cosine
   
   double precision, dimension (:,:,:), allocatable :: acsi,acsi_ens
   double precision, dimension (:), allocatable :: odfi,odfi_ens

   integer :: izero
   real :: fdummy

   !
   ! only used locally
   !
   allocate(acsi(isize,3,3),acsi_ens(isize,3,3))
   allocate(odfi(isize),odfi_ens(isize))


   izero = 0
   


!! initial indices and position of the tracers

   i1 = i1s ; i3 = i3s

!! LPO calculation along the streamline

   DO WHILE (step2 .GT. 0)

      X1i = X1s ; X3i = X3s
      fsei = fse
      odfi = odf ; acsi = acs
      odfi_ens = odf_ens ; acsi_ens = acs_ens
      step2 = step2 - 1

      CALL velocitycalc(X1i,X3i,U1i,U3i,i1,i3)
! this computes l, e, and epsnot
      CALL gradientcalc(X1i,X3i,i1,i3)
! this gets the derivatives
      CALL drex_deriv_ftrn(l,e,epsnot,isize,acsi,acsi_ens,dotacs, &
           dotacs_ens,odfi,odfi_ens,dotodf,dotodf_ens,&
           rt,rt_ens,tau,tau_ens,stressexp,alt, &
           Xol,lambda,Mob,chi,izero,fdummy)

      ds = 0.075d0*MIN(ABS(X1(i1+1)-X1(i1)),ABS(X3(i3+1)-X3(i3)))
      dt = MIN(ds/SQRT(U1i**2+U3i**2),1d-2/epsnot)
! vel
      kx1 = U1i*dt ; kz1 = U3i*dt
! FSE
      kfse1 = MATMUL(l,fsei)*dt

! odf
      kodf1 = dotodf*dt
      kodf1_ens = dotodf_ens*dt
! acs
      kac1 = dotacs*dt
      kac1_ens = dotacs_ens*dt

      X1i = X1s + 0.5d0*kx1
      X3i = X3s + 0.5d0*kz1
      fsei = fse + 0.5d0*kfse1
      odfi = odf + 0.5d0*kodf1
      acsi = acs + 0.5d0*kac1
      odfi_ens = odf_ens + 0.5d0*kodf1_ens
      acsi_ens = acs_ens + 0.5d0*kac1_ens

      call drex_check_phys_limits_ftrn(acsi,acsi_ens,odfi,odfi_ens,isize)

!!! find the UPPER-LEFT grid point closest to the point of calculation
      DO WHILE (X1(i1) .GT. X1i .AND. i1 .GT. 1) ; i1 = i1-1 ; ENDDO
      DO WHILE (X1(i1+1) .LT. X1i .AND. i1 .LT. nx1-1) ; i1 = i1+1 ; ENDDO

      DO WHILE (X3(i3) .GT. X3i .AND. i3 .GT. 1) ; i3 = i3-1 ; ENDDO
      DO WHILE (X3(i3+1) .LT. X3i .AND. i3 .LT. nx3-1) ; i3 = i3+1 ; ENDDO

      CALL velocitycalc(X1i,X3i,U1i,U3i,i1,i3)
      CALL gradientcalc(X1i,X3i,i1,i3)
      CALL drex_deriv_ftrn(l,e,epsnot,isize,acsi,acsi_ens,dotacs, &
           dotacs_ens,odfi,odfi_ens,dotodf,dotodf_ens,&
           rt,rt_ens,tau,tau_ens,stressexp,alt, &
           Xol,lambda,Mob,chi,izero,fdummy)
! velocities
      kx2 = U1i*dt 
      kz2 = U3i*dt
! finite strain
      kfse2 = MATMUL(l,fsei)*dt

! ODF
      kodf2 =     dotodf    *dt
      kodf2_ens = dotodf_ens*dt
! ACS
      kac2 =     dotacs    *dt
      kac2_ens = dotacs_ens*dt

      X1i = X1s + 0.5d0*kx2
      X3i = X3s + 0.5d0*kz2
      fsei = fse + 0.5d0*kfse2
      odfi = odf + 0.5d0*kodf2
      odfi_ens = odf_ens + 0.5d0*kodf2_ens
      acsi = acs + 0.5d0*kac2
      acsi_ens = acs_ens + 0.5d0*kac2_ens

      call drex_check_phys_limits_ftrn(acsi,acsi_ens,odfi,odfi_ens,isize)

!!! find the UPPER-LEFT grid point closest to the point of calculation
      DO WHILE (X1(i1) .GT. X1i .AND. i1 .GT. 1) ; i1 = i1-1 ; ENDDO
      DO WHILE (X1(i1+1) .LT. X1i .AND. i1 .LT. nx1-1) ; i1 = i1+1 ; ENDDO

      DO WHILE (X3(i3) .GT. X3i .AND. i3 .GT. 1) ; i3 = i3-1 ; ENDDO
      DO WHILE (X3(i3+1) .LT. X3i .AND. i3 .LT. nx3-1) ; i3 = i3+1 ; ENDDO

      CALL velocitycalc(X1i,X3i,U1i,U3i,i1,i3)
      CALL gradientcalc(X1i,X3i,i1,i3)
      CALL drex_deriv_ftrn(l,e,epsnot,isize,acsi,acsi_ens,dotacs, &
           dotacs_ens,odfi,odfi_ens,dotodf,dotodf_ens,&
           rt,rt_ens,tau,tau_ens,stressexp,alt, &
           Xol,lambda,Mob,chi,izero,fdummy)

      kx3 = U1i*dt ; kz3 = U3i*dt
      kfse3 = MATMUL(l,fsei)*dt

      kodf3 = dotodf*dt
      kodf3_ens = dotodf_ens*dt
      kac3 = dotacs*dt
      kac3_ens = dotacs_ens*dt

      X1i = X1s + kx3
      X3i = X3s + kz3
      fsei = fse + kfse3
      odfi = odf + kodf3
      odfi_ens = odf_ens + kodf3_ens
      acsi = acs + kac3
      acsi_ens = acs_ens + kac3_ens

      call drex_check_phys_limits_ftrn(acsi,acsi_ens,&
           odfi,odfi_ens,isize)

!!! find the UPPER-LEFT grid point closest to the point of calculation
      DO WHILE (X1(i1) .GT. X1i .AND. i1 .GT. 1) ; i1 = i1-1 ; ENDDO
      DO WHILE (X1(i1+1) .LT. X1i .AND. i1 .LT. nx1-1) ; i1 = i1+1 ; ENDDO

      DO WHILE (X3(i3) .GT. X3i .AND. i3 .GT. 1) ; i3 = i3-1 ; ENDDO
      DO WHILE (X3(i3+1) .LT. X3i .AND. i3 .LT. nx3-1) ; i3 = i3+1 ; ENDDO

      CALL velocitycalc(X1i,X3i,U1i,U3i,i1,i3)
      CALL gradientcalc(X1i,X3i,i1,i3)
      CALL drex_deriv_ftrn(l,e,epsnot,isize,acsi,acsi_ens,dotacs, &
           dotacs_ens,odfi,odfi_ens,dotodf,dotodf_ens,&
           rt,rt_ens,tau,tau_ens,stressexp,alt, &
           Xol,lambda,Mob,chi,izero,fdummy)

      kx4 = U1i*dt ; kz4 = U3i*dt
      kfse4 = MATMUL(l,fsei)*dt

      kodf4 = dotodf*dt
      kodf4_ens = dotodf_ens*dt
      kac4 = dotacs*dt
      kac4_ens = dotacs_ens*dt

      X1s = X1s + (kx1/2d0+kx2+kx3+kx4/2d0)/3d0
      X3s = X3s + (kz1/2d0+kz2+kz3+kz4/2d0)/3d0
      fse = fse + (kfse1/2d0+kfse2+kfse3+kfse4/2d0)/3d0
      acs = acs + (kac1/2d0+kac2+kac3+kac4/2d0)/3d0
      acs_ens = acs_ens + (kac1_ens/2d0+kac2_ens+kac3_ens+kac4_ens/2d0)/3d0
      odf = odf + (kodf1/2d0+kodf2+kodf3+kodf4/2d0)/3d0
      odf_ens = odf_ens + (kodf1_ens/2d0+kodf2_ens+kodf3_ens+kodf4_ens/2d0)/3d0

      call drex_check_phys_limits_ftrn(acs,acs_ens,&
           odf,odf_ens,isize)


!!! find the UPPER-LEFT grid point closest to the point of calculation
      DO WHILE (X1(i1) .GT. X1s .AND. i1 .GT. 1) ; i1 = i1-1 ; ENDDO
      DO WHILE (X1(i1+1) .LT. X1s .AND. i1 .LT. nx1-1) ; i1 = i1+1 ; ENDDO

      DO WHILE (X3(i3) .GT. X3s .AND. i3 .GT. 1) ; i3 = i3-1 ; ENDDO
      DO WHILE (X3(i3+1) .LT. X3s .AND. i3 .LT. nx3-1) ; i3 = i3+1 ; ENDDO

   END DO

   deallocate(odfi,odfi_ens,acsi,acsi_ens)

   RETURN

   END SUBROUTINE strain

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine check_veloc - matlab plot of velocity field                 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   SUBROUTINE check_veloc

   USE comvar

   IMPLICIT NONE
   
   INTEGER :: i3 ! loop counters

   OPEN(10,file='veloc.m',status='replace')

   WRITE(10,*)
   WRITE(10,*) "clear"
   WRITE(10,*) "figure"

   WRITE(10,*) "Ux=["
   DO i3 = 1,nx3
      WRITE(10,"(81(1PE12.4))") Ui(1,:,i3)
   END DO
   WRITE(10,*) "];"
   
   WRITE(10,*) "Uz=["
   DO i3 = 1,nx3
      WRITE(10,"(81(1PE12.4))") Ui(3,:,i3)
   END DO
   WRITE(10,*) "];"

   WRITE(10,*) "X=["
   WRITE(10,"(81(1PE12.4))") X1
   WRITE(10,*) "];"

   WRITE(10,*) "Z=["
   WRITE(10,"(81(1PE12.4))") X3
   WRITE(10,*) "];"

   WRITE(10,*) "quiver(X,Z,Ux,Uz)"
   WRITE(10,*) "axis([-4 4 0 3])"
   WRITE(10,*) "xlabel('X1, dimensionless distance from the ridge axis)')"
   WRITE(10,*) "ylabel('X3, dimensionless depth')"
   WRITE(10,*) "title('Velocity field')"

   WRITE(10,*) "axis square"
   WRITE(10,*) "axis ij"

   CLOSE(10)

   RETURN

   END SUBROUTINE check_veloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Subroutine pipar - Calculates GOL parameter at grid point              !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE pipar(i1,i3,X10,X30)

   USE comvar

   IMPLICIT NONE

   INTEGER :: i1,i3
   ! coordinates of the grid point

   DOUBLE PRECISION :: X10,X30,X1i,X3i
   ! initial and intermediate position of the point

   DOUBLE PRECISION :: U1i,U3i
   ! local velocity vector

   DOUBLE PRECISION, DIMENSION(3) :: veloc,isa
   ! local velocity vector and orientation of infinite strain axis

   DOUBLE PRECISION :: dt
   ! time step used to calculate GOL

   DOUBLE PRECISION :: thetaISA
   ! angle between ISA and flow direction

!!! RMQ if GOL not defined then GOL = -1d0 - 
!!! Max value of GOL set to 10
   GOL(i1,i3) = 1d1

   veloc = 0d0 ; isa = 0d0

   CALL velocitycalc(X10,X30,U1i,U3i,i1,i3)
   dt = 1d-2/SQRT(U1i**2+U3i**2)

!!! previous point on the streamline

   X1i = X10-dt*U1i
   X3i = X30-dt*U3i

   CALL velocitycalc(X1i,X3i,U1i,U3i,i1,i3)
   veloc(1) = U1i ; veloc(3) = U3i ; 
   veloc = veloc/SQRT(SUM(veloc**2))

!!! this computes the gradient, l, and the max eigenvalue
   CALL gradientcalc(X1i,X3i,i1,i3)
   l = l/epsnot

!!! calculation of the ISA
   CALL drex_isacalc_ftrn(isa,gol(i1,i3),l)
   IF (SUM(isa) .EQ. -3d0) isa=veloc

!!! angle between ISA and flow direction
   thetaISA = ACOS(SUM(veloc*isa))

!!! next point on the streamline

   X1i = X10+dt*U1i
   X3i = X30+dt*U3i

   CALL velocitycalc(X1i,X3i,U1i,U3i,i1,i3)
   veloc(1) = U1i ; veloc(3) = U3i ; 
   veloc = veloc/SQRT(SUM(veloc**2))

   CALL gradientcalc(X1i,X3i,i1,i3)
   l = l/epsnot

!!! calculation of the ISA
   CALL drex_isacalc_ftrn(isa,gol(i1,i3),l)
   IF (SUM(isa) .EQ. -3d0) isa=veloc

!!! Pi parameter
   CALL gradientcalc(X10,X30,i1,i3) ! why is this here?

   thetaISA = ABS(thetaISA-ACOS(SUM(veloc*isa)))
   IF (thetaISA .GT. ACOS(-1d0)) thetaISA = thetaISA-ACOS(-1d0)

   GOL(i1,i3) = MIN(GOL(i1,i3),thetaISA/2d0/dt/epsnot)

   RETURN

   END SUBROUTINE pipar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




