!
!
!  subroutines for D-REX by E. Kaminski as of 01/09/2003
!
!  minor changes by TWB
!
!  $Id: init.f90,v 1.5 2005/07/25 17:33:01 becker Exp $
!
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! subroutine INIT0, initialization of variables, arrays and parameters   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE init0

   USE comvar

   IMPLICIT NONE

   integer icode,irmode
   CHARACTER(30) :: fname
   ! name of the imput file containing X1, X3, Ui and Dij

   !
   ! initialize the drex-parameters
   !

   icode = 1 ! read para from file
   irmode = 0 ! built in random numbers
   
   call drex_init_para(isize,size3,Xol,&
        tau,tau_ens,Mob,chi,lambda,stressexp,&
        alt,del,ijkl,l1,l2,S0,S0_ens,icode)
   !
   ! allocation of general arrays
   !
   ALLOCATE(odf(isize),dotodf(isize))
   ALLOCATE(odf_ens(isize),dotodf_ens(isize))
   ALLOCATE(rt(isize),rt_ens(isize))
   ALLOCATE(acs(isize,3,3),dotacs(isize,3,3),acs0(isize,3,3))
   ALLOCATE(acs_ens(isize,3,3),dotacs_ens(isize,3,3))
   !
   ! initialize the initial angles of the orientation distribution
   ! function
   !
   call drex_init_acs_random_ftrn(isize,size3,acs0,irmode)
   !
   !
   ! initialize the parameters for the velocity field
   !
   call drex_init_vel_para(nx1,nx3,stepx1,stepx3,&
        i1first,i3first,i1last,i3last,fname)
   !
   ! allocate space for velocity arrays
   !
   ALLOCATE(X1(nx1),X3(nx3),Ui(3,nx1,nx3),Dij(3,3,nx1,nx3))
   ALLOCATE(Fij(3,3,nx1,nx3),grout(nx1,nx3),GOL(nx1,nx3))
   ALLOCATE(phi_fse(nx1,nx3),ln_fse(nx1,nx3),phi_a(nx1,nx3),perc_a(nx1,nx3))
   !
   ! initialize the velocities
   call drex_init_vel(fname,nx1,nx3,x1,x3,Ui,Dij,Fij)





   RETURN


   END SUBROUTINE init0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   


   !
   ! read in parameters for velocity model
   !
   !
   subroutine drex_init_vel_para(nx1,nx3,stepx1,stepx3,&
        i1first,i3first,i1last,i3last,fname)

     IMPLICIT NONE
     INTEGER :: nx1,nx3,stepx1,stepx3,&
          i1first,i3first,i1last,i3last
     CHARACTER(30) :: fname
     
     OPEN(15,file="input.vel.dat")
     ! file and other vel/gradient info
     READ(15,*) fname 
     READ(15,*) nx1,nx3,stepx1,stepx3 
     READ(15,*) i1first,i3first,i1last,i3last
     CLOSE(15)
   end subroutine drex_init_vel_para
   !
   ! read in velocities in velocity gradients
   !
   subroutine drex_init_vel(fname,nx1,nx3,x1,x3,Ui,Dij,Fij)
     
     IMPLICIT NONE
     INTEGER :: i1,i3,nx1,nx3
     CHARACTER(30) :: fname
     DOUBLE PRECISION, DIMENSION(nx1) :: X1
     ! grid points coordinates - distance from the ridge
     DOUBLE PRECISION, DIMENSION(nx3) :: X3
     ! grid points coordinates - depth
     DOUBLE PRECISION, DIMENSION(3,nx1,nx3):: Ui
     ! velocity vector - U1(x1,x3)=Ui(1,x1,x3)
     DOUBLE PRECISION, DIMENSION(3,3,nx1,nx3) :: Dij,Fij
     ! velocity gradient tensor - D11(x1,x3)=Dij(1,1,x1,x3)

!!! Initialization of grid, velocity field and velocity gradient
     
     X1 = 0d0 ; X3 = 0d0 ; Ui = 0d0 ; Dij = 0d0
     
     OPEN(8,file=fname)
   
     DO i1 = 1, nx1
        DO i3 = 1, nx3
           READ(8,1000) X1(i1),X3(i3),Ui(1,i1,i3),Ui(2,i1,i3),Ui(3,i1,i3), &
                Dij(1,1,i1,i3),Dij(1,3,i1,i3),Dij(2,3,i1,i3),Dij(3,1,i1,i3)
           Dij(3,3,i1,i3) = -Dij(1,1,i1,i3)-Dij(2,2,i1,i3)
        END DO
     END DO
1000 FORMAT(16(1pe14.6))

     CLOSE(8)
     !
     ! Initial deformation gradient tensor
     !
     Fij = 0d0 ; Fij(1,1,:,:) = 1d0 ; 
     Fij(2,2,:,:) = 1d0 ; Fij(3,3,:,:) = 1d0
   end subroutine drex_init_vel


