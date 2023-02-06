!!! Module of common variables                                             !!!
!
!
! $Id: comvar.f90,v 1.4 2005/07/25 17:33:01 becker Exp $
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   MODULE comvar

!!! Grid points - Velocity field - Velocity gradient tensor

   INTEGER :: nx1, nx3 !! grid size in X1 and X3 direction
   INTEGER :: stepx1, stepx3 
   !! interval of calculation of the LPO on the grid
   !! LPO calculates only every stepx1 and stepx3 point
   INTEGER :: i1first,i3first,i1last,i3last
   !! indices of the first and last points at which LPO is calculated

   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: X1
   ! grid points coordinates - distance from the ridge
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: X3
   ! grid points coordinates - depth

   DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: Ui
   ! velocity vector - U1(x1,x3)=Ui(1,x1,x3)

   DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: Dij
   ! velocity gradient tensor - D11(x1,x3)=Dij(1,1,x1,x3)

   DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: Fij
   ! finite strain tensor - F11(x1,x3)=Fij(1,1,x1,x3)

   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: grout
   ! dummy used to activate graphic output only for calculated points

!!! LPO calculation

   DOUBLE PRECISION, DIMENSION(3,3) :: l,e
   ! velocity gradient tensor and strain rate tensor

   DOUBLE PRECISION :: epsnot
   ! reference strain rate

!!! Dynamic recrystallization
   
   INTEGER :: size3, isize ! size = size3^3
   ! number of points in the (metric) Eulerian space

   DOUBLE PRECISION :: lambda, Mob, chi
   ! lambda = nucleation parameter
   ! Mob = grain mobility
   ! threshold volume fraction for activation of grain boundary sliding

   DOUBLE PRECISION :: Xol
   ! fraction of olivine in the aggregate

   DOUBLE PRECISION, DIMENSION(3,3,3) :: alt ! \epsilon_{ijk} tensor
   DOUBLE PRECISION, DIMENSION(3,3) :: del ! \delta_{ij} tensor
   INTEGER, DIMENSION(3,3) :: ijkl ! tensor of indices to form Cijkl from Sij
   INTEGER, DIMENSION(6) :: l1,l2 ! tensot of indices to form Sij from Cijkl

   DOUBLE PRECISION, DIMENSION(4) :: tau
   ! RSS for the 4 slip systems of olivine (only 3 are activated)
   DOUBLE PRECISION :: tau_ens ! RSS of enstatite slip system
   DOUBLE PRECISION :: stressexp ! stress exponent for olivine and enstatite

   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: odf,dotodf
   ! volume fraction of the olivine grains
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: odf_ens,dotodf_ens
   ! volume fraction of the enstatite grains

   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: rt,rt_ens
   ! dislocations density for olivine and enstatite

   DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: acs,dotacs,acs0
   DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: acs_ens,dotacs_ens
   !! matrix of direction cosine

!!! Output

   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: phi_fse
   ! orientation of the long axis of the FSE on the (X1,X3) plan

   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ln_fse
   ! finite strain = ln(longaxis/shortaxis)

   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: phi_a
   ! average orientation of a-axis on (X1,X3) plan

   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: perc_a
   ! percentage of S wave anisotropy

   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: GOL
   ! Grain Orientation Lag i.e. Pi parameter

!!! Elastic tensor 

   DOUBLE PRECISION, DIMENSION(6,6) :: S0,Sav,S0_ens
   ! stiffness matrix
   DOUBLE PRECISION, DIMENSION(3,3,3,3) :: Cav
   !Cijkl tensor at point of calculation

   END MODULE comvar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

