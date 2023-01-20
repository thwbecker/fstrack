! elastic_MM.f90
! Jules Browaeys
! July 2005
!
! Calculates the average elastic tensor of n minerals
! n < nmax
!
! ELASTIC_INPUT (subroutine) : input of elastic constants   

program ELASTIC_MM

use NRMOD
use TENSMOD
use DECMOD

implicit none




double precision, dimension(21) :: XH,XD
double precision, dimension(3,3) :: RT


integer, parameter :: nmax = 10,  iflout = 70
integer :: nmineral,n,i,nrot,j,k
double precision :: pi
double precision, dimension(6) :: KVAL
double precision, dimension(6,6) :: CM,KELM,VKEL,cmr
double precision, dimension(nmax) :: xvol
double precision, dimension(3,nmax) :: angle
double precision, dimension(3) :: rangle
double precision, dimension(6,6,nmax) :: CI
character(len=25), dimension(nmax) :: file_in
character(len=25) :: file_out

pi = acos(-1D0)
read(*,*) nmineral
xvol = 0D0
angle = 0D0
CI = 0D0
do n=1,nmineral
 read(*,'(a25)') file_in(n)
! write(*,*) file_in(n)
 read(*,*) xvol(n)
 read(*,*) (angle(i,n),i=1,3)
end do
read(*,'(a25)') file_out
!write(*,*) ' sum xvol =',sum(xvol)

do n=1,nmineral
 call ELASTIC_INPUT(n)
! call C66_OUT(CI(:,:,n))
end do
!write(*,*)

angle = angle
CM = 0D0
do n=1,nmineral
 RT = ROTMAT(angle(1,n),angle(2,n),angle(3,n))
 CM = CM + xvol(n)*CT4(RT4(T4C(CI(:,:,n)),RT))
end do

call C66_OUT(CM)
!stop

open(iflout,file='../data/'//file_out(1:len_trim(file_out)))
write(iflout,*) (CM(i,i),i=1,3)
write(iflout,*) CM(2,3),CM(1,3),CM(1,2)
write(iflout,*) (CM(i,i),i=4,6)
do i=1,3
 write(iflout,*) CM(i,4),CM(i,5),CM(i,6)
end do
write(iflout,*) CM(5,6),CM(4,6),CM(4,5)
close(iflout)

write(*,*)
write(*,*) ' STABILITY '

KELM = KELC(CM)
call JACOBI(KELM,KVAL,VKEL,nrot)
call EIGSRT(KVAL,VKEL)
call V6_OUT(KVAL)
!
!
write(*,*) ' DECOMPOSITION '

call DEC_CONT(cm)
call DEC_ISO(cm)
call DEC_SCCA(cm)
write(*,*) ' ANISOTROPIC PERCENTAGE =',PERCA
write(*,*) ' HEXAGONAL PERCENTAGE =',PERCT
write(*,*)


write(*,*) ' END PROGRAM OK'
write(*,*)

contains

subroutine ELASTIC_INPUT(im)
integer, parameter :: iflin = 40
integer :: im,ifl
double precision, dimension(6,6) :: C
ifl = iflin + im
open(ifl,file=file_in(im))
read(ifl,*) (C(i,i),i=1,3)
read(ifl,*) C(2,3),C(1,3),C(1,2)
read(ifl,*) (C(i,i),i=4,6)
do i=1,3
 read(ifl,*) C(i,4),C(i,5),C(i,6)
end do
read(ifl,*) C(5,6),C(4,6),C(4,5)
CI(:,:,im) = SYM66(C)
close(ifl)
end subroutine ELASTIC_INPUT

end program ELASTIC_MM
