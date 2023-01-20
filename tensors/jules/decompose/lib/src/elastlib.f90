! elastlib.f90
! Library of subroutines for elastic tensor conversion
!
! ELAST_INPUT (subroutine)     : elastic tensor input
! ELAST_MODEL (subroutine)     : calculation elastic model
! ELAST_SYM   (int.subroutine) : symmetry decomposition
! ELAST_LPO   (subroutine)     : elastic tensor LPO conversion
! ELAST_VOIGT (int.subroutine) : Voigt conversion
! ELAST_REUSS (int.subroutine) : Reuss conversion
! ELAST_LOG   (int.subroutine) : log conversion
! ELAST_MIN   (int.subroutine) : several minerals type conversion

subroutine ELAST_INPUT
use TENSMOD
use DREXMOD
implicit none
double precision :: d1
double precision, dimension(6,6) :: C1
!Mineral 1 : olivine
open(36,file='../../../files/elastic/Fo90Fa10.1500K')
 read(36,*) d1
 read(36,*) C1(1,1),C1(2,2),C1(3,3)
 read(36,*) C1(2,3),C1(1,3),C1(1,2)
 read(36,*) C1(4,4),C1(5,5),C1(6,6)
 read(36,*) C1(1,4),C1(1,5),C1(1,6)
 read(36,*) C1(2,4),C1(2,5),C1(2,6)
 read(36,*) C1(3,4),C1(3,5),C1(3,6)
 read(36,*) C1(5,6),C1(4,6),C1(4,5)
close(36)
OL%C=SYM66(C1)
!Mineral 2 : enstatite
open(36,file='../../../files/elastic/En.PTroom')
 read(36,*) d1
 read(36,*) C1(1,1),C1(2,2),C1(3,3)
 read(36,*) C1(2,3),C1(1,3),C1(1,2)
 read(36,*) C1(4,4),C1(5,5),C1(6,6)
 read(36,*) C1(1,4),C1(1,5),C1(1,6)
 read(36,*) C1(2,4),C1(2,5),C1(2,6)
 read(36,*) C1(3,4),C1(3,5),C1(3,6)
 read(36,*) C1(5,6),C1(4,6),C1(4,5)
close(36)
EN%C=SYM66(C1)
end subroutine ELAST_INPUT

subroutine ELAST_MODEL(ictype,isym,E21,pp)
use NRMOD
use TENSMOD
use DECMOD
implicit none
integer :: ictype,isym
double precision :: pp
double precision, dimension(21) :: E21
double precision, dimension(6,6) :: C
pp=0D0
call ELAST_LPO(ictype,C)
select case (isym)
case (0)
 E21=V21KEL(KELC(C))
case default
 call ELAST_SYM
end select

contains

subroutine ELAST_SYM
double precision :: dv
double precision, dimension(21) :: XH,XD
double precision, dimension(6,6) :: CE
call DEC_CONT(C)
call DEC_ISO(C)
call DEC_SCCA(C)
!Hexagonal proportion
!1 : hexagonal anisotropy
!0 : no hexagonal anisotropy
pp=PERCT/PERCA
call DEC_PROJ(X21,XH,XD,dv,isym)
CE=CKEL(KELV21(XH))
C=CT4(RT4(T4C(CE),transpose(SCC)))
E21=V21KEL(KELC(C))
end subroutine ELAST_SYM

end subroutine ELAST_MODEL

subroutine ELAST_LPO(ictype,C)
use NRMOD
use TENSMOD
use DREXMOD
implicit none
integer :: ictype,i,j,n
double precision, dimension(3,3) :: ROT
double precision, dimension(6,6) :: C,C1,C2
type(MINERAL), pointer :: PMX
logical, dimension(6,6) :: mk
C=0D0
C1=0D0
C2=0D0
PMX => OL
select case (ictype)
case(1)
 call ELAST_VOIGT(C1)
 if (nmx==2) then
  PMX => EN
  call ELAST_VOIGT(C2)
  call ELAST_MIN(1)
 else
  C=C1
 endif
case(2)
 call ELAST_REUSS(C1)
 if (nmx==2) then
  PMX => EN
  call ELAST_REUSS(C2)
  call ELAST_MIN(2)
 else 
  C=C1
 endif
case(3)
 mk=reshape((/(((i==j),i=1,6),j=1,6)/),(/6,6/))
 call ELAST_LOG(C1)
 if (nmx==2) then
  PMX => EN
  call ELAST_LOG(C2)
  call ELAST_MIN(3)
 else
  C=C1
 endif
end select

contains

subroutine ELAST_VOIGT(CE)
double precision, dimension(6,6) :: CE,CI
CI=PMX%C
CE=0D0
do n=1,NG
 ROT=transpose(PMX%F%acs(:,:,n))
 CE=CE+PMX%F%odf(n)*CT4(RT4(T4C(CI),ROT))
end do
end subroutine ELAST_VOIGT

subroutine ELAST_REUSS(CE)
double precision, dimension(6,6) :: CE,CI
call LUINV(KELC(PMX%C),6,CE)
CI=CKEL(CE)
CE=0D0
do n=1,NG
 ROT=transpose(PMX%F%acs(:,:,n))
 CE=CE+PMX%F%odf(n)*CT4(RT4(T4C(CI),ROT))
end do
call LUINV(KELC(CE),6,CI)
CE=CKEL(CI)
end subroutine ELAST_REUSS

subroutine ELAST_LOG(CE)
integer :: nrot
double precision, dimension(6) :: val
double precision, dimension(6,6) :: CE,CI,vv
mk=reshape((/(((i==j),i=1,6),j=1,6)/),(/6,6/))
CI=KELC(PMX%C)
call JACOBI(CI,val,vv,nrot)
CE=unpack(log(val),mask=mk,field=0D0)
CI=CT4(RT4(T4C(CKEL(CE)),vv))
CE=0D0
do n=1,NG
 ROT=transpose(PMX%F%acs(:,:,n))
 CE=CE+PMX%F%odf(n)*CT4(RT4(T4C(CI),ROT))
end do
CI=KELC(CE)
call JACOBI(CI,val,vv,nrot)
CI=unpack(exp(val),mask=mk,field=0D0)
CE=CT4(RT4(T4C(CKEL(CI)),vv))
end subroutine ELAST_LOG

subroutine ELAST_MIN(icm)
integer :: icm,nrot
double precision, dimension(6) :: val
double precision, dimension(6,6) :: CI,vv
select case (icm)
case(1)
 C=OL%x*C1+EN%x*C2
case(2)
 call LUINV(KELC(C1),6,CI)
 C1=CKEL(CI)
 call LUINV(KELC(C2),6,CI)
 C2=CKEL(CI)
 C=OL%x*C1+EN%x*C2
 call LUINV(KELC(C),6,CI)
 C=CKEL(CI)
case(3)
 CI=KELC(C1)
 call JACOBI(CI,val,vv,nrot)
 CI=unpack(log(val),mask=mk,field=0D0)
 C1=CT4(RT4(T4C(CKEL(CI)),vv))
 CI=KELC(C2)
 call JACOBI(CI,val,vv,nrot)
 CI=unpack(log(val),mask=mk,field=0D0)
 C2=CT4(RT4(T4C(CKEL(CI)),vv))
 C=OL%x*C1+EN%x*C2
 CI=KELC(C)
 call JACOBI(CI,val,vv,nrot)
 CI=unpack(exp(val),mask=mk,field=0D0)
 C=CT4(RT4(T4C(CKEL(CI)),vv))
end select
end subroutine ELAST_MIN

end subroutine ELAST_LPO
