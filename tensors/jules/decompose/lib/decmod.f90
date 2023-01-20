!  decmod.f90
!  Elastic symmetry coordinates system and decomposition
!
!  DEC_CONT (subroutine) : contractions of elastic tensor
!  DEC_ISO  (subroutine) : isotropic parameters and percentage
!  DEC_SCCA (subroutine) : symmetry cartesian coordinates axes
!  DEC_PROJ (subroutine) : symmetry projectors
!
!  KISO   = isotropic incompressibility modulus
!  GISO   = isotropic shear modulus
!  PERCA  = tensor norm percentage anisotropic
!  PERCT  = tensor norm percentage transverse isotropic
!  DI     = tensor dilatationnal
!  VO     = tensor Voigt
!  VDI    = eigenvectors tensor dilatationnal
!  VVO    = eigenvectors tensor Voigt
!  SCC    = Symmetry Cartesian Coordinates vectors
!  X21    = 21-D elastic vector in SCC system
!
!  Module dependency : nrmod,tensmod

      module DECMOD

       use NRMOD
       use TENSMOD
       implicit none
       double precision :: KISO,GISO,PERCA,PERCT
       double precision, dimension(3,3) :: DI,VO,VDI,VVO,SCC
       double precision, dimension(21) :: X21

      contains

      subroutine DEC_CONT(C)
       !Calculates contractions of elastic tensor
       !i.e. dilatationnal and Voigt tensors
       !and produces eigenvectors and eigenvalues
       integer :: i,nrot
       double precision, dimension(3) :: dval,vval
       double precision, dimension(3,3) :: dd,vv
       double precision, dimension(6,6) :: C
       !Dilatationnal 3x3 tensor
       DI=0D0
       do i=1,3
        DI(1,1)=C(1,i)+DI(1,1)
        DI(2,2)=C(2,i)+DI(2,2)
        DI(3,3)=C(3,i)+DI(3,3)
        DI(2,1)=C(6,i)+DI(2,1)
        DI(3,1)=C(5,i)+DI(3,1)
        DI(3,2)=C(4,i)+DI(3,2)
       end do
       DI(1,2)=DI(2,1)
       DI(1,3)=DI(3,1)
       DI(2,3)=DI(3,2)
       !Voigt 3x3 tensor
       VO=0D0
       VO(1,1)=C(1,1)+C(6,6)+C(5,5)
       VO(2,2)=C(6,6)+C(2,2)+C(4,4)
       VO(3,3)=C(5,5)+C(4,4)+C(3,3)
       VO(2,1)=C(1,6)+C(2,6)+C(4,5)
       VO(1,2)=VO(2,1)
       VO(3,1)=C(1,5)+C(3,5)+C(4,6)
       VO(1,3)=VO(3,1)
       VO(3,2)=C(2,4)+C(3,4)+C(5,6)
       VO(2,3)=VO(3,2)
       !Diagonalization
       dd=DI
       vv=VO
       call JACOBI(dd,dval,VDI,nrot)
       call EIGSRT(dval,VDI)
       call JACOBI(vv,vval,VVO,nrot)
       call EIGSRT(vval,VVO)

      end subroutine DEC_CONT

      subroutine DEC_ISO(C)
       !Isotropic parameters and anisotropic percentage
       integer :: i
       double precision, dimension(21) :: X,XH,XD
       double precision, dimension(6,6) :: C
       !Isotropic parameters
       KISO=0D0
       GISO=0D0
       do i=1,3
        KISO=KISO+DI(i,i)
        GISO=GISO+VO(i,i)
       end do
       KISO=KISO/9D0
       GISO=0.1D0*GISO-0.3D0*KISO
       !Isotropic vector
       XH=0D0
       XH(1)=KISO+4D0*GISO/3D0
       XH(2)=KISO+4D0*GISO/3D0
       XH(3)=KISO+4D0*GISO/3D0
       XH(4)=sqrt(2D0)*(KISO-2D0*GISO/3D0)
       XH(5)=sqrt(2D0)*(KISO-2D0*GISO/3D0)
       XH(6)=sqrt(2D0)*(KISO-2D0*GISO/3D0)
       XH(7)=2D0*GISO
       XH(8)=2D0*GISO
       XH(9)=2D0*GISO
       !Anisotropic percentage
       X=V21KEL(KELC(C))
       XD=X-XH
       PERCA=sqrt(dot_product(XD,XD))/sqrt(dot_product(X,X))*1D2
      end subroutine DEC_ISO

      subroutine DEC_SCCA(C)
       !Calculates SCC directions as bissectrix
       !of each close pair of dilatationnal and
       !Voigt eigenvectors
       !Higher symmetry axis = 3-axis
       !t   = angles VO(1st index)/DI(2nd index)
       !bss = bissectrix direction intermediate
       !SCT = SCC test
       !dvm = array deviation for permutations
       double precision, parameter :: pi=3.141592654D0
       integer :: i,j
       integer, dimension(1) :: pos
       integer, dimension(3) :: npos
       double precision, dimension(3) :: bss,dvm
       double precision, dimension(21) :: XH,XD
       double precision, dimension(3,3) :: t,SCT
       double precision, dimension(6,6) :: C
       !VO/DI eigenvectors angle matrix
       do i=1,3
       do j=1,3
        t(j,i)=dot_product(VDI(:,i),VVO(:,j))
       end do
       end do
       where (abs(t)>=1D0) t=sign(1D0,t)
       t=acos(t)
       where (t>0.5D0*pi) t=t-pi
       !Indice position
       !Association VO eigenvector to DI
       npos=minloc(abs(t),dim=1)
       !Calculates bissectrix vector & SCC
       do i=1,3
        bss=VDI(:,i)+sign(1D0,t(npos(i),i))*VVO(:,npos(i))
        SCC(:,i)=bss/sqrt(dot_product(bss,bss))
       end do
       !Transpose : basis transfer
       SCC=transpose(SCC)
       !Permutation 1,2,3
       npos=(/1,2,3/)
       X21=V21KEL(KELC(CT4(RT4(T4C(C),SCC))))
       call DEC_PROJ(X21,XH,XD,dvm(1),5)
       !Permutations 2,3,1 & 3,1,2
       do j=2,3
        npos=cshift(npos,shift=1)
        do i=1,3
         SCT(i,:)=SCC(npos(i),:)
        end do
        X21=V21KEL(KELC(CT4(RT4(T4C(C),SCT))))
        call DEC_PROJ(X21,XH,XD,dvm(j),5)
       end do
       !Permutation for minimum deviation
       pos=minloc(dvm)
       !Transverse isotropic percentage
       print *,dvm(pos(1))
       PERCT=PERCA-dvm(pos(1))
       !SCC coordinates system
       SCC=cshift(SCC,shift=pos(1)-1,dim=1)
       !21-D elastic vector in SCC axes
       X21=V21KEL(KELC(CT4(RT4(T4C(C),SCC))))
      end subroutine DEC_SCCA

      subroutine DEC_PROJ(X,XH,XD,dev,NSYM)
       !Symmetry projectors
       !X    =  input 21-D vector
       !XH   =  projected vector
       !XD   =  deviation vector
       !dev  =  percentage deviation
       !NSYM = 13 monoclinic
       !     =  9 orthorhombic
       !     =  6 tetragonal
       !     =  5 hexagonal (transverse isotropic)
       !     =  2 isotropic
       integer :: NSYM
       double precision :: dev,sq2,isq2,i15
       double precision, dimension(21) :: X,XH,XD
       sq2=sqrt(2D0)
       isq2=1D0/sq2
       i15=1D0/15D0
       XH=0D0
       XD=0D0
       dev=0D0
       if (NSYM==13) then
        XH=X
        XH(10:11)=0D0
        XH(13:14)=0D0
        XH(16:17)=0D0
        XH(19:20)=0D0
       endif
       if (NSYM<=9) then
        XH(1:9)=X(1:9)
       endif
       if (NSYM==6) then
        XH(1)=0.5D0*(X(1)+X(2))
        XH(2)=XH(1)
        XH(4)=0.5D0*(X(4)+X(5))
        XH(5)=XH(4)
        XH(7)=0.5D0*(X(7)+X(8))
        XH(8)=XH(7)
       endif
       if (NSYM==5) then
        XH(1)=0.375D0*(X(1)+X(2))+0.25D0*X(6)*isq2+0.25D0*X(9)
        XH(2)=XH(1)
        XH(4)=0.5D0*(X(4)+X(5))
        XH(5)=XH(4)
        XH(6)=0.25D0*(X(1)+X(2))*isq2+0.75D0*X(6)-0.5D0*X(9)*isq2
        XH(7)=0.5D0*(X(7)+X(8))
        XH(8)=XH(7)
        XH(9)=0.25D0*(X(1)+X(2))-0.5D0*X(6)*isq2+0.5D0*X(9)
       endif
       if (NSYM==2) then
        XH(1)=3D0*(X(1)+X(2)+X(3))+sq2*(X(4)+X(5)+X(6))
        XH(1)=XH(1)+2D0*(X(7)+X(8)+X(9))
        XH(1)=XH(1)*i15
        XH(2)=XH(1)
        XH(3)=XH(1)
        XH(4)=sq2*(X(1)+X(2)+X(3))+4D0*(X(4)+X(5)+X(6))
        XH(4)=XH(4)-sq2*(X(7)+X(8)+X(9))
        XH(4)=XH(4)*i15
        XH(5)=XH(4)
        XH(6)=XH(4)
        XH(7)=2D0*(X(1)+X(2)+X(3))-sq2*(X(4)+X(5)+X(6))
        XH(7)=XH(7)+3D0*(X(7)+X(8)+X(9))
        XH(7)=XH(7)*i15
        XH(8)=XH(7)
        XH(9)=XH(7)
       endif
       XD=X-XH
       dev=sqrt(dot_product(XD,XD))/sqrt(dot_product(X,X))*1D2
      end subroutine DEC_PROJ

      end module DECMOD
