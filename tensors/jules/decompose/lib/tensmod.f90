!  tensmod.f90
!  Independent tensor operations
!
!  DET33   (function)   : determinant matrix 3x3
!  SYM66   (function)   : symmetrization of 6x6 upper right triangle
!  T4C     (function)   : transform 6x6 matrix to 4th order tensor
!  NDELTA  (function)   : Kronecker function
!  CT4     (function)   : transform 4th order tensor to 6x6 matrix
!  RT4     (function)   : 4th order tensor rotation (RT4 = RRRR T4)
!  ROTMAT  (function)   : form 3D rotation matrix
!  KELC    (function)   : transform 6x6 matrix to 6x6 Kelvin
!  CKEL    (function)   : transform 6x6 Kelvin to 6x6 matrix
!  KELV21  (function)   : transform 21D vector to 6x6 Kelvin
!  V21KEL  (function)   : transform 6x6 Kelvin to 21D vector
!  V6_OUT  (subroutine) : screen output 6 vector
!  C66_OUT (subroutine) : screen output 6x6 matrix
!  C66_FILE(subroutine) : file output 6x6 matrix
!  C33_OUT (subroutine) : screen output 3x3 matrix
!  V21_OUT (subroutine) : screen output 21D vector

      module TENSMOD

      implicit none

      contains

      function DET33(MAT)
       !Determinant of 3x3 matrix
       double precision :: DET33
       double precision, dimension(3,3) :: MAT
       DET33=MAT(1,1)*(MAT(2,2)*MAT(3,3)-MAT(3,2)*MAT(2,3))
      end function DET33

      function SYM66(C)
       !Symmetrization of 6x6 upper right triangle
       integer :: i,j
       double precision, dimension(6) :: V
       double precision, dimension(6,6) :: C,SYM66
       logical, dimension(6,6) :: mk
       V=(/(C(i,i),i=1,6)/)
       mk=reshape((/(((i==j),i=1,6),j=1,6)/),(/6,6/))
       SYM66=C-unpack(V,mask=mk,field=0D0)+transpose(C)
      end function SYM66

      function T4C(C)
       !Transform symmetric 6x6 matrix to 4th order tensor
       integer :: i,j,k,l,p,q
       double precision, dimension(6,6) :: C
       double precision, dimension(3,3,3,3) :: T4C
       T4C=0D0
       do i=1,3
       do j=1,3
       do k=1,3
       do l=1,3
        p=NDELTA(i,j)*i+(1-NDELTA(i,j))*(9-i-j)
        q=NDELTA(k,l)*k+(1-NDELTA(k,l))*(9-k-l)
        T4C(i,j,k,l)=C(p,q)
       end do
       end do
       end do
       end do
      end function T4C

      function NDELTA(i,j)
       !Kronecker function
       integer :: i,j
       integer :: NDELTA
       NDELTA=0
       if (i==j) NDELTA=1
      end function NDELTA

      function CT4(T4)
       !Transform 4th order tensor to 6x6 matrix
       integer :: i
       double precision, dimension(6,6) :: CT4,C
       double precision, dimension(3,3,3,3) :: T4
       C=0D0
       !Quarter 11 diagonal
       do i=1,3
        C(i,i)=T4(i,i,i,i)
       end do
       !Quarter 11 upper right triangle
       C(2,3)=0.5D0*(T4(2,2,3,3)+T4(3,3,2,2))
       C(1,3)=0.5D0*(T4(1,1,3,3)+T4(3,3,1,1))
       C(1,2)=0.5D0*(T4(1,1,2,2)+T4(2,2,1,1))
       !Quarter 12
       do i=1,3
        C(i,4)=0.25D0* &
        (T4(i,i,2,3)+T4(i,i,3,2)+T4(2,3,i,i)+T4(3,2,i,i))
       end do
       do i=1,3
        C(i,5)=0.25D0* &
        (T4(i,i,1,3)+T4(i,i,3,1)+T4(1,3,i,i)+T4(3,1,i,i))
       end do
       do i=1,3
        C(i,6)=0.25D0* &
        (T4(i,i,1,2)+T4(i,i,2,1)+T4(1,2,i,i)+T4(2,1,i,i))
       end do
       !Quarter 22 diagonal
       C(4,4)=0.25D0* &
        (T4(2,3,2,3)+T4(2,3,3,2)+T4(3,2,2,3)+T4(3,2,3,2))
       C(5,5)=0.25D0* &
        (T4(1,3,1,3)+T4(1,3,3,1)+T4(3,1,1,3)+T4(3,1,3,1))
       C(6,6)=0.25D0* &
        (T4(2,1,2,1)+T4(2,1,1,2)+T4(1,2,2,1)+T4(1,2,1,2))
       !Quarter 22 upper right triangle
       C(4,5)=0.125D0* &
        (T4(2,3,1,3)+T4(2,3,3,1)+T4(3,2,1,3)+T4(3,2,3,1)+ &
         T4(1,3,2,3)+T4(1,3,3,2)+T4(3,1,2,3)+T4(3,1,3,2))
       C(4,6)=0.125D0* &
        (T4(2,3,1,2)+T4(2,3,2,1)+T4(3,2,1,2)+T4(3,2,2,1)+ &
         T4(1,2,2,3)+T4(1,2,3,2)+T4(2,1,2,3)+T4(2,1,3,2))
       C(5,6)=0.125D0* &
        (T4(1,3,1,2)+T4(1,3,2,1)+T4(3,1,1,2)+T4(3,1,2,1)+ &
         T4(1,2,1,3)+T4(1,2,3,1)+T4(2,1,1,3)+T4(2,1,3,1))
       !Transform symmetric
       CT4=SYM66(C)
      end function CT4

      function RT4(T4,R)
       !Rotation of 4th order tensor
       !R = rotation matrix & transfer matrix
       !RT4(ijkl)=R(im) R(jn) R(kp) R(lq) T4(mnpq)
       integer :: i1,i2,i3,i4,j1,j2,j3,j4
       double precision, dimension(3,3) :: R
       double precision, dimension(3,3,3,3) :: T4,RT4
       RT4=0D0
       !Rotated tensor index loop
       do i1=1,3
       do i2=1,3
       do i3=1,3
       do i4=1,3
         !Summation index loop
         do j1=1,3
         do j2=1,3
         do j3=1,3
         do j4=1,3
           RT4(i1,i2,i3,i4)=RT4(i1,i2,i3,i4)+ &
           R(i1,j1)*R(i2,j2)*R(i3,j3)*R(i4,j4)*T4(j1,j2,j3,j4)
         end do
         end do
         end do
         end do
       end do
       end do
       end do
       end do
      end function RT4

      function ROTMAT(phi1,theta,phi2)
       !Forms 3D rotation matrix
       !Euler angles in degree
       !Transfer matrix XYZ to E1,E2,E3
       !COLUMNS = E1,E2,E3 coordinates in XYZ
       double precision :: phi1,theta,phi2,pi180
       double precision, dimension(3,3) :: ROTMAT
       pi180=acos(-1D0)/180D0
       phi1=pi180*phi1
       theta=pi180*theta
       phi2=pi180*phi2
       !Euler angles rotation matrix
       !E1 vector in external coordinates
       ROTMAT(1,1)=cos(phi1)*cos(phi2)-sin(phi1)*sin(phi2)*cos(theta)
       ROTMAT(2,1)=sin(phi1)*cos(phi2)+cos(phi1)*sin(phi2)*cos(theta)
       ROTMAT(3,1)=sin(phi2)*sin(theta)
       !E2 vector in external coordinates
       ROTMAT(1,2)=-cos(phi1)*sin(phi2)-sin(phi1)*cos(phi2)*cos(theta)
       ROTMAT(2,2)=-sin(phi1)*sin(phi2)+cos(phi1)*cos(phi2)*cos(theta)
       ROTMAT(3,2)= cos(phi2)*sin(theta)
       !E3 vector in external coordinates
       ROTMAT(1,3)= sin(phi1)*sin(theta)
       ROTMAT(2,3)=-cos(phi1)*sin(theta)
       ROTMAT(3,3)= cos(theta)
      end function ROTMAT

      function KELC(C)
       !Transform 6x6 matrix to 6x6 Kelvin
       double precision, dimension(6,6) :: C,KELC
       KELC(1:3,1:3)=C(1:3,1:3)
       KELC(1:3,4:6)=C(1:3,4:6)*sqrt(2D0)
       KELC(4:6,1:3)=C(4:6,1:3)*sqrt(2D0)
       KELC(4:6,4:6)=C(4:6,4:6)*2D0
      end function KELC

      function CKEL(KEL)
       !Transform 6x6 Kelvin to 6x6 matrix
       double precision, dimension(6,6) :: KEL,CKEL
       CKEL(1:3,1:3)=KEL(1:3,1:3)
       CKEL(1:3,4:6)=KEL(1:3,4:6)/sqrt(2D0)
       CKEL(4:6,1:3)=KEL(4:6,1:3)/sqrt(2D0)
       CKEL(4:6,4:6)=KEL(4:6,4:6)/2D0
      end function CKEL

      function KELV21(XK)
       !Transform 21D vector to 6x6 Kelvin tensor
       integer :: i,j
       double precision, dimension(6) :: V33,V66
       double precision, dimension(9) :: V36
       double precision, dimension(21) :: X,XK
       double precision, dimension(6,6) :: KELV21,KEL
       logical, dimension(6,6) :: mk33,mk36,mk66
       KEL=0D0
       !Preprocess factors
       X=XK
       X(4:6)=XK(4:6)/sqrt(2D0)
       X(10:21)=XK(10:21)/sqrt(2D0)
       !Quarter 11 upper right triangle
       V33=(/X(1),X(6),X(2),X(5),X(4),X(3)/)
       mk33=reshape((/(((j>=i).and.(j<=3),i=1,6),j=1,6)/),(/6,6/))
       KEL=unpack(V33,mask=mk33,field=KEL)
       !Quarter 12
       V36=(/X(10),X(16),X(13),X(14),X(11),X(17),X(18),X(15),X(12)/)
       mk36=reshape((/(((j>3).and.(i<=3),i=1,6),j=1,6)/),(/6,6/))
       KEL=unpack(V36,mask=mk36,field=KEL)
       !Quarter 22 upper right triangle
       V66=(/X(7),X(21),X(8),X(20),X(19),X(9)/)
       mk66=reshape((/(((j>=i).and.(i>3),i=1,6),j=1,6)/),(/6,6/))
       KEL=unpack(V66,mask=mk66,field=KEL)
       !Transform symmetric
       KELV21=SYM66(KEL)
      end function KELV21 

      function V21KEL(KEL)
       !Transform 6x6 Kelvin tensor to 21D vector
       double precision, dimension(6,6) :: KEL
       double precision, dimension(21) :: V21KEL
          
       V21KEL(1)  = KEL(1,1)
       V21KEL(2)  = KEL(2,2)
       V21KEL(3)  = KEL(3,3)
       V21KEL(4)  = sqrt(2D0)*KEL(2,3)
       V21KEL(5)  = sqrt(2D0)*KEL(1,3)
       V21KEL(6)  = sqrt(2D0)*KEL(1,2)
       V21KEL(7)  = KEL(4,4)
       V21KEL(8)  = KEL(5,5)
       V21KEL(9)  = KEL(6,6)
       V21KEL(10) = sqrt(2D0)*KEL(1,4)
       V21KEL(11) = sqrt(2D0)*KEL(2,5)
       V21KEL(12) = sqrt(2D0)*KEL(3,6)
       V21KEL(13) = sqrt(2D0)*KEL(3,4)
       V21KEL(14) = sqrt(2D0)*KEL(1,5)
       V21KEL(15) = sqrt(2D0)*KEL(2,6)
       V21KEL(16) = sqrt(2D0)*KEL(2,4)
       V21KEL(17) = sqrt(2D0)*KEL(3,5)
       V21KEL(18) = sqrt(2D0)*KEL(1,6)
       V21KEL(19) = sqrt(2D0)*KEL(5,6)
       V21KEL(20) = sqrt(2D0)*KEL(4,6)
       V21KEL(21) = sqrt(2D0)*KEL(4,5)
      end function V21KEL

      subroutine V6_OUT(V)
       !Screen output 6 vector
       integer :: i
       double precision, dimension(6) :: V
       write(*,*)
       write(*, &
       '(1x,e12.4,1x,e12.4,1x,e12.4,1x,e12.4,1x,e12.4,1x,e12.4)') &
       (V(i),i=1,6)
       write(*,*)
      end subroutine V6_OUT

      subroutine C66_OUT(C)
       !Screen output 6x6 matrix
       integer :: i
       double precision, dimension(6,6) :: C
       write(*,*)
       do i=1,6
        write(*, &
         '(1x,e12.4,1x,e12.4,1x,e12.4,1x,e12.4,1x,e12.4,1x,e12.4)') &
         C(i,1),C(i,2),C(i,3),C(i,4),C(i,5),C(i,6)
       end do
       write(*,*)
      end subroutine C66_OUT

      subroutine C66_FILE(C,ifile)
       !File output 6x6 matrix
       integer :: i,ifile
       double precision, dimension(6,6) :: C
       write(ifile,*)
       do i=1,6
       write(ifile, &
         '(1x,e12.4,1x,e12.4,1x,e12.4,1x,e12.4,1x,e12.4,1x,e12.4)') &
         C(i,1),C(i,2),C(i,3),C(i,4),C(i,5),C(i,6)
       end do
       write(ifile,*)
      end subroutine C66_FILE

      subroutine C33_OUT(C)
       !Screen output 3x3 matrix
       integer :: i
       double precision, dimension(3,3) :: C
       write(*,*)
       do i=1,3
        write(*, &
         '(1x,e12.4,1x,e12.4,1x,e12.4)') C(i,1),C(i,2),C(i,3)
       end do
       write(*,*)
      end subroutine C33_OUT

      subroutine V21_OUT(X)
       !Screen output 21D vector
       double precision, dimension(21) :: X
       write(*,*)
       write(*,'(1x,e12.4,1x,e12.4,1x,e12.4)') X(1),X(2),X(3)
       write(*,'(1x,e12.4,1x,e12.4,1x,e12.4)') X(4),X(5),X(6)
       write(*,'(1x,e12.4,1x,e12.4,1x,e12.4)') X(7),X(8),X(9)
       write(*,'(1x,e12.4,1x,e12.4,1x,e12.4)') X(10),X(11),X(12)
       write(*,'(1x,e12.4,1x,e12.4,1x,e12.4)') X(13),X(14),X(15)
       write(*,'(1x,e12.4,1x,e12.4,1x,e12.4)') X(16),X(17),X(18)
       write(*,'(1x,e12.4,1x,e12.4,1x,e12.4)') X(19),X(20),X(21)
       write(*,*)
      end subroutine V21_OUT

      end module TENSMOD
