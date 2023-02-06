!
! d-rex matrix utility routines
!
! $Id: drex_matrix.f90,v 1.3 2011/04/12 06:18:38 becker Exp $
!
! based on Drex by Kaminski, and further routines by Jules 
! Browaeys as of July 2005
!
! modification by Thorsten Becker
!
! contains:
! - norms and low level routines
! - conversion for 6x6 matrices
! - conversion from 6x6 to 3,3,3,3 cijkl
! - conversion from 6x6 to [21]
! - rotation routines
!
!
!
!
! norms and low level routines
!
!
!
!****************************************************************
!
!  NDELTA   :   Kronecker function
function NDELTA(i,j)

  implicit none

  integer :: i,j
  integer :: NDELTA

  NDELTA=0
  if (i==j) NDELTA=1

end function NDELTA

!
! compute norm^2 of 4th order tensor
!
subroutine drex_tens4_norm2_ftrn(c4, norm2)
  implicit none
  integer :: i,j,k,l
  double precision, intent(in), dimension(3,3,3,3) :: C4
  double precision, intent(out) :: norm2

  norm2 = 0.0d0
  do i=1,3
     do j=1,3
        do k=1,3
           do l=1,3
              norm2 = norm2 + c4(i,j,k,l)*c4(i,j,k,l)
           end do
        end do
     end do
  end do
end subroutine drex_tens4_norm2_ftrn

!
!
! conversion between Voigt and kelvin 6x6 matrix format and 
! [21] vector
!
!
!
!****************************************************************
!
!  drex_d21v       : form 6x6 Voigt matrix from [21] 
!                    vector(independent)
!

subroutine drex_d21v_ftrn(x,c)

  implicit none
  ! input
  double precision, intent(in),dimension(21) :: X
  ! output
  double precision, intent(out),dimension(6,6) :: C
  double precision :: one_over_sqrt_two
  parameter (one_over_sqrt_two = 0.70710678118654752440084436210485d0)

  c = 0.0d0

  C(1,1) = x(1)
  C(2,2) = x(2)
  C(3,3) = x(3)
  c(2,3) = X(4)* one_over_sqrt_two 
  c(1,3) = X(5)* one_over_sqrt_two 
  c(1,2) = X(6)* one_over_sqrt_two 
  c(4,4) = x(7)/2d0
  c(5,5) = x(8)/2d0
  c(6,6) = x(9)/2d0
  c(1,4) = x(10)/2d0
  c(2,5) = x(11)/2d0
  c(3,6) = x(12)/2d0
  c(3,4) = x(13)/2d0
  c(1,5) = x(14)/2d0
  c(2,6) = x(15)/2d0
  c(2,4) = x(16)/2d0
  c(3,5) = x(17)/2d0
  c(1,6) = x(18)/2d0
  c(5,6) = x(19)/2d0* one_over_sqrt_two
  c(4,6) = x(20)/2d0* one_over_sqrt_two
  c(4,5) = x(21)/2d0* one_over_sqrt_two
  ! fill in other elements
  call drex_fullsym6_ftrn(c)

  return

end subroutine DREX_D21V_FTRN

!
!****************************************************************
!
!  drex_v21d       : form 21D vector from 6x6 Voigt 
!                    matrix (independent)
!
subroutine drex_v21d_ftrn(C,X)

  implicit none

  double precision, intent(in),dimension(6,6) :: C
  double precision, intent(out),dimension(21) :: X
  double precision :: sqrt_two
  parameter (sqrt_two = 1.41421356237309504d0)
  X  = 0D0
  X(1)  = C(1,1)
  X(2)  = C(2,2)
  X(3)  = C(3,3)
  X(4)  = sqrt_two*C(2,3)
  X(5)  = sqrt_two*C(1,3)
  X(6)  = sqrt_two*C(1,2)
  X(7)  = 2D0*C(4,4)
  X(8)  = 2D0*C(5,5)
  X(9)  = 2D0*C(6,6)
  X(10) = 2D0*C(1,4)
  X(11) = 2D0*C(2,5)
  X(12) = 2D0*C(3,6)
  X(13) = 2D0*C(3,4)
  X(14) = 2D0*C(1,5)
  X(15) = 2D0*C(2,6)
  X(16) = 2D0*C(2,4)
  X(17) = 2D0*C(3,5)
  X(18) = 2D0*C(1,6)
  X(19) = 2D0*sqrt_two*C(5,6)
  X(20) = 2D0*sqrt_two*C(4,6)
  X(21) = 2D0*sqrt_two*C(4,5)

  return

end subroutine DREX_V21D_FTRN

!
! Transform 6x6 Kelvin tensor to 21D vector
!
subroutine drex_v21kel(kel, x)
  double precision, intent(in),dimension(6,6) :: KEL
  double precision, intent(out),dimension(21) :: x
  double precision :: sqrt_two
  parameter (sqrt_two = 1.41421356237309504d0)
  
  X(1)  = KEL(1,1)
  X(2)  = KEL(2,2)
  X(3)  = KEL(3,3)
  X(4)  = sqrt_two*KEL(2,3)
  X(5)  = sqrt_two*KEL(1,3)
  X(6)  = sqrt_two*KEL(1,2)
  X(7)  = KEL(4,4)
  X(8)  = KEL(5,5)
  X(9)  = KEL(6,6)
  X(10) = sqrt_two*KEL(1,4)
  X(11) = sqrt_two*KEL(2,5)
  X(12) = sqrt_two*KEL(3,6)
  X(13) = sqrt_two*KEL(3,4)
  X(14) = sqrt_two*KEL(1,5)
  X(15) = sqrt_two*KEL(2,6)
  X(16) = sqrt_two*KEL(2,4)
  X(17) = sqrt_two*KEL(3,5)
  X(18) = sqrt_two*KEL(1,6)
  X(19) = sqrt_two*KEL(5,6)
  X(20) = sqrt_two*KEL(4,6)
  X(21) = sqrt_two*KEL(4,5)
end subroutine drex_v21kel


!
! low level operations on 6x6
!
!
!****************************************************************
!
!  drex_fullsym6 :   upper right 6x6 matrix into
!                    full symmetric 6x6 Voigt matrix
!
subroutine drex_fullsym6_ftrn(C)

  implicit none

  double precision, dimension(6,6) :: C
  C(2,1)=C(1,2)
  C(3,1)=C(1,3)
  C(4,1)=C(1,4)
  C(5,1)=C(1,5)
  C(6,1)=C(1,6)
  C(3,2)=C(2,3)
  C(4,2)=C(2,4)
  C(5,2)=C(2,5)
  C(6,2)=C(2,6)
  C(4,3)=C(3,4)
  C(5,3)=C(3,5)
  C(6,3)=C(3,6)
  C(5,4)=C(4,5)
  C(6,4)=C(4,6)
  C(6,5)=C(5,6)
  return

end subroutine DREX_FULLSYM6_FTRN
!
! Transform a 6x6 Voigt matrix to 6x6 Kelvin
!
function drex_kelc(C)
  double precision, intent(in), dimension(6,6) :: c
  double precision, dimension(6,6) :: drex_kelc
  double precision :: sqrt2
  parameter (sqrt2 = 1.41421356237309504d0)
  
  drex_kelc(1:3,1:3)=C(1:3,1:3)
  drex_kelc(1:3,4:6)=C(1:3,4:6)*sqrt2
  drex_kelc(4:6,1:3)=C(4:6,1:3)*sqrt2
  drex_kelc(4:6,4:6)=C(4:6,4:6)*2d0
end function DREX_KELC



!
! conversion between 6x6 voigt matrix and Cijkl [3,3,3,3]
! tensor
!

!****************************************************************
!
!  drex_TENS4_ftrn    :   transform 6x6 Voigt matrix
!                         into 4th order tensor
!
subroutine drex_tens4_ftrn(c,c4)

  implicit none

  integer :: i,j,k,l
  integer :: p,q
  integer :: NDELTA
  double precision, intent(in),dimension(6,6) :: C
  double precision, intent(out),dimension(3,3,3,3) :: C4

  C4=0D0

  do i=1,3
     do j=1,3
        do k=1,3
           do l=1,3

              p=NDELTA(i,j)*i+(1-NDELTA(i,j))*(9-i-j)
              q=NDELTA(k,l)*k+(1-NDELTA(k,l))*(9-k-l)
              C4(i,j,k,l)=C(p,q)

           end do
        end do
     end do
  end do
end subroutine DREX_TENS4_FTRN

!
!
! same if you have the ijkl index tensor
!
subroutine drex_tens4_ijkl_ftrn(c,ijkl,c4)

  implicit none

  integer :: i,j,k,l
  double precision, intent(in),dimension(6,6) :: C
  integer, intent(in),dimension(3,3) :: ijkl ! tensor of indices to form Cijkl from Sij
 
  double precision, intent(out),dimension(3,3,3,3) :: C4

  c4=0d0
  DO i = 1 , 3 
     DO j = 1 , 3 
        DO k = 1 , 3 
           DO l = 1 , 3
              c4(i,j,k,l) = c(ijkl(i,j),ijkl(k,l))
           END DO
        END DO
     END DO
  END DO

end subroutine DREX_TENS4_IJKL_FTRN
!
!****************************************************************
!
!  d42mat6     :   transform 4th order tensor
!                  into 6x6 Voigt matrix
!
subroutine drex_c42mat6(C4,C)

  implicit none

  integer :: i
  ! input
  double precision, intent(in), dimension(3,3,3,3) :: C4
  ! output
  double precision, intent(out), dimension(6,6) :: C

  C = 0D0

  do i=1,3
     C(i,i)=C4(i,i,i,i)
  end do
  do i=2,3
     C(1,i)=(C4(1,1,i,i)+C4(i,i,1,1))/2D0
     C(i,1)=C(1,i)
  end do
  C(2,3)=(C4(2,2,3,3)+C4(3,3,2,2))/2D0
  C(3,2)=C(2,3)

  do i=1,3
     C(i,4)=(C4(i,i,2,3)+C4(i,i,3,2)+                       &
          C4(2,3,i,i)+C4(3,2,i,i))/4D0
     C(4,i)=C(i,4)
  end do
  do i=1,3
     C(i,5)=(C4(i,i,1,3)+C4(i,i,3,1)+                       &
          C4(1,3,i,i)+C4(3,1,i,i))/4D0
     C(5,i)=C(i,5)
  end do
  do i=1,3
     C(i,6)=(C4(i,i,1,2)+C4(i,i,2,1)+                       &
          C4(1,2,i,i)+C4(2,1,i,i))/4D0
     C(6,i)=C(i,6)
  end do

  C(4,4)=(C4(2,3,2,3)+C4(2,3,3,2)+                          &
       C4(3,2,2,3)+C4(3,2,3,2))/4D0
  C(5,5)=(C4(1,3,1,3)+C4(1,3,3,1)+                          &
       C4(3,1,1,3)+C4(3,1,3,1))/4D0
  C(6,6)=(C4(2,1,2,1)+C4(2,1,1,2)+                          &
       C4(1,2,2,1)+C4(1,2,1,2))/4D0
  C(4,5)=(C4(2,3,1,3)+C4(2,3,3,1)+                          &
       C4(3,2,1,3)+C4(3,2,3,1)+                          &
       C4(1,3,2,3)+C4(1,3,3,2)+                          &
       C4(3,1,2,3)+C4(3,1,3,2))/8D0

  C(5,4)=C(4,5)
  C(4,6)=(C4(2,3,1,2)+C4(2,3,2,1)+                          &
       C4(3,2,1,2)+C4(3,2,2,1)+                          &
       C4(1,2,2,3)+C4(1,2,3,2)+                          &
       C4(2,1,2,3)+C4(2,1,3,2))/8D0
  C(6,4)=C(4,6)
  C(5,6)=(C4(1,3,1,2)+C4(1,3,2,1)+                          &
       C4(3,1,1,2)+C4(3,1,2,1)+                          &
       C4(1,2,1,3)+C4(1,2,3,1)+                          &
       C4(2,1,1,3)+C4(2,1,3,1))/8D0
  C(6,5)=C(5,6)

  return

end subroutine DREX_C42MAT6


!
!
!
!
!  matrix rotation
!
!
!
!
!
! get a rotation matrix for euler angles alpha, beta, gamma [rad]
! as in Dahlen and Tromp, p. 920
!
subroutine drex_calc_rotmat_cart(rot,alpha,beta,gamma)
  implicit none
  double precision, dimension(3,3) :: rot
  double precision :: sin_alpha,cos_alpha,alpha,beta,gamma
  double precision :: sin_beta,cos_beta
  double precision :: sin_gamma,cos_gamma

  sin_alpha=sin(alpha);cos_alpha=cos(alpha)
  sin_beta =sin(beta);cos_beta  =cos(beta)
  sin_gamma=sin(gamma);cos_gamma=cos(gamma)

  rot(1,1) = cos_alpha*cos_beta*cos_gamma - sin_alpha*sin_gamma; 
  rot(1,2) = sin_alpha*cos_beta*cos_gamma + cos_alpha*sin_gamma;
  rot(1,3) = -sin_beta*cos_gamma;
  rot(2,1) = -cos_alpha*cos_beta*sin_gamma -sin_alpha*cos_gamma;
  rot(2,2) = -sin_alpha*cos_beta*sin_gamma +cos_alpha*cos_gamma;
  rot(2,3) = sin_beta*sin_gamma;
  rot(3,1) = cos_alpha*sin_beta;
  rot(3,2) = sin_alpha*sin_beta;
  rot(3,3) = cos_beta;
end subroutine drex_calc_rotmat_cart
!
!
! rotate a sav(6,6) Voigt matrix by Euler angles alpha,beta,gamma [GIVEN IN RADIANS]
!
!
subroutine drex_rotate_6x6_rad_ftrn(sav,savr,alpha,beta,gamma)
  implicit none
  double precision, dimension(6,6) :: sav,savr
  double precision, dimension(3,3,3,3) :: cv,cvr
  double precision, dimension(3,3) :: rot
  double precision :: alpha,beta,gamma
  
  call drex_calc_rotmat_cart(rot,alpha,beta,gamma) !get the rotation matrix

  call drex_fullsym6_ftrn(sav)        ! fill in symmetric elements in [6,6] matrix
  call drex_tens4_ftrn(sav,cv)       !convert to 4th order tensor
  call drex_rot4(cv,rot,cvr)         !rotate 4th order tensor
  call drex_c42mat6(cvr,savr)        !convert 4th order tensor into 6,6 voigt

end subroutine drex_rotate_6x6_rad_ftrn
!
! same in degree
!
subroutine drex_rotate_6x6_deg_ftrn(sav,savr,alpha,beta,gamma)
  implicit none
  double precision, dimension(6,6) :: sav,savr
  double precision :: alpha,beta,gamma,ar,br,gr,pif
  parameter(pif=57.29577951308232d0)
  ar = alpha/pif
  br = beta/pif
  gr = gamma/pif

  call drex_rotate_6x6_rad_ftrn(sav,savr,ar,br,gr)
end subroutine drex_rotate_6x6_deg_ftrn

!
!
! rotate a sav(6,6) Voigt matrix by rotation matrix rot
!
!
subroutine drex_rotate_6x6_rot_ftrn(sav,rot,savr)
  implicit none
  double precision, intent(in), dimension(6,6) :: sav
  double precision, intent(out), dimension(6,6) :: savr
  double precision, intent(in), dimension(3,3) :: rot
  double precision, dimension(3,3,3,3) :: cv,cvr
  
  call drex_fullsym6_ftrn(sav)        ! fill in symmetric elements in [6,6] matrix
  call drex_tens4_ftrn(sav,cv)       !convert to 4th order tensor
  call drex_rot4(cv,rot,cvr)         !rotate 4th order tensor
  call drex_c42mat6(cvr,savr)        !convert 4th order tensor into 6,6 voigt

end subroutine drex_rotate_6x6_rot_ftrn
 

 
!
! rotate a c[21] vector using the rotation matrix rot
!
subroutine drex_rotate21(xin,rot,xout)
  implicit none
  double precision, dimension(21), intent(in) :: xin
  double precision, dimension(3,3), intent(in) :: rot
  double precision, dimension(21), intent(out) :: xout
  ! local
  double precision, dimension(6,6) :: x66
  double precision, dimension(3,3,3,3) :: xijin,xijout

  call drex_d21v_ftrn(xin,x66)  !form a [6,6] Voigt matrix
  call drex_tens4_ftrn(x66,xijin)  ! convert to 4th order tensor
  call drex_rot4(xijin,rot,xijout) !rotate 4th order
  call drex_c42mat6(xijout,x66)    !convert to [6,6] Voigt
  call drex_v21d_ftrn(x66,xout)    !form a [21] vector
end subroutine drex_rotate21
!
!****************************************************************
!
!  ROT4     :   rotation of 4th order tensor 
!  using a (3,3) rotation matrix
!
!
subroutine drex_rot4(C4,R,C4C)

  implicit none

  integer :: i1,i2,i3,i4,j1,j2,j3,j4
  double precision, intent(in),dimension(3,3,3,3) :: C4
  double precision, intent(in),dimension(3,3) :: R

  double precision, intent(out),dimension(3,3,3,3) :: C4C
  C4C = 0D0

  do i1=1,3
     do i2=1,3
        do i3=1,3
           do i4=1,3
              do j1=1,3
                 do j2=1,3
                    do j3=1,3
                       do j4=1,3
                          C4C(i1,i2,i3,i4) = C4C(i1,i2,i3,i4) +               &
                               R(i1,j1)*R(i2,j2)*&
                               R(i3,j3)*R(i4,j4)*C4(j1,j2,j3,j4)
                       end do
                    end do
                 end do
              end do
           end do
        end do
     end do
  end do


  return

end subroutine DREX_ROT4


