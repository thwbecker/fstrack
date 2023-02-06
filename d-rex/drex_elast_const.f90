!
! routines dealing with actual elastic constants
!
!
!
! constants as used in the original D-REX
!
!!! Stiffness matrix for Olivine (GigaPascals)
!
! this is olivine at room temperature, pressure, as in Abramson et al. (1997) 
! or Isaak (1992), roughly
!
subroutine  drex_get_olivine0_sav_ftrn(s0)
  implicit none
  double precision, dimension (6,6) :: s0
  s0 = 0.0d0
  S0(1,1) = 320.71d0 ; S0(1,2) = 69.84d0  ; S0(1,3) = 71.22d0
  S0(2,1) = S0(1,2)  ; S0(2,2) = 197.25d0 ; S0(2,3) = 74.8d0
  S0(3,1) = S0(1,3)  ; S0(3,2) = S0(2,3)  ; S0(3,3) = 234.32d0
  S0(4,4) = 63.77d0  ; S0(5,5) = 77.67d0  ; S0(6,6) = 78.36d0
end subroutine drex_get_olivine0_sav_ftrn
!
!!! Stiffness matrix for Enstatite (GPa)
!
! this is also at room temperature, pressure, from Chai et al. (1997)
!
! (the original D-REX had a typo, (1,2) was 79.6
subroutine  drex_get_enstatite0_sav_ftrn(s0)
  implicit none
  double precision, dimension (6,6) :: s0
  s0 = 0.0d0
  S0(1,1) = 236.9d0 ;      S0(1,2) = 79.9d0  ;      S0(1,3) = 63.2d0
  S0(2,1) = S0(1,2)  ;     S0(2,2) = 180.5d0 ;      S0(2,3) = 56.8d0
  S0(3,1) = S0(1,3)  ;     S0(3,2) = S0(2,3)  ;     S0(3,3) = 230.4d0
  S0(4,4) = 84.3d0  ;      S0(5,5) = 79.4d0  ;      S0(6,6) = 80.1d0
end subroutine drex_get_enstatite0_sav_ftrn
!
! values for olivine from Anderson & Isaak (1995) at 1500 K
! Forsterite 90 - Fayalite 10
!
subroutine  drex_get_olivine1_sav_ftrn(s0)
  implicit none
  double precision, dimension (6,6) :: s0
  s0 = 0.0d0
  S0(1,1) = 272d0 ;    S0(1,2) = 56.d0  ;   S0(1,3) = 60.d0
  S0(2,1) = S0(1,2)  ; S0(2,2) = 160d0 ;    S0(2,3) = 66d0
  S0(3,1) = S0(1,3)  ; S0(3,2) = S0(2,3)  ; S0(3,3) = 192d0
  S0(4,4) = 49d0  ;    S0(5,5) = 62d0  ;    S0(6,6) = 60d0
end subroutine drex_get_olivine1_sav_ftrn
!
! values for enstatite from Bass (1995) at room pressure, temperature
!
subroutine  drex_get_enstatite1_sav_ftrn(s0)
  implicit none
  double precision, dimension (6,6) :: s0
!!! Stiffness matrix for Olivine (GigaPascals)
  s0 = 0.0d0
  S0(1,1) = 225d0 ; S0(1,2) = 72.d0  ; S0(1,3) = 54.d0
  S0(2,1) = S0(1,2)  ; S0(2,2) = 178d0 ; S0(2,3) = 53d0
  S0(3,1) = S0(1,3)  ; S0(3,2) = S0(2,3)  ; S0(3,3) = 214d0
  S0(4,4) = 78d0  ; S0(5,5) = 76d0  ; S0(6,6) = 82d0
end subroutine drex_get_enstatite1_sav_ftrn
!
!
! varying functions
!

!
! values for olivine as function of temperature and pressure
! from Estey and Douglas (1986)
! give pressures in GPa
! give temperature in K, reference is 298
!
! this has been rotated to be consistent with the original KR
! and abramson as well as Isaak
!
! the indices 2 and 3, were flipped, kindof,
! this is equivalent to an Euler rotation alpha=90,beta=90,gamma=90
!
#define O_LOC_11 1,1
#define O_LOC_12 1,3
#define O_LOC_13 1,2
#define O_LOC_21 3,1
#define O_LOC_22 3,3
#define O_LOC_23 3,2
#define O_LOC_31 2,1
#define O_LOC_32 2,3
#define O_LOC_33 2,2
#define O_LOC_44 4,4
#define O_LOC_55 6,6
#define O_LOC_66 5,5
subroutine  drex_get_olivine_sav_ftrn(s0, T, p)
  implicit none
  double precision, dimension (6,6) :: s0,sp,st
  double precision :: p, T
!!! Stiffness matrix for Olivine (GigaPascals)
  s0 = 0.0d0 ; sp = 0.0d0 ; st = 0.0d0
  S0(O_LOC_11) = 323.7d0;      S0(O_LOC_12) =71.6d0  ;      S0(O_LOC_13) =66.4d0
  S0(O_LOC_21) = S0(O_LOC_12) ;S0(O_LOC_22) = 235.1d0  ;    S0(O_LOC_23) = 75.6d0
  S0(O_LOC_31) = S0(O_LOC_13) ;S0(O_LOC_32) = S0(O_LOC_23); S0(O_LOC_33) = 197.6d0
  S0(O_LOC_44) =64.62d0;       S0(O_LOC_55) =79.04d0;       S0(O_LOC_66) =78.05d0
! pressure derivatives
  SP(O_LOC_11) =7.98d0  ;      SP(O_LOC_12) =4.48d0  ;      SP(O_LOC_13) =4.74d0
  SP(O_LOC_21) = SP(O_LOC_12); SP(O_LOC_22) = 6.38d0  ;     SP(O_LOC_23) =3.76d0 
  SP(O_LOC_31) = SP(O_LOC_13); SP(O_LOC_32) = SP(O_LOC_23); SP(O_LOC_33) = 6.37d0 
  SP(O_LOC_44) =2.17d0 ;       SP(O_LOC_55) =2.31d0 ;       SP(O_LOC_66) =1.64d0 
! those are all negative, really
  ST(O_LOC_11) =0.0340d0;      ST(O_LOC_12) =.0094d0 ;      ST(O_LOC_13) =0.0105d0
  ST(O_LOC_21) = ST(O_LOC_12); ST(O_LOC_22) =0.0286d0  ;    ST(O_LOC_23) =0.0051d0
  ST(O_LOC_31) = ST(O_LOC_13); ST(O_LOC_32) = ST(O_LOC_23); ST(O_LOC_33) =0.0285d0
  ST(O_LOC_44) =0.0128d0;      ST(O_LOC_55) =0.0157d0;      ST(O_LOC_66) =0.013d0
! total tensor
  s0 = s0 + p * sp - (t-298d0) * st
end subroutine drex_get_olivine_sav_ftrn
#undef O_LOC_11 
#undef O_LOC_12 
#undef O_LOC_13 
#undef O_LOC_21 
#undef O_LOC_22 
#undef O_LOC_23 
#undef O_LOC_31 
#undef O_LOC_32 
#undef O_LOC_33 
#undef O_LOC_44 
#undef O_LOC_55 
#undef O_LOC_66 

!
! this is taking the pressure derivatives and M0 from Abramson et al (1997)
! and T derivatives from RAM smaple of Isaak (JGR, 1992)
!
subroutine  drex_get_olivine_new_sav_ftrn(s0, T, p)
  implicit none
  double precision, dimension (6,6) :: s0,sp,st
  double precision :: p, T
!!! Stiffness matrix for Olivine (GigaPascals)
  s0 = 0.0d0 ; sp = 0.0d0 ; st = 0.0d0
  S0(1,1) = 320.5d0; S0(1,2) =68.1d0 ;S0(1,3) =71.6d0
  S0(2,1) = S0(1,2) ;S0(2,2) =196.5d0;S0(2,3) =76.8d0
  S0(3,1) = S0(1,3) ;S0(3,2) =S0(2,3);S0(3,3) =233.5d0
  S0(4,4) = 64.0d0;  S0(5,5) =77.0d0; S0(6,6) =78.7d0
! pressure derivatives
  SP(1,1) = 6.54d0; SP(1,2) =3.86d0 ;SP(1,3) =3.57d0
  SP(2,1) = SP(1,2) ;SP(2,2) =5.38d0;SP(2,3) =3.37d0
  SP(3,1) = SP(1,3) ;SP(3,2) =SP(2,3);SP(3,3) =5.51d0
  SP(4,4) = 1.67d0;  SP(5,5) =1.81d0; SP(6,6) =1.93d0
!
! those are temperature derivatives, from RAM sample of Isaak 
!
  ST(1,1) =0.0402d0;ST(1,2) =0.0114d0 ;  ST(1,3) =0.0096d0
  ST(2,1) = ST(1,2); ST(2,2) =0.0310d0  ; ST(2,3) =0.0072d0
  ST(3,1) = ST(1,3); ST(3,2) =ST(2,3);    ST(3,3) =0.0353d0
  ST(4,4) =0.0126d0;ST(5,5) =0.0130d0;   ST(6,6) =0.0156d0
! total tensor
  s0 = s0 + p * sp - (t-298d0) * st
end subroutine drex_get_olivine_new_sav_ftrn

!
!
! values for enstatite as function of depth
! i flipped indices 1 and 3
! (this way, things are consistent with Chai et al)
!
! this is equivalent to a rotation with Euler angles of 0,90,0
!
#define E_LOC_11 3,3
#define E_LOC_12 3,2
#define E_LOC_13 3,1
#define E_LOC_21 2,3
#define E_LOC_22 2,2
#define E_LOC_23 2,1
#define E_LOC_31 1,3
#define E_LOC_32 1,2
#define E_LOC_33 1,1
#define E_LOC_44 6,6
#define E_LOC_55 5,5
#define E_LOC_66 4,4
!
! from Estey & Douglas
!
subroutine  drex_get_enstatite_sav_ftrn(s0, T, p)
  implicit none
  double precision, dimension (6,6) :: s0,sp,st
  double precision :: T, p 
!!! Stiffness matrix for orhtopyroxene, also from Estey and Douglas
  s0 = 0.0d0 ; sp = 0.0d0 ; st = 0.0d0
  S0(E_LOC_11) = 210.4d0;      S0(E_LOC_12) = 46.0d0  ;     S0(E_LOC_13) =54.8d0
  S0(E_LOC_21) = S0(E_LOC_12); S0(E_LOC_22) = 160.5d0 ;     S0(E_LOC_23) =71.00d0
  S0(E_LOC_31) = S0(E_LOC_13); S0(E_LOC_32) = S0(E_LOC_23); S0(E_LOC_33) = 228.6d0
  S0(E_LOC_44) = 77.66d0;      S0(E_LOC_55) = 75.48d0;      S0(E_LOC_66) =81.75d0
! pressure derivatives
  SP(E_LOC_11) = 16.42d0 ;     SP(E_LOC_12) =8.73d0  ;      SP(E_LOC_13) =9.09d0
  SP(E_LOC_21) = SP(E_LOC_12); SP(E_LOC_22) =9.190d0 ;      SP(E_LOC_23) =6.97d0 
  SP(E_LOC_31) = SP(E_LOC_13); SP(E_LOC_32) = SP(E_LOC_23); SP(E_LOC_33) =11.04d0 
  SP(E_LOC_44) = 2.75d0 ;      SP(E_LOC_55) =2.92d0 ;       SP(E_LOC_66) =2.38d0 
! those are all negative, really
  ST(E_LOC_11) =.05160d0;     ST(E_LOC_12) = 0.0107d0 ;     ST(E_LOC_13) =.03180d0
  ST(E_LOC_21) = ST(E_LOC_12);ST(E_LOC_22) = 0.0328d0 ;     ST(E_LOC_23) =0.0212d0
  ST(E_LOC_31) = ST(E_LOC_13);ST(E_LOC_32) = ST(E_LOC_23);  ST(E_LOC_33) =0.0352d0
  ST(E_LOC_44) = 0.0145d0;    ST(E_LOC_55) = 0.0138d0;      ST(E_LOC_66) =0.0131d0
! total tensor
  s0 = s0 + p * sp - (t-298d0) * st
  

end subroutine drex_get_enstatite_sav_ftrn



!
! En_0,80 Fs0,20 sample from Chai et al (JGR 1997)
!
subroutine  drex_get_enstatite_new_sav_ftrn(s0, T, p)
  implicit none
  double precision, dimension (6,6) :: s0,sp,st
  double precision :: T, p 
  !
!!! Stiffness matrix for orthopyroxene
  !
  s0 = 0.0d0 ; sp = 0.0d0 ; st = 0.0d0
  S0(1,1) = 231.0d0; S0(1,2) = 78.9d0  ;S0(1,3) =61.4d0
  S0(2,1) = S0(1,2); S0(2,2) = 169.8d0 ;S0(2,3) =49.1d0
  S0(3,1) = S0(1,3); S0(3,2) = S0(2,3); S0(3,3) =215.7d0
  S0(4,4) = 82.8d0; S0(5,5) = 76.5d0; S0(6,6) =78.1d0
  ! pressure derivatives
  SP(1,1) = 11.0d0; SP(1,2) =7.8d0  ; SP(1,3) =14.0d0
  SP(2,1) = SP(1,2); SP(2,2) =10.7d0 ; SP(2,3) =8.5d0 
  SP(3,1) = SP(1,3); SP(3,2) = SP(2,3); SP(3,3) =16.1d0
  SP(4,4) = 2.26d0 ; SP(5,5) =2.65d0 ;  SP(6,6) =2.76d0 
  ! those are all negative, really
  ST(E_LOC_11) =.05160d0;     ST(E_LOC_12) = 0.0107d0 ;     ST(E_LOC_13) =.03180d0
  ST(E_LOC_21) = ST(E_LOC_12);ST(E_LOC_22) = 0.0328d0 ;     ST(E_LOC_23) =0.0212d0
  ST(E_LOC_31) = ST(E_LOC_13);ST(E_LOC_32) = ST(E_LOC_23);  ST(E_LOC_33) =0.0352d0
  ST(E_LOC_44) = 0.0145d0;    ST(E_LOC_55) = 0.0138d0;      ST(E_LOC_66) =0.0131d0
  ! total tensor
  s0 = s0 + p * sp - (t-298d0) * st


end subroutine drex_get_enstatite_new_sav_ftrn



#undef E_LOC_11 
#undef E_LOC_12 
#undef E_LOC_13 
#undef E_LOC_21 
#undef E_LOC_22 
#undef E_LOC_23 
#undef E_LOC_31 
#undef E_LOC_32 
#undef E_LOC_33 
#undef E_LOC_44 
#undef E_LOC_55 
#undef E_LOC_66 
