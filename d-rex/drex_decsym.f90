!
!  subroutines from D-REX by E. Kaminski as of 01/11/2003
!
!  added code from Jules Browaeys as of 07/2005
!
!  minor changes by TWB
!
!  bug fix by sebastien chevrot 04/2011
!
!  bug fix of old version of symmetry axes by Manabu, 07/2011
!
!  $Id: drex_decsym.f90,v 1.19 2011/04/12 06:18:38 becker Exp becker $
!
! Equivalent transverse isotropic media calculated using the theory of   
! Jules Browaeys. 
!
!

!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  subroutine for elastic decomposition
!  into transverse isotropy
!
!  drex_decsym   :   decomposition into transverse
!                     isotropy symmetry
!
!  Library of independent subroutines
!  for elasticity and tensor
!

!
!****************************************************************
!
!  DECSYM   :   decomposition into transverse
!               isotropy symmetry
subroutine drex_decsym_ftrn(c, kmod, gmod, symm_frac, tiaxis, &
     dc_iso, dc_hex, dc_tet, dc_ort, dc_mon, dc_tri, vel, &
     dc_orient, scc_rot,scca_old_mode)
  !
  ! input:  ced(6,6) elastic matrix in Voigt format
  !
  ! for finding best axes:
  ! scca_old_mode: 1: only find best-fit, and align
  !                0: find best and worst fit, and use those as a coordinate system
  !
  ! output: 
  !
  ! kmod, gmod: K and G moduli of isotropic component
  ! symm_frac[6]
  !
  ! fractional norms for:
  ! 1: aniisotropy total
  ! 2: anisotropy due to hex
  ! 3: anisotropy due to tet
  ! 4: anisotropy due to ort
  ! 5: anisotropy due to mon
  ! 6: anisotropy due to tri
  !
  ! tiaxis(6): normalized symmetry axis vector for the best fit hexagonal (TI) axis
  !            the first three are the best-fit, the last three the 
  !            worst fit (for slow symmetries)
  !
  ! dc_iso, dc_hex, dc_tet, dc_ort, dc_mon, dc_tri:
  !
  ! tensor components as (6,6) Voigt matrices, dc_orient is the c tensor in the hexagonal SCC system
  ! for isotropic, hexagonal, tetragonal, orthorhombic, monoclinic, and triclinic parts
  !
  ! vel(7) velocitiy factors
  !
  ! vel(1) C_11 for ani in hex sytem for p velocities
  ! vel(2) C_33
  !
  ! vel(3) = c_{33}^0,pseudo v_p^2\rho
  ! vel(4) = c_{44}^0 pseudo v_S^2\rho
  !
  !
  ! vel(5) = epsilon*(c_33^0 or c_11^0) = (C_11-C-33)/2, for P wave anisotropy
  ! vel(6) = gamma * c_44^0             = (C_66-C_44)/2, for S wave anisotropy
  ! vel(7) = delta * (c_33^0 or c_11^0)   = (C_13 - C_33+2 C_44), for ellipticity
  !         vel(5-7) have to be normalized, depending on choices from 
  !         Thomsen 1988 (c_33 for epsilon and gamma) or Mensch T and Rasolofosaon P (1997) (c_11)
  !
  ! vel(8), vel(9): vs1*sqrt(rho), vs2*sqrt(rho)
  !
 
  implicit none
  
  ! input
  double precision, intent(in),dimension(6,6) :: c
  integer, intent(in) :: scca_old_mode

  ! output
  double precision, intent(out) :: kmod, gmod ! isotropic tensors
  double precision, intent(out),dimension(9) :: vel 

  double precision, intent(out),dimension(6) :: symm_frac 

  double precision, intent(out),dimension(6) :: tiaxis ! VTI axis, best, worst
  double precision, intent(out),dimension(6,6) :: dc_iso,dc_hex,dc_tet,dc_ort,dc_mon,&
       dc_tri,dc_orient
  double precision, intent(out), dimension (3,3) :: scc_rot ! rotation matrix
  !
  ! local
  double precision, dimension(21) :: x,xp,xd,xiso,xhex,xani
  double precision :: xn,dev,devo
!  double precision :: dev5,dev2,dev13,dev9,dev6
  double precision, dimension(3,3) :: di, vo,vvo,vdi
  double precision, dimension(3) :: dval, vval
  ! init
  symm_frac = 0d0
  !
  !
  call drex_fullsym6_ftrn(c)! fill in symmetric elements in [6,6] tensor
  

!  print *,c
!  stop

  !
  ! get the contraction of the Voigt matrix
  !
  call drex_dec_cont(c,di,dval,vdi,vo,vval,vvo)
  !
  ! convert the Voigt matrix to [21] format
  call drex_v21d_ftrn(c,x)
  !
  ! norm of total tensor
  !
  xn = sqrt(dot_product(x,x))
  !
  ! get isotropic component xiso[21] and anisotropic vector xd[21]
  ! symm frac will be the anisotropic norm (of xd[])
  ! i.e. the difference between the full tensor and 
  ! the isotropic one 
  call drex_dec_iso(x,xiso,xd,symm_frac(1),kmod,gmod,di,vo)
  !
  call drex_d21v_ftrn(xiso,dc_iso)
  !
  ! rotate into Scc system, this takes the voigt matrix
  !
  ! 
  ! old_mode: 1: only find best-fit, and align
  !           0: find best and worst fit, and use those as a coordinate system
  call drex_dec_scca(c,vdi,vvo,scc_rot,x,scca_old_mode) ! new or old version?
  

  
  ! original tensor in rotated frame
  call drex_d21v_ftrn(x,dc_orient)
  !
  ! fast VTI axis
  !
  tiaxis(1:3)=scc_rot(3,:)          !best fit 
  tiaxis(4:6)=scc_rot(1,:)          !worst fit
  ! normalize
  tiaxis(1:3)=tiaxis(1:3)/sqrt(sum(tiaxis(1:3)*tiaxis(1:3)))
  tiaxis(4:6)=tiaxis(4:6)/sqrt(sum(tiaxis(4:6)*tiaxis(4:6)))
  !
  ! remove isotropic (which does not depend on this rotation)
  xani = x - xiso               ! anisotropic tensor in hex system
  !
  ! hexagonal projection
  !
  ! this will compute the xd difference between the 
  ! anisotropic tensor, in the scca system, and the 
  ! hexagonal part, xp, and return the norm of the difference xd
  !
  call drex_dec_proj(xani,xp,xd,dev,5)
  ! this is the hexagonala component, i.e. the difference of 
  ! anisotropic total norm and the non-hexagonal part

  symm_frac(2) = symm_frac(1) - dev


  call drex_d21v_ftrn(xp,dc_hex)
  xhex = xiso + xp ! best fit hexagonal = isotropic + hex component
  !
  ! use best-fit hexagonal to compute eps, gamma, delta factors
  !
  !
  ! fast and slow P velocities from hexagonal system tensor
  !
  ! normally, vp2 > vp1, but could be reversed if the symmetry axis 
  ! is slow
  !
  vel(1) = sqrt(xhex(3))                  ! vp1*sqrt(dens)
  vel(2) = sqrt(xhex(1))                  ! vp2*sqrt(dens)
  !
  ! shear wave
  !
  vel(8)=sqrt(xhex(7)/2.0d0)           ! vs1 = sqrt(c_hex(4,4)/dens)
  vel(9)=sqrt(xhex(9)/2.0d0)           ! vs2 = sqrt(c_hex(6,6)/dens)
  ! 
  ! reference values, based on best-fit hex
  !
  ! mean of c11 and c11, c33ref, P velocity
  ! (better to use C^iso_11)
  !
  vel(3) = (xhex(1)+xhex(3))/2.0d0 ! c_33^0 = (c_11+c_33)/2; c_11=x(1);   c_33=x(3)
  ! mean of c44 and c66, S velocity
  ! (better to use C^iso_44
  !
  vel(4) = (xhex(7)+xhex(9))/4.0d0 ! c_44^0 = (c_44+c_66)/2; c_44=x(7)/2; c_66=x(9)/2
 
  !
  ! epsilon, gamma, delta
  ! use C_XX = C_33 for Thomsen (1988) definition and
  !            C_11 for Mensch T and Rasolofosaon P (1997) 
  !
  ! epsilon = (C_11-C_33)/2C_XX^0; vel(3) is eps*c_XX^0
  ! gamma   = (C_66-C_44)/2C_44^0; vel(4) is gamma*c_44^0
  ! delta   = (C_13-C_33+2 C_44)/(C_XX^0); C_13 = x(5)/sqrt(2); vel(5)=delta*c_XX^0
  !
  vel(5) = (xhex(1)-xhex(3))/2.0d0 ! non-normalized eps
  vel(6) = (xhex(9)-xhex(7))/4.0d0 ! non-normalize gamma
  vel(7) = (xhex(5)/sqrt(2.0d0)-xhex(3)+xhex(7)) ! delta
  !
  ! tetragonal fraction
  devo = dev
  call drex_dec_proj(xd,xp,xd,dev,6)
  symm_frac(3) = devo - dev
  call drex_d21v_ftrn(xp,dc_tet)
  ! orthorhombic
  devo=dev
  call drex_dec_proj(xd,xp,xd,dev,9)
  symm_frac(4) = devo-dev
  call drex_d21v_ftrn(xp,dc_ort)
  ! monoclinic
  devo=dev
  call drex_dec_proj(xd,xp,xd,dev,13)
  symm_frac(5) = devo-dev
  call drex_d21v_ftrn(xp,dc_mon)
  ! triclinic
  symm_frac(6) = dev
  call drex_d21v_ftrn(xd,dc_tri)

  !
  ! normalize the norms, these are now symmetry fractions
  symm_frac = symm_frac/xn
  
  
  !call drex_dec_proj(x,xp,xd,dev13,13) !monoclinic
  !call drex_dec_proj(x,xp,xd,dev9,9) !orth
  !call drex_dec_proj(x,xp,xd,dev6,6) !tetra
  !call drex_dec_proj(x,xp,xd,dev5,5) !hex
  !call drex_dec_proj(x,xp,xd,dev2,2) !iso
  !write(*,*) ' monoclinic % =',(dev9-dev13)/xn
  !write(*,*) ' orthorhombic % =',(dev6-dev9)/xn
  !write(*,*) ' tetragonal % =',(dev5-dev6)/xn
  !write(*,*) ' hexagonal % =',(dev2-dev5)/xn
  !write(*,*) ' anisotropic % =',dev2/xn

end  subroutine drex_decsym_ftrn




!
!****************************************************************
!
!  EIGSRT   :   order eigenvectors with
!               increasing eigenvalues

subroutine drex_eigsrt(d,v,n,np)

  ! Order eigenvalues and eigenvectors
  ! 1 : max
  ! 2 : mid
  ! 3 : min

  implicit none

  integer :: np,n
  integer :: i,j,k
  double precision, dimension(np) :: d
  double precision, dimension(np,np) :: v
  double precision :: p

  do i=1,n-1
     k=i
     p=d(i)
     do j=i+1,n
        if (d(j)>=p) then
           k=j
           p=d(j)
        end if
     end do
     if (k/=i) then
        d(k)=d(i)
        d(i)=p
        do j=1,n
           p=v(j,i)
           v(j,i)=v(j,k)
           v(j,k)=p
        end do
     end if
  end do

  return

end subroutine DREX_EIGSRT

!
!****************************************************************
!
!  PERMUT   :   permutation of index 123

subroutine permut(INDEX,PERM)

  implicit none

  integer,intent(in) :: INDEX
  integer, intent(out), dimension(3) :: PERM

  if (INDEX==1) then
     PERM(1)=1
     PERM(2)=2
     PERM(3)=3
  end if
  if (INDEX==2) then
     PERM(1)=2
     PERM(2)=3
     PERM(3)=1
  end if
  if (INDEX==3) then 
     PERM(1)=3
     PERM(2)=1
     PERM(3)=2
  endif

  return

end subroutine PERMUT



subroutine drex_dec_cont(c,di,dval,vdi,vo,vval,vvo)
  !
  ! Calculates contractions of elastic cij voigt matrix
  ! i.e. dilatationnal and Voigt tensors
  ! and produces eigenvectors and eigenvalues
  !
  use drex_nrmod
  implicit none
  integer :: i,nrot
  double precision, intent(out),dimension(3) :: dval,vval
  double precision, intent(out),dimension(3,3) :: di,vo,vdi,vvo
  double precision, intent(in), dimension(6,6) :: C
  double precision, dimension(3,3) :: dicopy,vocopy
  !
  ! Dilatationnal 3x3 tensor
  !
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
  ! Diagonalization
  ! make copies, jacobi is destructive
  dicopy=di
  vocopy=vo
  ! get eigensystem
  call drex_nr_jacobi(dicopy,dval,VDI,nrot)
  call drex_nr_eigsrt(dval,VDI)
  call drex_nr_jacobi(vocopy,vval,VVO,nrot)
  call drex_nr_eigsrt(vval,VVO)

end subroutine DREX_DEC_CONT

!
! compute isotropic projection and elasticity constants
!
subroutine drex_dec_iso(x,xp,xd,dev,kiso,giso,di,vo)
  implicit none
  !Isotropic parameters and anisotropic percentage
  integer :: i
  double precision, intent(in),dimension(21) :: x
  double precision, intent(in), dimension(3,3) :: di,vo
  double precision, intent(out) :: kiso,giso,dev
  double precision, intent(out),dimension(21) :: xp,xd
  double precision :: vsum, dsum
  !
  ! isotropic parameters
  !
  vsum=0d0
  dsum=0d0
  do i=1,3
     dsum=dsum+di(i,i)
     vsum=vsum+vo(i,i)
  end do
  kiso=dsum/9d0
  giso=0.1d0*vsum - dsum/30.0d0

  !Isotropic projection vector
  xp=0d0
  xp(1)=KISO+4D0*GISO/3D0
  XP(2)=xp(1)
  XP(3)=xp(1)
  XP(4)=sqrt(2D0)*(KISO-2D0*GISO/3D0)
  XP(5)=xp(4)
  XP(6)=xp(4)
  XP(7)=2D0*GISO
  XP(8)=xp(7)
  XP(9)=xp(7)
  ! difference 
  xd = x - xp
  dev = sqrt(dot_product(xd,xd))

  
end subroutine DREX_DEC_ISO
!
! get scc system, rotate, and obtain hex fraction
!
! this version is adapted from the DREX v.1 and v.2 distributions
! (not tested well)
!
subroutine drex_dec_scca_drex(c,vdi,vvo,scc,x21r)
  implicit none
  !
  ! Calculates SCC directions as bissectrix
  ! of each close pair of dilatationnal and
  ! Voigt eigenvectors
  ! Higher symmetry axis = 3-axis
  !
  !t   = angles VO(1st index)/DI(2nd index)
  !SCT = SCC test
  !dvm = array deviation for permutations
  !
  double precision, parameter :: pi=3.14159265358979323846d0
  !
  ! input voigt matrix
  !
  double precision, intent(in),dimension(6,6) :: c
  !
  ! input eigensystem of contractions
  !
  double precision, intent(in),dimension(3,3) :: vdi,vvo
  ! output
  double precision, intent(out),dimension(3,3) :: scc
  double precision, intent(out),dimension(21) :: x21r
  ! local
  integer :: i1,i2,ndvc
  integer, dimension(3) :: ihs
  double precision :: dev,sdv,adv,advc,scn
  double precision, dimension(3,3) :: vecdi,vecvo
  double precision, dimension(6,6) :: cr
  double precision, dimension(21) :: xp,xd
  !
  !VO/DI eigenvectors angle matrix
  !
  vecdi = vdi                   ! copies
  vecvo = vvo
  
  print *,'drex version of scca is not tested well'

  do i1=1,3
     ndvc=0
     advc=10D0
     scn=0D0
     do i2=1,3
        sdv=dot_product(vecdi(:,i1),vecvo(:,i2))
        if (abs(sdv)>=1D0) sdv=sign(1D0,sdv)
        adv=acos(sdv)
        if (sdv<0d0) adv=acos(-1d0)-adv
        if (adv<advc) then
           ndvc=sign(1d0,sdv)*i2
           advc=adv
        endif
     end do
     
     vecdi(:,i1)=(vecdi(:,i1)+ndvc*vecvo(:,abs(ndvc)))/2d0
     scn=sqrt(vecdi(1,i1)**2d0+vecdi(2,i1)**2d0+            &
          vecdi(3,i1)**2d0)
     vecdi(:,i1)=vecdi(:,i1)/scn
  end do
  ! Higher symmetry axis
  scc = transpose(vecdi)

  call drex_v21d_ftrn(c,x21r)
  sdv = sqrt(dot_product(x21r,x21r))

  do i1=1,3
     call permut(i1,ihs)
     do i2=1,3
        vecdi(i2,:)=scc(ihs(i2),:)
     end do
     call drex_rotate_6x6_rot_ftrn(c,vecdi,cr)
     call drex_v21d_ftrn(cr,x21r)
     call drex_dec_proj(x21r,xp,xd,dev,5)
     if (dev < sdv) then         
        sdv=dev
        ndvc=i1
     endif
  end do
  
  vecdi = scc
  call permut(ndvc,ihs)
  do i1=1,3
     scc(i1,:)=vecdi(ihs(i1),:)
  end do
  !
  !21-D elastic vector in SCC axes
  !
  call drex_rotate_6x6_rot_ftrn(c,scc,cr)
  call drex_v21d_ftrn(cr,x21r)
end subroutine drex_dec_scca_drex


subroutine drex_dec_scca(c,vdi,vvo,scc,x21r,old_mode)
  implicit none
  !
  ! this version is from Jules Browaeys decmod routine
  ! as originally implemented in drex as of fstrack
  ! 
  ! modified by bug fix by Sebastian Chevrot 04/2011
  !
  ! additions due to Manabu, 07/2011
  !
  ! Calculates SCC directions as bissectrix
  ! of each close pair of dilatationnal and
  ! Voigt eigenvectors
  ! Higher symmetry axis = 3-axis
  !
  !t   = angles VO(1st index)/DI(2nd index)
  !bss = bissectrix direction intermediate
  !SCT = SCC test
  !dvm = array deviation for permutations
  !
  ! old_mode: 1: only find best-fit, and align
  ! old_mode: 0: find best and worst fit, and use those as a coordinate system
  !
  double precision, parameter :: pi=3.14159265358979323846d0
  !
  ! input voigt matrix
  !
  double precision, intent(in),dimension(6,6) :: c
  !
  ! input eigensystem of contractions
  !
  double precision, intent(in),dimension(3,3) :: vdi,vvo
  integer, intent(in) :: old_mode
  ! output
  double precision, intent(out),dimension(3,3) :: scc
  double precision, intent(out),dimension(21) :: x21r
  ! local
  integer :: i,j
  integer, dimension(2) :: pos
  integer, dimension(3) :: npos
  double precision :: dev
  double precision, dimension(3) :: bss,dvm
  double precision, dimension(3,3) :: t,SCT
  double precision, dimension(6,6) :: cr
  double precision, dimension(21) :: xp,xd
  !
  !VO/DI eigenvectors angle matrix
  !
  do i=1,3
     do j=1,3
        t(j,i) = dot_product(vdi(:,i),vvo(:,j))
     end do
  end do
  
  where (abs(t)>=1D0) t=sign(1D0,t)
  t=acos(t)
  where (t>0.5D0*pi) t=t-pi
  !Indice position
  !Association VO eigenvector to DI
  npos=minloc(abs(t),dim=1)
  !
  ! Calculates bissectrix vector & SCC
  !
  ! old version
  !do i=1,3
  !   bss=VDI(:,i)+sign(1D0,t(npos(i),i))*vvo(:,npos(i))
  !   scc(:,i)=bss/sqrt(dot_product(bss,bss))
  !end do
  !
  ! fix by sebastien chevrot, 04/2011
  !
  do i=1,3
     bss = VDI(:,1)+sign(1D0,t(npos(1),1))*vvo(:,npos(1))
  enddo
  scc(:,1) = bss/dsqrt(dot_product(bss,bss))
  bss = VDI(:,2)+sign(1D0,t(npos(2),2))*vvo(:,npos(2))
  bss = bss - dot_product(scc(:,1),bss)*scc(:,1)
  scc(:,2) = bss/dsqrt(dot_product(bss,bss))
  bss(1) = scc(2,1)*scc(3,2)-scc(3,1)*scc(2,2)
  bss(2) = scc(3,1)*scc(1,2)-scc(1,1)*scc(3,2)
  bss(3) = scc(1,1)*scc(2,2)-scc(2,1)*scc(1,2)
  scc(:,3) = bss/dsqrt(dot_product(bss,bss))
  
  
  !Transpose : basis transfer
  scc=transpose(scc)
  !Permutation 1,2,3
  npos=(/1,2,3/)
  !
  ! rotate c and convert to [21]
  call drex_rotate_6x6_rot_ftrn(c,scc,cr)
  call drex_v21d_ftrn(cr,x21r)
  ! project 
  call drex_dec_proj(x21r,xp,xd,dvm(1),5)
  !
  !Permutations 2,3,1 & 3,1,2
  !
  do j=2,3
     npos=cshift(npos,shift=1)
     do i=1,3
        SCT(i,:)=SCC(npos(i),:)
     end do
     call drex_rotate_6x6_rot_ftrn(c,sct,cr)
     call drex_v21d_ftrn(cr,x21r)
     ! project this permutation
     call drex_dec_proj(x21r,xp,xd,dvm(j),5)
  end do
  !
  ! choose the permutation for minimum deviation (this the old version)
  ! only
  !
  if(old_mode .eq. 1)then
     pos(1) = minloc(dvm,dim=1)
  ! rest of transverse isotropic fraction
     dev = dvm(pos(1))
     ! SCC coordinates system
     scc=cshift(scc,shift=pos(1)-1,dim=1)
  else
  !
  ! this the new version by Manabu
  !
     pos(1) = minloc(dvm,dim=1) ! for best fit case
     pos(2) = maxloc(dvm,dim=1) ! for worst fit case
     ! rest of transverse isotropic fraction
     dev = dvm(pos(1))
     !SCC coordinates system
     SCT=SCC
     !
     ! best fit axis
     !
     if(pos(1).eq.1) then
        SCC(3,:)=SCT(3,:)
     elseif (pos(1).eq.2) then
        SCC(3,:)=SCT(1,:)
     else
        SCC(3,:)=SCT(2,:)
     end if
     !
     ! worst fit axis
     !
     if(pos(2).eq.1) then
        SCC(1,:)=SCT(3,:)
     elseif (pos(2).eq.2) then
        SCC(1,:)=SCT(1,:)
     else
        SCC(1,:)=SCT(2,:)
     end if
     !
     ! intermediate fit axis ... defined as the cross product of SCC(1,:) (worst, x) and SCC(3,:) (best, z) 
     ! such that the x, worst, 1 - y, intermediate, 2 - z, best, 3 axis form a right-handed coordinate
     ! system
     !
     SCC(2,1)=SCC(3,2)*SCC(1,3)-SCC(3,3)*SCC(1,2)
     SCC(2,2)=SCC(3,3)*SCC(1,1)-SCC(3,1)*SCC(1,3)
     SCC(2,3)=SCC(3,1)*SCC(1,2)-SCC(3,2)*SCC(1,1)
  endif
  !
  !21-D elastic vector in SCC axes
  !
  call drex_rotate_6x6_rot_ftrn(c,scc,cr)
  call drex_v21d_ftrn(cr,x21r)
end subroutine drex_dec_scca



!
! compute symmetry projections
!
subroutine drex_dec_proj(x,xp,xd,dev,nsym)
  implicit none
  !
  !Symmetry projectors
  !
  !X    =  input 21-D vector
  !Xp   =  projected vector
  !XD   =  deviation vector
  !xdn  =  norm of deviation
  !NSYM = 13 monoclinic
  !     =  9 orthorhombic
  !     =  6 tetragonal
  !     =  5 hexagonal (transverse isotropic)
  !     =  2 isotropic
  integer, intent(in) :: nsym
  double precision, intent(in),dimension(21) :: x
  double precision, intent(out),dimension(21) :: xp,xd
  double precision, intent(out) :: dev
  ! local
  double precision :: sq2,isq2,i15
  parameter (sq2 = 1.41421356237309504d0) ! sqrt(2)
  isq2=1D0/sq2
  i15=1D0/15D0
  xp=0d0
  if (NSYM==13) then            ! monoclinic
     XP=X
     XP(10:11)=0D0
     XP(13:14)=0D0
     XP(16:17)=0D0
     XP(19:20)=0D0
  endif
  if (NSYM<=9) then
     XP(1:9)=X(1:9)
  endif
  if (NSYM==6) then
     XP(1)=0.5D0*(X(1)+X(2))
     XP(2)=XP(1)
     XP(4)=0.5D0*(X(4)+X(5))
     XP(5)=XP(4)
     XP(7)=0.5D0*(X(7)+X(8))
     XP(8)=XP(7)
  endif
  if (NSYM==5) then
     XP(1)=0.375D0*(X(1)+X(2))+0.25D0*X(6)*isq2+0.25D0*X(9)
     XP(2)=XP(1)
     XP(4)=0.5D0*(X(4)+X(5))
     XP(5)=XP(4)
     XP(6)=0.25D0*(X(1)+X(2))*isq2+0.75D0*X(6)-0.5D0*X(9)*isq2
     XP(7)=0.5D0*(X(7)+X(8))
     XP(8)=XP(7)
     XP(9)=0.25D0*(X(1)+X(2))-0.5D0*X(6)*isq2+0.5D0*X(9)
  endif
  if (NSYM==2) then
     XP(1)=3D0*(X(1)+X(2)+X(3))+sq2*(X(4)+X(5)+X(6))
     XP(1)=XP(1)+2D0*(X(7)+X(8)+X(9))
     XP(1)=XP(1)*i15
     XP(2)=XP(1)
     XP(3)=XP(1)
     XP(4)=sq2*(X(1)+X(2)+X(3))+4D0*(X(4)+X(5)+X(6))
     XP(4)=XP(4)-sq2*(X(7)+X(8)+X(9))
     XP(4)=XP(4)*i15
     XP(5)=XP(4)
     XP(6)=XP(4)
     XP(7)=2D0*(X(1)+X(2)+X(3))-sq2*(X(4)+X(5)+X(6))
     XP(7)=XP(7)+3D0*(X(7)+X(8)+X(9))
     XP(7)=XP(7)*i15
     XP(8)=XP(7)
     XP(9)=XP(7)
  endif
  xd = x - xp                       ! deviation
  ! norm of rest of projection
  dev = sqrt(dot_product(xd,xd))

end subroutine DREX_DEC_PROJ

