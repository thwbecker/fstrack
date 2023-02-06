subroutine drex_dec_scca(c,vdi,vvo,scc,x21r)
  implicit none
  !
  ! this version is from Jules Browaeys decmod routine
  ! as originally implemented in drex as of fstrack
  ! 
  ! modified by bug fix by Sebastian Chevrot 04/2011
  !
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
  integer :: i,j
  integer, dimension(2) :: pos ! 1 for best fit case and 2 for worst fit case
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

  pos(1) = minloc(dvm) ! for best fit case
  pos(2) = maxloc(dvm) ! for worst fit case
! rest of transverse isotropic fraction
  dev = dvm(pos(1))
!SCC coordinates system
  SCT=SCC
! best fit axis
  if(pos(1)==1) then
	SCC(3,:)=SCT(3,:)
  else(pos(1)==2) then
	SCC(3,:)=SCT(1,:)
  else
	SCC(3,:)=SCT(2,:)
  end if
! worst fit axis
  if(pos(2)==1) then
	SCC(1,:)=SCT(3,:)
  else(pos(2)==2) then
	SCC(1,:)=SCT(1,:)
  else
	SCC(1,:)=SCT(2,:)
  end if
! intermediate fit axis ... defined as the cross product of SCC(3,:) and SCC(1,:)
	SCC(2,1)=SCC(3,2)*SCC(1,3)-SCC(3,3)*SCC(1,2)
	SCC(2,2)=SCC(3,3)*SCC(1,1)-SCC(3,1)*SCC(1,3)
	SCC(2,3)=SCC(3,1)*SCC(1,2)-SCC(3,2)*SCC(1,1)
  ! 
  !21-D elastic vector in SCC axes
  !
  call drex_rotate_6x6_rot_ftrn(c,scc,cr)
  call drex_v21d_ftrn(cr,x21r)
end subroutine drex_dec_scca
