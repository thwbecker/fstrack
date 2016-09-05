!
! utility routines to convert stiffness tensors to 
! splitting observations
!
! code from Vera Schulte-Pelkum as of 02/2005
!
! minor modification by TWB 
!
! $Id: vera_util.f90,v 1.8 2010/12/31 18:59:54 becker Exp becker $
!
!
! main routines: vera_slowness and vera__split_from_tensor
!

!
! obtain the fast S wave splitting from a stiffness tensor
! at incidence and back_azimuth angles [deg] with layer thickness 
! laythick[km]
!
! output: faz, dt  (fast azimuth [deg], delay time [s])
!         sfastx,sfasty: x and y components of normalized 
!                         slowness vector
!         vsphase(2): fast and slow S wave velocities
!
subroutine vera_layer_split_from_tensor_ftrn(cijkl,incidence,azimuth,&
     laythick,sfaste,sfastn,faz,dt,vsphase)
  implicit none
  ! INPUT:
  double precision, intent(in) :: cijkl(3,3,3,3) ! stiffness tensor normalized by
                                 ! density
  double precision, intent(in) :: incidence, azimuth ! of wave in [deg]
  double precision, intent(in) :: laythick  ! layer thickness in [km]
  ! OUTPUT:
  double precision, intent(out) :: faz,dt,sfaste,sfastn
  double precision, intent(out) :: vsphase(2)
  ! INTERNAL
  double precision :: pif,qh,rinc,razi,cos_inc,&
       sfastlen
  ! intermediate velocity parameters
  double precision, dimension (3):: slowmag,groupmag,qinc
  double precision, dimension (3,3) :: slowvec,polvec,groupvec
  ! other
  integer :: fastcomp,slowcomp
  ! parameters
  parameter (pif = 57.2957795130823d0) ! 180/pi
  !
  !	    build unit vector of incident slowness
  !	    length q = 1; sin inc = horizontal projection of q
  !	    minus signs since we are looking at incoming waves from center
  !	    of unit sphere
  !
  
  rinc =  incidence/pif           ! convert to radians
  razi =  azimuth/pif        
  
  qh = sin(rinc)
  cos_inc = cos(rinc)           ! cos(incidence)
  !
  ! the following components are for the negative
  ! incidence vector
  !
  !	    z component
  qinc(3) = - cos_inc
  !	    x component (i.e. South) for azimuth azi is -cos(azi)
  qinc(1) =   qh * cos(razi)
  !	    y component (i.e. East) is sin(azi)
  qinc(2) = - qh * sin(razi)
  !	    solve for phase, group velocities, polarizations
  !	    results are sorted as (component 1-3, phase qP, qSV, qSH)


  call vera_slowness(cijkl,qinc,slowmag,slowvec,polvec,&
       groupmag,groupvec)
  !
  !	    find fast S wave (qSV or qSH)
  !
  if(slowmag(2) < slowmag(3)) then
     fastcomp = 2
     slowcomp = 3
  else
     fastcomp = 3
     slowcomp = 2
  end if
  
  !	    save x and y direction of fast S polarization
  sfastlen = sqrt(polvec(1,fastcomp)**2+polvec(2,fastcomp)**2)
  if(abs(sfastlen).lt.1e-8)then
     print *,'vera_layer_split_from_tensor: error'
     print *,'fast velocity is near zero:',sfastlen
     stop
  endif
  ! 
  ! changed from Vera's convention to expected Cartesian 
  ! system where x: south, y: east, z: up
  !
  sfaste =  polvec(2,fastcomp)/sfastlen
  sfastn = -polvec(1,fastcomp)/sfastlen
  
  !	    fast and slow S phase velocities
  !
  vsphase(1) = 1.0d0/slowmag(fastcomp) !fast
  vsphase(2) = 1.0d0/slowmag(slowcomp) !slow
  
  !	    calculate delay time in layer
  if(cos_inc.eq.0)then
     print *,'vera_layer_split_from_tensor: error'
     print *,'incidence: ',incidence
     stop
  endif
  ! delay time is travel time slow - travel time fast
  dt = laythick / cos_inc * (1.d0/vsphase(2) - 1.0d0/vsphase(1))
  
  !         fast azimuth, in deg
  faz = pif * atan2(sfaste,sfastn)
  if(faz.lt.0)then 
     faz = faz + 360.d0
  endif
  ! limit to [0;180]
  if(faz.gt.180)then
     faz = faz - 180.0
  endif
end subroutine vera_layer_split_from_tensor_ftrn

!
!	    solve for phase, group velocities, polarizations
!	    results are sorted as (component 1-3, phase qP, qSV, qSH)
! input: c(3,3,3,3) stiffness / density tensor
!        qhat(3) incident slowness vector
! 
!
!
! Subroutine vera_slowness solves for the slowness, polarization,
! and group velocity vectors for a specified phase direction.
! For the case of general anisotropy, this is solved as an
! eigenvalue problem.
!  Given:
!    c(3,3,3,3) = elastic tensor
!    qhat(3)    = phase (slowness) direction.  This is a unit vector.
!  Returned:
!    qmag(3) = slowness magnitudes (eigenvalues)
!    q(3,3)  = slowness vectors
!    a(3,3)  = polarization vectors (eigenvectors)
!    umag(3) = group velocity magnitudes
!    u(3,3)  = group velocity vectors
! All eigenvalues and eigenvectors are sorted as
!      1) qP
!      2) qSV
!      3) qSH
!
!
! x(1) = east, x(2) = north, x(3) up
!
subroutine vera_slowness(c,qhat,qmag,q,a,umag,u)
  implicit none 
  double precision :: c(3,3,3,3),qhat(3),qmag(3),q(3,3),a(3,3),&
       umag(3),u(3,3),err,dotcor,emax,amax,&
       m(3,3),h(3,3),eigen(3,3),ev(3)
  integer :: index(3),i,j,k,l,n

  parameter(err=1e-7)          ! error condition for Jacobi
                                ! (this used to be 1e-5)
  do  i=1,3
     do k=1,3
        m(i,k)=0.d0
        do  j=1,3
           do l=1,3
              m(i,k)=m(i,k)+c(i,j,k,l)*qhat(j)*qhat(l)
           enddo
           h(i,k)=m(i,k)
        enddo
     enddo
  enddo

  !
  ! get eigensystem, eigenvalues are the diagonal of h()
  !
  call vera_jacobi(h,eigen,3,3,err)
  !

  ! Now sort by polarization
  ! qP is largest eigenvalue (smallest slowness)
  emax=0.d0
  do i=1,3
     if (h(i,i).gt.emax) then
        emax=h(i,i)
        index(1)=i
     end if
  enddo
  ! qSV has most vertical (001) polarization
  amax=0.d0
  do i=1,3
     if (i.ne.index(1)) then
        if (abs(eigen(3,i)).ge.amax) then
           amax=abs(eigen(3,i))
           index(2)=i
        end if
     endif
  enddo
  
  ! qSH is the remaining polarization
  do i=1,3
     if (i.ne.index(1).and.i.ne.index(2)) then
        index(3)=i
     end if
  enddo

  do i=1,3
     ev(i) = h(index(i),index(i))
     if (ev(i).lt.0.) then
        print *,'vera_slowness: We are going to bomb!!!!!'
        print *,'vera_slowness: i,index= ',i,index
        print *,'vera_slowness: h= ',h
     end if
     if (ev(i).ne.0.) then
        qmag(i)=1.d0/sqrt(ev(i))
     else
        qmag(i)=0.d0
     end if
     do j=1,3
        q(j,i)=qhat(j)*qmag(i)
        a(j,i)=eigen(j,index(i))
     enddo
  enddo
!      

! Now find group velocity vectors
  do  i=1,3
     dotcor=0.d0
     do  j=1,3
        u(j,i)=0.d0
        do  k=1,3
           do  l=1,3
              do  n=1,3
                 u(j,i)=u(j,i)+c(j,k,l,n)*a(k,i)*q(l,i)*a(n,i)
              enddo
           enddo
        enddo
        dotcor=dotcor+u(j,i)*q(j,i)
     enddo
     umag(i)=0.d0
     do j=1,3
        if (dotcor.ne.0.) u(j,i)=u(j,i)/dotcor
        umag(i)=umag(i)+u(j,i)**2
     enddo
     umag(i)=sqrt(umag(i))
  enddo
  
  return
end subroutine vera_slowness
    

!
! convert a non-normalized, upper triangular stiffness tensor
! sav[6,6] in C format (not all elements filled) into a normalized 
! Cijkl tensor [3,3,3,3] which will be returned Fortran style
!
! input: sav, rho
! output: sij
! sav is unchanged
!
subroutine vera_sav_to_cijkl_ftrn(sav,rho,cij)
  implicit none
  double precision :: sav(36),c(6,6),cij(3,3,3,3), &
       rho
  integer :: i,j
! resort from C to fortran style, upper right hand
  do i=1,6
     do j=i,6
        c(i,j) = sav((i-1)*6+j)
     enddo
  enddo
! fill in symmetric parts and normalize
  call vera_fill_and_normalize_c2(c,rho)
! convert to cijkl
  call vera_c2c4(c,cij)
end subroutine vera_sav_to_cijkl_ftrn


! print splitting parameters to filep
!
subroutine vera_print_splitting_ftrn(inc,az,sfastx,sfasty,faz,dt,filep)
  implicit none
  double precision :: inc,az,faz,dt,sfastx,sfasty
  integer :: filep
  !	    write inc, az, sfastx, sfasty, dt(layer) to file 
  write(filep,'(f5.1,1x,f5.1,1x,f10.7,1x,f10.7,1x,f8.3,1x,f8.5)') &
       inc,az,sfastx,sfasty,faz,dt
end subroutine vera_print_splitting_ftrn

! given a c[6,6] tensor where the upper triangle of elements
! is filled in with elements, fill in the symmetric lower half
! and normalize by density
!
! based on code anisotropy written by Dan Raymer, Leeds University.
! portion to put Cij into Cijkl format
! 30.10.2000
!
! this version normalises by density (e.g. 3.353 g/cm^3 - olivine)
! input: Cij in GPa, not normalised by density; order 11,22,..,12,13,..,23,..
! output: cijkl in GPa/density[g/cm^3] for input in christoffel equation solver
!         tensor_phase
!
subroutine vera_fill_and_normalize_c2(c,rho)
  implicit none
  double precision :: c(6,6),rho
  integer :: i,j
  do i=1,6
     do j=i,6
        c(i,j)=c(i,j)/rho
        c(j,i)=c(i,j)
     enddo
  enddo
  return
end subroutine vera_fill_and_normalize_c2


!
! convert a c[6,6] matrix into a c[3,3,3,3] Cijkl
! tensor
!
subroutine vera_c2c4(c,cc)   
  implicit none
  double precision ::  c(6,6),cc(3,3,3,3)
  
  cc = 0.d0

  cc(1,1,1,1) = c(1,1)
  cc(2,2,2,2) = c(2,2)
  cc(3,3,3,3) = c(3,3)
  cc(2,3,2,3) = c(4,4)
  cc(3,2,3,2) =cc(2,3,2,3)
  cc(2,3,3,2) =cc(2,3,2,3)
  cc(3,2,2,3) =cc(2,3,2,3)
  cc(1,3,1,3) = c(5,5)
  cc(3,1,1,3) =cc(1,3,1,3)
  cc(1,3,3,1) =cc(1,3,1,3)
  cc(3,1,3,1) =cc(1,3,1,3)
  cc(1,1,2,2) = c(1,2)
  cc(2,2,1,1) =cc(1,1,2,2)
  cc(1,1,3,3) = c(1,3)
  cc(3,3,1,1) =cc(1,1,3,3)
  cc(1,1,2,3) = c(1,4)
  cc(1,1,3,2) =cc(1,1,2,3)
  cc(2,3,1,1) =cc(1,1,2,3)
  cc(3,2,1,1) =cc(1,1,2,3)
  cc(1,1,1,3) = c(1,5)
  cc(1,1,3,1) =cc(1,1,1,3)
  cc(1,3,1,1) =cc(1,1,1,3)
  cc(3,1,1,1) =cc(1,1,1,3)
  cc(1,1,1,2) = c(1,6)
  cc(1,1,2,1) =cc(1,1,1,2)
  cc(1,2,1,1) =cc(1,1,1,2)
  cc(2,1,1,1) =cc(1,1,1,2)
  cc(2,2,3,3) = c(2,3)
  cc(3,3,2,2) =cc(2,2,3,3)
  cc(2,2,2,3) = c(2,4)
  cc(2,2,3,2) =cc(2,2,2,3)
  cc(2,3,2,2) =cc(2,2,2,3)
  cc(3,2,2,2) =cc(2,2,2,3)
  cc(2,2,1,3) = c(2,5)
  cc(2,2,3,1) =cc(2,2,1,3)
  cc(1,3,2,2) =cc(2,2,1,3)
  cc(3,1,2,2) =cc(2,2,1,3)
  cc(2,2,1,2) = c(2,6)
  cc(2,2,2,1) =cc(2,2,1,2)
  cc(1,2,2,2) =cc(2,2,1,2)
  cc(2,1,2,2) =cc(2,2,1,2)
  cc(3,3,2,3) = c(3,4)
  cc(3,3,3,2) = cc(3,3,2,3)
  cc(2,3,3,3) = cc(3,3,2,3)
  cc(3,2,3,3) = cc(3,3,2,3)
  cc(3,3,1,3) = c(3,5)
  cc(3,3,3,1) = cc(3,3,1,3)
  cc(1,3,3,3) = cc(3,3,1,3)
  cc(3,1,3,3) = cc(3,3,1,3)
  cc(3,3,1,2) = c(3,6)
  cc(3,3,2,1) = cc(3,3,1,2)
  cc(1,2,3,3) = cc(3,3,1,2)
  cc(2,1,3,3) = cc(3,3,1,2)
  cc(2,3,1,3) = c(4,5)
  cc(3,2,1,3) =cc(2,3,1,3)
  cc(1,3,3,2) =cc(2,3,1,3)
  cc(1,3,2,3) =cc(2,3,1,3)
  cc(2,3,3,1) =cc(2,3,1,3)
  cc(3,2,3,1) =cc(2,3,1,3)
  cc(3,1,2,3) =cc(2,3,1,3)
  cc(3,1,3,2) =cc(2,3,1,3)
  cc(2,3,1,2) = c(4,6)
  cc(3,2,1,2) =cc(2,3,1,2)
  cc(1,2,2,3) =cc(2,3,1,2)
  cc(1,2,3,2) =cc(2,3,1,2)
  cc(2,3,2,1) =cc(2,3,1,2)
  cc(3,2,2,1) =cc(2,3,1,2)
  cc(2,1,2,3) =cc(2,3,1,2)
  cc(2,1,3,2) =cc(2,3,1,2)
  cc(1,3,1,2) = c(5,6)
  cc(3,1,1,2) =cc(1,3,1,2)
  cc(1,2,1,3) =cc(1,3,1,2)
  cc(1,2,3,1) =cc(1,3,1,2)
  cc(1,3,2,1) =cc(1,3,1,2)
  cc(3,1,2,1) =cc(1,3,1,2)
  cc(2,1,1,3) =cc(1,3,1,2)
  cc(2,1,3,1) =cc(1,3,1,2)
  cc(1,2,1,2) = c(6,6)
  cc(2,1,1,2) =cc(1,2,1,2)
  cc(1,2,2,1) =cc(1,2,1,2)
  cc(2,1,2,1) =cc(1,2,1,2)
  return
end subroutine vera_c2c4

!
! convert a c[3,3,3,3] Cijkl tensor
! to c[6,6] matrix
!
subroutine vera_c4c2(cc,c)   
  implicit none
  double precision ::  c(6,6),cc(3,3,3,3)
  INTEGER :: i,j
  
  c(1,1) = cc(1,1,1,1) 
  c(2,2) = cc(2,2,2,2) 
  c(3,3) = cc(3,3,3,3) 
  c(4,4) = cc(2,3,2,3) 
  c(5,5) = cc(1,3,1,3) 
  c(1,2) = cc(1,1,2,2) 
  c(1,3) = cc(1,1,3,3) 
  c(1,4) = cc(1,1,2,3) 
  c(1,5) = cc(1,1,1,3) 
  c(1,6) = cc(1,1,1,2) 
  c(2,3) = cc(2,2,3,3) 
  c(2,4) = cc(2,2,2,3) 
  c(2,5) = cc(2,2,1,3) 
  c(2,6) = cc(2,2,1,2) 
  c(3,4) = cc(3,3,2,3) 
  c(3,5) = cc(3,3,1,3) 
  c(3,6) = cc(3,3,1,2) 
  c(4,5) = cc(2,3,1,3) 
  c(4,6) = cc(2,3,1,2) 
  c(5,6) = cc(1,3,1,2) 
  c(6,6) = cc(1,2,1,2) 
  do i=1,6
     do j=i,6
        c(j,i) = c(i,j)
     enddo
  enddo
  return
end subroutine vera_c4c2

! print 6,6 in thorsten's sav format
subroutine vera_print_c2_ftrn(c)   
  implicit none
  double precision ::  c(6,6)
  INTEGER :: i,j
  write(6,'(24(f5.1,1x))')0,0,0,((c(i,j),&
       j=i,6),i=1,6)
  return 
end subroutine vera_print_c2_ftrn

!
! read Cijkl tensor from file pointer fp
!
subroutine vera_read_cijkl(c,fp)
  implicit none
  integer fp,i,j,k,l
  double precision, dimension (3,3,3,3) :: c
  read(fp,*) &
     ((((c(i,j,k,l),i=1,3),j=1,3),k=1,3),l=1,3)

end subroutine vera_read_cijkl

      

!
! print Cijkl tensor to filep
!
subroutine vera_print_cijkl_ftrn(cij,filep)
  implicit none
  double precision :: cij(3,3,3,3)
  integer :: i,j,k,l,filep
  write(filep,'(1p,9(e14.6,1x))') ((((cij(i,j,k,l),&
       i=1,3),j=1,3),k=1,3),l=1,3)
  
end subroutine vera_print_cijkl_ftrn

! print to file
subroutine vera_print_cijkl_file_ftrn(cc)   
  implicit none
  double precision ::  cc(3,3,3,3)
  integer :: fp
  fp = 99
  open(fp,file = 'tmp.cijkl.dat')
  call vera_print_cijkl_ftrn(cc,fp)
  close(fp)
end subroutine vera_print_cijkl_file_ftrn
