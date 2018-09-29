!       
!       aniosotropic reflectivity package from Vera Schulte-Pelkum in March 2005
!       based on code by Peter Shearer and various others, see individual READMEs 
!       and source code comments for details and description
!       
!       minor modifications by Thorsten Becker
!       
!       $Id: spectoseis.f90,v 1.3 2010/03/04 21:51:31 becker Exp becker $
!       
!       
!	program aniso_transient
!
!       takes output of anianza (eigenvectors and reflection/transmission matrices),
!       calculates time series
!       3.7.1998

!       added lh filter (same as used in pol)
!       17.3.2000

!       this version also calculates S response
!       26.10.2000

!       changed to pick upgoing eigenvector (second solution, negative vertical
!       slowness since z is positive down in anicake)
!       22.5.2001

!       added SKS frequency window and band command line argument 
!       4.6.2001

!       output file in ascii contains :
!       time, incident P - x,y,z, incident SV - x,y,z, incident SH - x,y,z
!       
!       x component is radial, -y component is transverse, since incidence in anicake
!       is in x-z plane

!
! this file has the main routines
!
! $Id: spectoseis.f90,v 1.3 2010/03/04 21:51:31 becker Exp becker $
!

subroutine spec_read_file(specfile,qp,qsv,qsh,qpi,qsvi,qshi,&
     ppmat,psvmat,pshmat,svpmat,svsvmat,svshmat,shpmat,shsvmat,shshmat,&
     imaxfreq,freq,nfreq)
  implicit none

  integer,intent(in) ::  imaxfreq,specfile ! max sample number, file pointer for input
  complex,intent(out) :: qp(3),qsv(3),qsh(3) !incidence vectors
  complex,intent(out) :: qpi(3),qsvi(3),qshi(3)
  integer,intent(out) :: nfreq !flip flags, number of samples


  !
  ! output 
  !
  complex, intent(out),dimension(1:imaxfreq) ::   ppmat,svpmat,shpmat,&
       psvmat,svsvmat,shsvmat,pshmat,svshmat,shshmat

  real, intent(out),dimension(1:imaxfreq) ::	freq
  !
  ! local
  !
  
  character (len=256)::		muell
  integer :: i,status
  real :: real1,imag1,real2,imag2,real3,imag3

  !	skip header and downgoing eigenvectors
  do i=1,13
     read(specfile,*)
  end do

  !	read upgoing eigenvectors (negative vertical slowness solutions)
  !	order is qP (re(x), im(x), re(y), im(y), re(z), im(z)), next line
  !	qS1, next line qS2
  read(specfile,*)  real1,imag1,real2,imag2,real3,imag3
  qp(1) = cmplx(real1,imag1)
  qp(2) = cmplx(real2,imag2)
  qp(3) = cmplx(real3,imag3)
  read(specfile,*)  real1,imag1,real2,imag2,real3,imag3
  qsv(1) = cmplx(real1,imag1)
  qsv(2) = cmplx(real2,imag2)
  qsv(3) = cmplx(real3,imag3)
  read(specfile,*)  real1,imag1,real2,imag2,real3,imag3
  qsh(1) = cmplx(real1,imag1)
  qsh(2) = cmplx(real2,imag2)
  qsh(3) = cmplx(real3,imag3)

  !	skip second layer eigenvectors and bottom layer downgoing eigenvectors
  do i = 1,11
     read(specfile,*)
  end do

  !	read upgoing initial eigenvectors (negative vertical slowness solutions)
  read(specfile,*)  real1,imag1,real2,imag2,real3,imag3
  qpi(1) = cmplx(real1,imag1)
  qpi(2) = cmplx(real2,imag2)
  qpi(3) = cmplx(real3,imag3)
  read(specfile,*)  real1,imag1,real2,imag2,real3,imag3
  qsvi(1) = cmplx(real1,imag1)
  qsvi(2) = cmplx(real2,imag2)
  qsvi(3) = cmplx(real3,imag3)
  read(specfile,*)  real1,imag1,real2,imag2,real3,imag3
  qshi(1) = cmplx(real1,imag1)
  qshi(2) = cmplx(real2,imag2)
  qshi(3) = cmplx(real3,imag3)

  !	skip travel times, free surface ref_u
  do i = 1,6
     read(specfile,*)
  end do


  nfreq = 0

  !	read reflection/transmission coefficients (incident/outgoing phase)
  do
     read(specfile,'(a5,f26.16)',iostat=status) muell, freq(nfreq+1)
     if (status<0) exit

     read(specfile,*) real1,imag1,real2,imag2,real3,imag3
     ppmat(nfreq+1)   = cmplx(real1,imag1)
     svpmat(nfreq+1)  = cmplx(real2,imag2)
     shpmat(nfreq+1)  = cmplx(real3,imag3)

     read(specfile,*) real1,imag1,real2,imag2,real3,imag3
     psvmat(nfreq+1)  = cmplx(real1,imag1)
     svsvmat(nfreq+1) = cmplx(real2,imag2)
     shsvmat(nfreq+1) = cmplx(real3,imag3)

     read(specfile,*) real1,imag1,real2,imag2,real3,imag3
     pshmat(nfreq+1)  = cmplx(real1,imag1)
     svshmat(nfreq+1) = cmplx(real2,imag2)
     shshmat(nfreq+1) = cmplx(real3,imag3)

     nfreq = nfreq+1
     if(nfreq .ge.imaxfreq)then
        print *,'spec_read_file: error: increase imaxfreq'
        stop
     endif
  end do

  return
end subroutine spec_read_file



!
!
! given the relevant matrices, compute seismograms and filter them
!
subroutine spec_compute_seis(bandphase,qp,qsv,qsh,qpi,qsvi,qshi,&
     ppmat,psvmat,pshmat,svpmat,svsvmat,svshmat,shpmat,shsvmat,shshmat,&
     pxtrace,pytrace,pztrace,svxtrace,svytrace,svztrace,shxtrace,shytrace,shztrace,&
     imaxfreq,freq,nfreq,ntime,dt,&
     timefile,print_orig_seis)
  implicit none
  ! input
  logical,intent(in) :: print_orig_seis !print unfiltered seismograms to timefile
  integer,intent(in) ::  nfreq,imaxfreq,timefile ! samples, max samples, file I/O
  character (len=256), intent(in) :: bandphase !filter type: P, SKS, SK2, SI
  ! incidence vectors
  complex,intent(in) ::	qp(3),qsv(3),qsh(3)
  complex,intent(in) ::	qpi(3),qsvi(3),qshi(3)
  ! matriced
  complex, intent(in),dimension(1:imaxfreq) ::   ppmat,svpmat,shpmat,&
       psvmat,svsvmat,shsvmat,pshmat,svshmat,shshmat
  ! frequence content
  real, dimension(1:imaxfreq),intent(in) ::	freq
  ! output 
  real, intent(out) :: dt
  integer, intent(out):: ntime
  real, dimension(1:imaxfreq), intent(out) ::	pxtrace,pytrace,pztrace,svxtrace,svytrace,svztrace,shxtrace,shytrace,shztrace
  !
  ! local
  !
  !	taper in f domain (1) or not (0) 
  integer, parameter :: taper = 1
  integer :: i,ierr
  ! flags
  integer :: dilatation,svflip,shflip

  !
  complex, dimension(1:imaxfreq) ::  pxspect,pyspect,pzspect
  complex, dimension(1:imaxfreq) ::  svxspect,svyspect,svzspect
  complex, dimension(1:imaxfreq) ::  shxspect,shyspect,shzspect
  !
  real, dimension(1:imaxfreq) :: af
  real ::				highpasslow,highpasshigh,df
  real ::				lowpasslow,lowpasshigh
  integer ::			highpasslength,lowpasslength


  !	check whether initial P is dilatation or compression
  dilatation = 0
  if(real(qpi(3))>0.) dilatation = 1

  !	check SV polarity
  svflip = 0
  if(real(qsvi(1))<0.) svflip = 1

  !	check SH polarity
  shflip = 0
  if(real(qshi(2))<0.) shflip = 1
  

  !
  df = freq(2)

  !	calculate spectra (eigenvectors times matrix coeffs)
  !	incident P wave
  pxspect = qp(1)*ppmat+qsv(1)*psvmat+qsh(1)*pshmat
  pyspect = qp(2)*ppmat+qsv(2)*psvmat+qsh(2)*pshmat
  pzspect = qp(3)*ppmat+qsv(3)*psvmat+qsh(3)*pshmat
  !	incident SV
  svxspect = qp(1)*svpmat+qsv(1)*svsvmat+qsh(1)*svshmat
  svyspect = qp(2)*svpmat+qsv(2)*svsvmat+qsh(2)*svshmat
  svzspect = qp(3)*svpmat+qsv(3)*svsvmat+qsh(3)*svshmat
  !	incident SH
  shxspect = qp(1)*shpmat+qsv(1)*shsvmat+qsh(1)*shshmat
  shyspect = qp(2)*shpmat+qsv(2)*shsvmat+qsh(2)*shshmat
  shzspect = qp(3)*shpmat+qsv(3)*shsvmat+qsh(3)*shshmat

  !	getts returns time series, # time points, dt
  call getts(pxspect, nfreq,df,taper,pxtrace, ntime,dt)
  call getts(pyspect, nfreq,df,taper,pytrace, ntime,dt)
  call getts(pzspect, nfreq,df,taper,pztrace, ntime,dt)
  call getts(svxspect,nfreq,df,taper,svxtrace,ntime,dt)
  call getts(svyspect,nfreq,df,taper,svytrace,ntime,dt)
  call getts(svzspect,nfreq,df,taper,svztrace,ntime,dt)
  call getts(shxspect,nfreq,df,taper,shxtrace,ntime,dt)
  call getts(shyspect,nfreq,df,taper,shytrace,ntime,dt)
  call getts(shzspect,nfreq,df,taper,shztrace,ntime,dt)

  !	change to z positive up coordinate system
  call flipsign(pytrace, ntime,imaxfreq)

  call flipsign(pztrace, ntime,imaxfreq)
  call flipsign(svytrace,ntime,imaxfreq)

  call flipsign(svztrace,ntime,imaxfreq)
  call flipsign(shytrace,ntime,imaxfreq)
  call flipsign(shztrace,ntime,imaxfreq)	


  !	flip all P components if initial P was dilatation
  if(dilatation.eq.1) then
     call flipsign(pxtrace,ntime,imaxfreq)
     call flipsign(pytrace,ntime,imaxfreq)
     call flipsign(pztrace,ntime,imaxfreq)
  endif

  !	correct SV polarity if necessary
  if(svflip.eq.1) then
     call flipsign(svxtrace,ntime,imaxfreq)
     call flipsign(svytrace,ntime,imaxfreq)
     call flipsign(svztrace,ntime,imaxfreq)
  endif

  !	correct SH polarity if necessary
  if(shflip.eq.1) then
     call flipsign(shxtrace,ntime,imaxfreq)
     call flipsign(shytrace,ntime,imaxfreq)
     call flipsign(shztrace,ntime,imaxfreq)
  endif
  if(print_orig_seis)then
  !       
  !       output of general time domain seismograms
  !       
     do i = 1,ntime
        write(timefile,('(10f26.16)')) (i-1)*dt, &
             pxtrace(i),pytrace(i),pztrace(i),&
             svxtrace(i),svytrace(i),svztrace(i),&  
             shxtrace(i),shytrace(i),shztrace(i)	  
     end do
  endif
  !	print*, 'nsamp ',ntime,' samprate ',dt
  ! 	apply lh filter (same as in pol)
  !	15 - 33 sec bandpass
  !	highpasslow = 0.03
  !	highpasslow = 0.01
  !	highpasshigh = 0.5/dt
  !	highpasslength = 121/dt+1
  !	lowpasslow = 0.0
  !	lowpasshigh = 0.067
  !	lowpasslength = 15/dt+1

  ! FILTERING

  if(bandphase(1:1)=='P') then
     !       write(0,*)'using P  filter'
     !       10 - 100 sec bandpass
     highpasslow = 0.01
     highpasshigh = 2.5
     highpasslength = 601
     lowpasslow = 0.0
     lowpasshigh = 0.1
     lowpasslength = 51
  else if(bandphase(1:3)=='SKS') then
     !       
     !       3 + 20 sec corners bandpass (center freq 7 sec, SKS band)
     !       
     !       0.05 Hz to 0.3 Hz, center 0.14 Hz
     !       highpasslow = 0.05
     !       highpasshigh = 2.5
     !       highpasslength = 601
     !       
     !       
     !       write(0,*)'using SKS  filter'
     highpasslow = 0.0
     highpasshigh = 0.3	! 3.3s
     highpasslength = 11
     lowpasslow = 0.0
     lowpasshigh = 0.1    ! 10s
     lowpasslength = 51
  else if(bandphase(1:2)=='SI') then
     !       
     !       modified SKS splitting, center period 12.5s, used for splitting intensity
     !       
     !       
     !	   write(0,*)'using splitting intensity filter'
     highpasslow = 0.066
     highpasshigh = 0.01
     highpasslength = 7.5/dt+1
     lowpasslow = 0.0
     lowpasshigh = 0.1
     lowpasslength = 51


     !	  highpasslow = 0.0
     !	  highpasshigh = 0.15
     !	  highpasslength = 3.5/dt+1
     !	  lowpasslow = 0.0
     !	  lowpasshigh = 0.05
     !	  lowpasslength = 3.5/dt+1
  else if(bandphase(1:3)=='SK2') then
     !
     !       modified SKS splitting, center period 12.5s
     !
     highpasslow = 0.0
     highpasshigh = 0.33 !  3s
     highpasslength = 3.5/dt+1
     lowpasslow = 0.0
     lowpasshigh = 0.05  ! 20s
     lowpasslength = 3.5/dt+1
  else if(bandphase(1:3)=='SK3') then
     !
     !       modified SKS splitting, center period 15s
     !
     highpasslow = 0.0
     highpasshigh = 0.1 !  10s
     highpasslength = 3.5/dt+1
     lowpasslow = 0.0
     lowpasshigh = 0.05  ! 20s
     lowpasslength = 3.5/dt+1
  else  
     return ! back to main if filter undefined
  endif
  !       check for odd half-lengths
  if(mod(highpasslength,2).eq.0) then
     highpasslength = highpasslength + 1
  endif
  if(mod(lowpasslength,2).eq.0) then
     lowpasslength = lowpasslength + 1
  endif

  call filter(pxtrace,af,ntime,dt,highpasslow,highpasshigh,highpasslength,ierr)
  call filter(af,pxtrace,ntime,dt,lowpasslow, lowpasshigh, lowpasslength, ierr)
  call filter(pytrace,af,ntime,dt,highpasslow,highpasshigh,highpasslength,ierr)
  call filter(af,pytrace,ntime,dt,lowpasslow, lowpasshigh, lowpasslength, ierr)
  call filter(pztrace,af,ntime,dt,highpasslow,highpasshigh,highpasslength,ierr)
  call filter(af,pztrace,ntime,dt,lowpasslow, lowpasshigh, lowpasslength, ierr)
  call filter(svxtrace,af,ntime,dt,highpasslow,highpasshigh,highpasslength,ierr)
  call filter(af,svxtrace,ntime,dt,lowpasslow, lowpasshigh, lowpasslength, ierr)
  call filter(svytrace,af,ntime,dt,highpasslow,highpasshigh,highpasslength,ierr)
  call filter(af,svytrace,ntime,dt,lowpasslow, lowpasshigh, lowpasslength, ierr)
  call filter(svztrace,af,ntime,dt,highpasslow,highpasshigh,highpasslength,ierr)
  call filter(af,svztrace,ntime,dt,lowpasslow, lowpasshigh, lowpasslength, ierr)
  call filter(shxtrace,af,ntime,dt,highpasslow,highpasshigh,highpasslength,ierr)
  call filter(af,shxtrace,ntime,dt,lowpasslow, lowpasshigh, lowpasslength, ierr)
  call filter(shytrace,af,ntime,dt,highpasslow,highpasshigh,highpasslength,ierr)
  call filter(af,shytrace,ntime,dt,lowpasslow, lowpasshigh, lowpasslength, ierr)
  call filter(shztrace,af,ntime,dt,highpasslow,highpasshigh,highpasslength,ierr)
  call filter(af,shztrace,ntime,dt,lowpasslow, lowpasshigh, lowpasslength, ierr)
  
end subroutine spec_compute_seis


subroutine flipsign(trace,npoints,imaxfreq)

  implicit none
  integer, intent(in) :: imaxfreq,npoints
  integer	::			ipoint
  real, dimension(1:imaxfreq) ::	trace


  do ipoint = 1, npoints
     trace(ipoint) = -trace(ipoint)
  end do

  return
end subroutine flipsign
	

!
! print seismograms to filtfile, and transverse/radial component to sifile if 
!
subroutine spec_print_seis(dt,imaxfreq,pxtrace,pytrace,pztrace,svxtrace,svytrace,svztrace,&
     shxtrace,shytrace,shztrace,ntime,filtfile)
  implicit none
  ! input
  integer, intent(in) :: ntime,filtfile,imaxfreq
  real, intent(in) :: dt
  real, intent(in),dimension(1:imaxfreq) ::	pxtrace,pytrace,pztrace
  real, intent(in),dimension(1:imaxfreq) ::	svxtrace,svytrace,svztrace
  real, intent(in),dimension(1:imaxfreq) ::	shxtrace,shytrace,shztrace
  
  integer :: i

  !       
  !       output of filtered seismograms
  !       
  do i = 1,ntime
     ! don't change the output format here, skssplit_xcorr relies on it
     write(filtfile,('(10f26.16)')) (i-1)*dt, pxtrace(i),pytrace(i),pztrace(i),&
          svxtrace(i),svytrace(i),svztrace(i),	  &
          shxtrace(i),shytrace(i),shztrace(i)	  
  end do
  
end subroutine spec_print_seis
!
! print  transverse/radial component 
!
subroutine spec_print_rt_seis(dt,ntime,imaxfreq,&
     svxtrace,svytrace,svztrace,sifile)
  implicit none
  ! input
  integer, intent(in) :: ntime,sifile,imaxfreq
  real, intent(in) :: dt
  real, intent(in),dimension(1:imaxfreq) ::	svxtrace,svytrace,svztrace
  integer :: i

  ! save the transverse and radial component for splitting intensity work
  write(sifile,('(2(e26.14,1x))'))dt,dt
  do i = 1,ntime
     !       transverse   radial
     write(sifile,('(2(e26.14,1x))'))-svytrace(i),-svxtrace(i)
  enddo
  close(sifile)
end subroutine spec_print_rt_seis
