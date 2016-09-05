c
c aniosotropic reflectivity package from Vera Schulte-Pelkum in March 2005
c based on code by Peter Shearer and various others, see individual READMEs 
c and source code comments for details and description
c  
c minor modifications by Thorsten Becker
c
c $Id: spectoseis.f,v 1.5 2005/07/28 22:06:37 becker Exp $
c
c
c	program aniso_transient
	
* takes output of anianza (eigenvectors and reflection/transmission matrices),
* calculates time series
* 3.7.1998

* added lh filter (same as used in pol)
* 17.3.2000

* this version also calculates S response
* 26.10.2000

* changed to pick upgoing eigenvector (second solution, negative vertical
* slowness since z is positive down in anicake)
* 22.5.2001

* added SKS frequency window and band command line argument 
* 4.6.2001

* output file in ascii contains :
* time, incident P - x,y,z, incident SV - x,y,z, incident SH - x,y,z
* x component is radial, y component is transverse, since incidence in anicake
* is in x-z plane

	implicit none
	
	complex ::			qpx,qsvx,qshx,qpy,qsvy,qshy,qpz,qsvz,qshz
	complex ::			qpxi,qsvxi,qshxi,qpyi,qsvyi,qshyi
	complex ::			qpzi,qsvzi,qshzi
	complex, dimension(1:20000) ::	ppmat,psvmat,pshmat,svpmat,svsvmat,svshmat
	complex, dimension(1:20000) ::  shpmat,shsvmat,shshmat
	complex, dimension(1:20000) ::  pxspect,pyspect,pzspect
	complex, dimension(1:20000) ::  svxspect,svyspect,svzspect
	complex, dimension(1:20000) ::  shxspect,shyspect,shzspect
	real ::				df,dt,real1,imag1,real2,imag2,real3,imag3
	real ::				highpasslow,highpasshigh
	real ::				lowpasslow,lowpasshigh
	real, dimension(1:20000) ::	freq,pxtrace,pytrace,pztrace,af
	real, dimension(1:20000) ::	svxtrace,svytrace,svztrace
	real, dimension(1:20000) ::	shxtrace,shytrace,shztrace
	integer ::			highpasslength,lowpasslength
	integer	::			nfreq,taper,ntime,i,ierr,iargc
	integer ::			specfile,status,timefile,filtfile,sifile
	integer ::			dilatation,svflip,shflip
	character (len=256)::		muell,infile, bandphase,save_rt_string
	logical ::                      save_rt
#ifdef USE_COMMAND_LINE_AND_FILES
!	command line arguments
	if(iargc().lt.2) then
	   print*, 'usage: spectoseis spectrumfile bandphase [save_rT]'
	   print*, '       calculates qP, qSV and qSH synthetics in stick (comb_out)'
	   print*, '       and filtered (aniso_filt_out) format'
	   print*, '       spectrumfile is output of anicake (cake_out)'
	   print*, '       bandphase is SKS (pulse width 6s) or P (band 10-100s)'
	   print*, '       save_rT will save radial and transverse component into file'
	   stop
	end if 	
	call getarg(1,infile)
	call getarg(2,bandphase)
	if(iargc().gt.2)then
	   call getarg(3,save_rt_string)
	   read(save_rt_string,*)i
	   if(i.ne.0)then 
	      save_rt = .true. 
	   else 
	      save_rt = .false.
	   endif
	else
	   save_rt = .false.
	endif
C       input 
 	specfile = 10
	open(specfile,file = infile)
c       time domain output
 	timefile = 11
	open(timefile,file = 'comb_out')
c       filtered output
	filtfile = 12
	open(filtfile,file = 'aniso_filt_out')
#else
	if(iargc().lt.1) then
	   print*, 'usage: spectoseis_steam bandphase'
	   print*, '       calculates qP, qSV and qSH synthetics in stick (comb_out)'
	   print*, '       and filtered (aniso_filt_out) format'
	   print*, '       output of anicake (cake_out) is read from stdin'
	   print*, '       bandphase is SKS (pulse width 6s) or P (band 10-100s)'
	   stop
	end if 	
	call getarg(1,bandphase)
	if(iargc().gt.1)then
	   call getarg(2,save_rt_string)
	   read(save_rt_string,*)i
	   if(i.ne.0)then 
	      save_rt = .true. 
	   else 
	      save_rt = .false.
	   endif
	else
	   save_rt = .false.
	endif

C       defaults:
C       get input from stdin
	specfile = 5		! stdin
!       timefile will be omitted,
!       filtered time file is stdout
	filtfile=6
#endif
	if(save_rt)then
	   sifile=20		! splitting intensity
	   open(sifile,file = 'sirt.out')
	endif

!	taper in f domain (1) or not (0) 
	taper = 1

	
*	skip header and downgoing eigenvectors
	do i=1,13
	  read(specfile,*)
	end do
	
*	read upgoing eigenvectors (negative vertical slowness solutions)
*	order is qP (re(x), im(x), re(y), im(y), re(z), im(z)), next line
*	qS1, next line qS2
	read(specfile,*)  real1,imag1,real2,imag2,real3,imag3
	qpx = cmplx(real1,imag1)
	qpy = cmplx(real2,imag2)
	qpz = cmplx(real3,imag3)
	read(specfile,*)  real1,imag1,real2,imag2,real3,imag3
	qsvx = cmplx(real1,imag1)
	qsvy = cmplx(real2,imag2)
	qsvz = cmplx(real3,imag3)
	read(specfile,*)  real1,imag1,real2,imag2,real3,imag3
	qshx = cmplx(real1,imag1)
	qshy = cmplx(real2,imag2)
	qshz = cmplx(real3,imag3)
	
*	skip second layer eigenvectors and bottom layer downgoing eigenvectors
	do i = 1,11
	  read(specfile,*)
	end do
	
*	read upgoing initial eigenvectors (negative vertical slowness solutions)
	read(specfile,*)  real1,imag1,real2,imag2,real3,imag3
	qpxi = cmplx(real1,imag1)
	qpyi = cmplx(real2,imag2)
	qpzi = cmplx(real3,imag3)
	read(specfile,*)  real1,imag1,real2,imag2,real3,imag3
	qsvxi = cmplx(real1,imag1)
	qsvyi = cmplx(real2,imag2)
	qsvzi = cmplx(real3,imag3)
	read(specfile,*)  real1,imag1,real2,imag2,real3,imag3
	qshxi = cmplx(real1,imag1)
	qshyi = cmplx(real2,imag2)
	qshzi = cmplx(real3,imag3)
	
*	check whether initial P is dilatation or compression
	dilatation = 0
	if(real(qpzi)>0.) dilatation = 1
	
*	check SV polarity
	svflip = 0
	if(real(qsvxi)<0.) svflip = 1
	
*	check SH polarity
	shflip = 0
	if(real(qshyi)<0.) shflip = 1
	
*	skip travel times, free surface ref_u
	do i = 1,6
	  read(specfile,*)
	end do
	

	nfreq = 0

*	read reflection/transmission coefficients (incident/outgoing phase)
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
	end do
!
!       END INPUT
!	

	df = freq(2)
	
*	calculate spectra (eigenvectors times matrix coeffs)
*	incident P wave
	pxspect = qpx*ppmat+qsvx*psvmat+qshx*pshmat
	pyspect = qpy*ppmat+qsvy*psvmat+qshy*pshmat
	pzspect = qpz*ppmat+qsvz*psvmat+qshz*pshmat
*	incident SV
	svxspect = qpx*svpmat+qsvx*svsvmat+qshx*svshmat
	svyspect = qpy*svpmat+qsvy*svsvmat+qshy*svshmat
	svzspect = qpz*svpmat+qsvz*svsvmat+qshz*svshmat
*	incident SH
	shxspect = qpx*shpmat+qsvx*shsvmat+qshx*shshmat
	shyspect = qpy*shpmat+qsvy*shsvmat+qshy*shshmat
	shzspect = qpz*shpmat+qsvz*shsvmat+qshz*shshmat
	
*	getts returns time series, # time points, dt
	call getts(pxspect, nfreq,df,taper,pxtrace, ntime,dt)
	call getts(pyspect, nfreq,df,taper,pytrace, ntime,dt)
	call getts(pzspect, nfreq,df,taper,pztrace, ntime,dt)
	call getts(svxspect,nfreq,df,taper,svxtrace,ntime,dt)
	call getts(svyspect,nfreq,df,taper,svytrace,ntime,dt)
	call getts(svzspect,nfreq,df,taper,svztrace,ntime,dt)
	call getts(shxspect,nfreq,df,taper,shxtrace,ntime,dt)
	call getts(shyspect,nfreq,df,taper,shytrace,ntime,dt)
	call getts(shzspect,nfreq,df,taper,shztrace,ntime,dt)
	
*	change to z positive up coordinate system
	call flipsign(pytrace, ntime)
	call flipsign(pztrace, ntime)
	call flipsign(svytrace,ntime)
	call flipsign(svztrace,ntime)
	call flipsign(shytrace,ntime)
	call flipsign(shztrace,ntime)	
	
*	flip all P components if initial P was dilatation
	if(dilatation.eq.1) then
	  call flipsign(pxtrace,ntime)
	  call flipsign(pytrace,ntime)
	  call flipsign(pztrace,ntime)
	endif
	
*	correct SV polarity if necessary
	if(svflip.eq.1) then
	  call flipsign(svxtrace,ntime)
	  call flipsign(svytrace,ntime)
	  call flipsign(svztrace,ntime)
	endif
	
*	correct SH polarity if necessary
	if(shflip.eq.1) then
	  call flipsign(shxtrace,ntime)
	  call flipsign(shytrace,ntime)
	  call flipsign(shztrace,ntime)
	endif
#ifdef USE_COMMAND_LINE_AND_FILES
!
! only for command line opration
!
! output of general time domain seismograms
!
	do i = 1,ntime
	  write(timefile,'(10f26.16)') (i-1)*dt, pxtrace(i),pytrace(i),
     .          pztrace(i),
     .		svxtrace(i),svytrace(i),svztrace(i),	  
     .		shxtrace(i),shytrace(i),shztrace(i)	  
	end do
#endif		
c	print*, 'nsamp ',ntime,' samprate ',dt
* 	apply lh filter (same as in pol)
*	15 - 33 sec bandpass
*	highpasslow = 0.03
c	highpasslow = 0.01
*	highpasshigh = 0.5/dt
*	highpasslength = 121/dt+1
*	lowpasslow = 0.0
*	lowpasshigh = 0.067
*	lowpasslength = 15/dt+1

	if(bandphase(1:1)=='P') then
*	  10 - 100 sec bandpass
	  highpasslow = 0.01
	  highpasshigh = 2.5
	  highpasslength = 601
	  lowpasslow = 0.0
	  lowpasshigh = 0.1
	  lowpasslength = 51
	elseif(bandphase(1:3)=='SKS') then
*	  3 + 20 sec corners bandpass (center freq 7 sec, SKS band)
*	  0.05 Hz to 0.3 Hz, center 0.14 Hz
c	  highpasslow = 0.05
c	  highpasshigh = 2.5
c	  highpasslength = 601
*	  try lowpass; corner at 6 s
	  highpasslow = 0.0
	  highpasshigh = 0.3
	  highpasslength = 11
	  lowpasslow = 0.0
	  lowpasshigh = 0.1
	  lowpasslength = 51
	else
	  print*, 'specify P or SKS for filter band in commmand line'
	  stop
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

!
!       output of filtered seismograms
!
	do i = 1,ntime
	   ! don't change the output format here, skssplit_xcorr relies on it
	   write(filtfile,'(10f26.16)') (i-1)*dt, pxtrace(i),pytrace(i),pztrace(i),
     .		svxtrace(i),svytrace(i),svztrace(i),	  
     .		shxtrace(i),shytrace(i),shztrace(i)	  
	end do

	if(save_rt)then		! save the transverse and radial component for splitting intensity work
	   write(sifile,'(2(e26.14,1x))')dt,dt
	   do i = 1,ntime
C                                           transverse   radial
	      write(sifile,'(2(e26.14,1x))')-svytrace(i),svxtrace(i)
	   enddo
	   close(sifile)
	endif
	  
c	stop
	end


	subroutine flipsign(trace,npoints)
	
	implicit none
	
	real, dimension(1:20000) ::	trace
	integer	::			npoints, ipoint
	
	do ipoint = 1, npoints
	  trace(ipoint) = -1. * trace(ipoint)
	end do
	
	return
	end
	
	
