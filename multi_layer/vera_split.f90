!
!
! COMPUTE SPLITTING OF RADIAL/TRANSVERSE COMPONENT TIMESERIES BY MEANS OF CROSS-CORRELATION 
!
!
!
! aniosotropic reflectivity package from Vera Schulte-Pelkum in March 2005
! based on code by Peter Shearer and various others, see individual READMEs 
! and source code comments for details and description
!  
!
! $Id: vera_split.f90,v 1.2 2010/03/04 21:51:31 becker Exp becker $
!
! reads anicake-aniso_t_S output and calculates splitting parameters:
! grid search over rotation/delay using cross correlation
!
! looks for rotation angle and delay time that make orthogonal
! components identical mod scale (for single layer) or most similar (multiple
! layer splitting) - see levin & park 1999
!
! error calculation could be added by looking at width of cross correlation
! peak in time and az (not tested)
!
! cross correlation is done over window length winsec, currently hardwired to
! 10 s, centered on peak of radial component
!
! xcorr routine provided by Peter Shearer

! mar 2005 vsp
!
! minor modifications by Thorsten Becker
!
!
!


subroutine vera_split(t,svr,svt,npts,fastaz,deltat,misfit,null)
  implicit none
  !
  ! input
  !
  integer,intent(in) :: npts
  real, intent(in),dimension(1:npts) :: t,svr,svt ! time, radial, transverse components

  !
  ! output
  !
  real,intent(out) :: misfit,fastaz,deltat ! splitting parameters, fastaz in degree, delta in s
  integer,intent(out) :: null   ! null flag

  integer ::		i,ipeakr,ipeakt
  integer ::		nwin,nwinhalf,allnshiftmax,nshift
  integer ::		nshiftmax,nshiftmin
  real    ::		svrpeak,svtpeak,svrpeaktime,svtpeaktime,dt
  real    ::		winsec,alpha,alphainc,allcoeffmax,maxalpha
  real    ::		coeffmax,scalemax,coeffmin,scalemin,abscoeffmax
  real    ::		pi,dum,rms(2)
  real    ::            sin_alpha, cos_alpha,pi_fac,pi_half
  parameter 		(pi=3.14159265358979323846264d0,&
            pi_fac = 57.295779513082320876798154814105d0,&
            pi_half = 1.57079632679489661d0)
  real, dimension(:), allocatable :: fast,slow
  real, dimension(-1000:1000) :: xcorrfn
  
  if(npts.lt.40)then
     print *,'vera_xcorr: less than 40 samples?',npts
     stop
  endif
  
  allocate(fast(npts))
  allocate(slow(npts))


  svrpeak = 0.
  svtpeak = 0.
  svrpeaktime = -99.
  svtpeaktime = -99.
  ipeakr = 0
  ipeakt = 0
  
  rms = 0
  
  do i=1,npts
     if (abs(svr(i)) > svrpeak) then
        svrpeak = abs(svr(i))
        svrpeaktime = t(i)
        ipeakr = i
     end if
     if (abs(svt(i)) > svtpeak) then
        svtpeak = abs(svt(i))
        svtpeaktime = t(i)
        ipeakt = i
     end if
     rms(1) = rms(1) + svr(i)**2
     rms(2) = rms(2) + svt(i)**2
     
     if(i.eq.2)then
        dt = t(2) - t(1)
     else if(i.ne.1)then
        if(abs((t(i)-t(i-1))-dt).gt.1e-4)then
           print *,'vera_xcorr: dt error',t(i),t(i-1),dt,i,npts
           stop
        endif
     endif
  enddo
  rms = sqrt(rms/npts)
  

  !     test for null splitting (only spurious energy on transverse)
  null = 0
  if (svtpeak/svrpeak < 0.05) then
     null = 1
  end if
  
  
  !
  ! window length in seconds
  !
  winsec = 10.
  nwin = nint(winsec/dt)
  nwinhalf = nint(winsec/2 / dt)
  !      
  !     rotate components and cross correlate
  !
  alpha = -pi_half
  alphainc = 0.5/pi_fac     ! increment in radian
  allcoeffmax = 0.
  allnshiftmax = 0
  maxalpha = -99.
  do while (alpha < pi_half)
     sin_alpha = sin(alpha)
     cos_alpha = cos(alpha)
     slow = svr*sin_alpha+svt*cos_alpha;
     fast = svr*cos_alpha-svt*sin_alpha;
     call vera_getxcor(fast,slow,npts,npts,ipeakr-nwinhalf,&
          ipeakr-nwinhalf,nwin,nwinhalf,1,xcorrfn,coeffmax,nshiftmax,&
          scalemax,coeffmin,nshiftmin,scalemin)
     abscoeffmax = amax1(abs(coeffmax),abs(coeffmin))
     if (abs(coeffmax) > abs(coeffmin)) then
        nshift = nshiftmax
     else if (abs(coeffmax) < abs(coeffmin)) then
        nshift = nshiftmin
     end if
     if (abscoeffmax > allcoeffmax .and. nshift > 0) then
        allcoeffmax = abscoeffmax
        allnshiftmax = nshift
        maxalpha = alpha
     end if
!     write(*,'(e15.7,e15.7,i10,e15.7,e15.7,i10,e15.7)') 
!     &        alpha*180.0d0/pi,coeffmax,nshiftmax,
!     &        scalemax,coeffmin,nshiftmin,scalemin
     alpha = alpha + alphainc
  end do

  misfit=1.0d0-allcoeffmax
  if(allnshiftmax < 2)then
     null=1
     fastaz = 0.0d0
     deltat = 0.0d0
  else
     !     split
     fastaz = maxalpha * pi_fac
     deltat = allnshiftmax * dt
  endif

  deallocate(slow)
  deallocate(fast)

end subroutine vera_split
