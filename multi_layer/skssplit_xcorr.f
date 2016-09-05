c
c aniosotropic reflectivity package from Vera Schulte-Pelkum in March 2005
c based on code by Peter Shearer and various others, see individual READMEs 
c and source code comments for details and description
c  
c
c $Id: skssplit_xcorr.f,v 1.3 2005/04/04 02:24:29 becker Exp becker $
c
c
c      program skssplit_xcorr
      
* reads anicake-aniso_t_S output and calculates splitting parameters:
* grid search over rotation/delay using cross correlation
c
* looks for rotation angle and delay time that make orthogonal
* components identical mod scale (for single layer) or most similar (multiple
* layer splitting) - see levin & park 1999
c
* error calculation could be added by looking at width of cross correlation
* peak in time and az (not tested)
c
* cross correlation is done over window length winsec, currently hardwired to
* 10 s, centered on peak of radial component
c
* xcorr routine provided by Peter Shearer

* mar 2005 vsp

c
c minor modifications by Thorsten Becker
c
c in regular version (compiled with -DUSE_COMMAND_LINE_AND_FILES)
c
* input is aniso_filt_out
* output in skssplit_out
c
c in stream_ref version, compiled by NOT SETTING 
c USE_COMMAND_LINE_AND_FILES
c
* input is stdin (aniso_filt_out)
* output in stdout (skssplit_out)
c
c
c
c
c
#define NMAX_SAMPLE 16384
      implicit none
      
      integer ::		iread,iaz,iout,i,ipeakr,ipeakt,iargc
      integer ::		nwin,nwinhalf,allnshiftmax,nshift
      integer ::		nshiftmax,nshiftmin,npts
      integer ::		status,null
      real    ::		svrpeak,svtpeak,svrpeaktime,svtpeaktime,dt
      real    ::		winsec,alpha,alphainc,allcoeffmax,maxalpha
      real    ::		coeffmax,scalemax,coeffmin,scalemin,abscoeffmax
      real    ::		fastaz,deltat,rayp,az,pi,dum
      real    ::                sin_alpha, cos_alpha,pi_fac,pi_half
      character*200  rayp_string,az_string
      real,dimension(1:NMAX_SAMPLE) ::	t,svr,svt,xcorrfn,fast,slow
      parameter 		(pi=3.14159265358979323846264d0,
     &     pi_fac = 57.295779513082320876798154814105d0,
     &     pi_half = 1.57079632679489661d0)

#ifdef USE_COMMAND_LINE_AND_FILES
! input file
      iread = 10
      open(iread,file='aniso_filt_out')
! output file
      iout  = 12
      open(iout,file='skssplit_out',position='append') 
! read ray parameter and azimuth from file (those are just echoed)                        
      iaz   = 11
      open(iaz,file='incaz.out')
      read(iaz,*) rayp, az
      close(iaz)
#else
! stream version, get paramters from STDIN
      if(iargc()/=2) then
         print*, 'usage: skssplit_xcorr_stream rayp az'
         print*, '       calculates SKS splitting'
         stop
      end if 	
      call getarg(1,rayp_string)
      call getarg(2,az_string) 
      read(rayp_string,*)rayp
      read(az_string,*)az

      iread = 5             ! stdin
      iout = 6                  ! stdout
#endif
c      
*     read to end of synthetics file; read horizontal components for
*     incident SV (SKS/SKKS); determine peak times and amplitudes
c
      svrpeak = 0.
      svtpeak = 0.
      svrpeaktime = -99.
      svtpeaktime = -99.
      ipeakr = 0
      ipeakt = 0

      i = 1
      do
         ! need's to be same as in spectoseis output format
        read(iread,'(10f26.16)',iostat=status) t(i),dum,dum,dum,svr(i),
     &        svt(i)
	if (status<0) exit
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
	i = i + 1
      end do
      npts = i
      if(npts.gt.NMAX_SAMPLE)then
         print *,'skssplit_xcorr: too many samples'
         print *,npts,NMAX_SAMPLE
         stop
      endif
*     test for null splitting (only spurious energy on transverse)
      null = 0
      if (svtpeak/svrpeak < 0.05) then
        null = 1
      end if
c      print*, svtpeaktime,svrpeaktime,svtpeak,svrpeak,null
            
c      print*, svtpeaktime,svrpeaktime

      dt = t(2)
      winsec = 10.
      nwin = nint(winsec/dt)
      nwinhalf = nint(winsec/2 / dt)
c      
*     rotate components and cross correlate
c
      alpha = -pi_half
      alphainc = 0.5/pi_fac     ! increment in radian
      allcoeffmax = 0.
      allnshiftmax = 0
      maxalpha = -99.
      do while (alpha < pi_half .and. null.eq.0)
         sin_alpha = sin(alpha)
         cos_alpha = cos(alpha)
         slow = svr*sin_alpha+svt*cos_alpha;
         fast = svr*cos_alpha-svt*sin_alpha;
         call getxcor(fast,slow,npts,npts,ipeakr-nwinhalf,
     &        ipeakr-nwinhalf,
     &        nwin,nwinhalf,1,xcorrfn,coeffmax,nshiftmax,
     &        scalemax,
     &        coeffmin,nshiftmin,scalemin)
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
c     write(*,'(e15.7,e15.7,i10,e15.7,e15.7,i10,e15.7)') 
c     &        alpha*180.0d0/pi,coeffmax,nshiftmax,
c     &        scalemax,coeffmin,nshiftmin,scalemin
         alpha = alpha + alphainc
      end do

c
c
      if (null.eq.1 .or. allnshiftmax<2) then 
C     null
c         write(iout,'(f6.4,1x,f6.1,1x,a7,1x,a5)') rayp,az,'  null ',
c     &        ' 0.0 '
         write(iout,'(f6.4,1x,f6.1,1x,a7,1x,a5)') rayp,az,' nan ',
     &        ' 0.0 '
      else
C     split
         fastaz = maxalpha * pi_fac
         deltat = allnshiftmax * dt
         write(iout,'(f6.4,1x,f6.1,1x,f6.1,1x,f6.2)') 
     &        rayp,az,fastaz,deltat      
         
      end if
      if(iout.ne.6)then
         close(iout)      
      endif
c     stop
      end
