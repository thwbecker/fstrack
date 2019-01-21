c
c aniosotropic reflectivity package from Vera Schulte-Pelkum in March 2005
c based on code by Peter Shearer and various others, see individual READMEs 
c and source code comments for details and description
c  
c
c $Id: skssplit_xcorr_driver.F,v 1.1 2011/04/12 06:19:51 becker Exp becker $
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
      integer ::		status,npts
      integer ::  null
      real    ::		svrpeak,svtpeak,svrpeaktime,svtpeaktime,dt,misfit
      real    ::		fastaz,delta,rayp,az,dum
      character*200  rayp_string,az_string
      real,dimension(1:NMAX_SAMPLE) ::	t,svr,svt,xcorrfn

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

      i = 1
      do
         ! need's to be same as in spectoseis output format
         read(iread,*,iostat=status) t(i),dum,dum,dum,svr(i),
     &        svt(i)
         if (status<0) exit
         i=i+1
      end do
      npts = i-1

      !
      ! compute cross correlation
      !
      
      call vera_split(t,svr,svt,npts,fastaz,delta,misfit,null)
c
      if (null.eq.1) then 
C     null
c         write(iout,'(f7.4,1x,f6.1,1x,a7,1x,a5)') rayp,az,'  null ',
c     &        ' 0.0 '
         write(iout,'(f7.4,1x,f6.1,1x,a7,1x,a5,1x,e14.7)') rayp,az,
     &        ' nan ',' 0.0 ',misfit
      else
C     split
         write(iout,'(f7.4,1x,f6.1,1x,f6.1,1x,f6.2,1x,e14.7)') 
     &        rayp,az,fastaz,delta,misfit
         
      end if
      if(iout.ne.6)then
         close(iout)      
      endif
c     stop
      end
