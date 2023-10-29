!       
!       aniosotropic reflectivity package from Vera Schulte-Pelkum in March 2005
!       based on code by Peter Shearer and various others, see individual READMEs 
!       and source code comments for details and description
!       
!       minor modifications by Thorsten Becker
!       
!       $Id: spectoseis_driver.f90,v 1.2 2005/11/18 19:51:19 becker Exp $
!       
!       
!	program aniso_transient

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
!
!       modularized by twb
!
!
!     main routines are in spectoseis.f90, this file has the drivers and can be compiled 
!     to yield a strem (stdin/stdout) or file I/O based version
!
!    $Id: spectoseis_driver.f90,v 1.2 2005/11/18 19:51:19 becker Exp $ 
!
#include "../fstrack_ftrn.h"

program main
!
integer,parameter :: imaxfreq=20000 ! max number of samples
!
complex,dimension(3) ::			qp,qsv,qsh,qpi,qshi,qsvi
complex, dimension(1:imaxfreq) ::  ppmat,psvmat,pshmat,svpmat,svsvmat,svshmat
complex, dimension(1:imaxfreq) ::  shpmat,shsvmat,shshmat
real ::				df,dt
real, dimension(1:imaxfreq) ::	freq,pxtrace,pytrace,pztrace
real, dimension(1:imaxfreq) ::	svxtrace,svytrace,svztrace
real, dimension(1:imaxfreq) ::	shxtrace,shytrace,shztrace
integer	::			nfreq,ntime,i,iargc
integer ::			specfile,timefile,filtfile,sifile
character (len=256)::		infile, bandphase,save_rt_string
logical ::                      save_rt,print_orig_seis
#ifdef USE_COMMAND_LINE_AND_FILES
!	command line arguments
if(iargc().lt.2) then
   print*, 'usage: spectoseis spectrumfile bandphase [save_rT]'
   print*, '       calculates qP, qSV and qSH synthetics in stick (comb_out)'
   print*, '       and filtered (aniso_filt_out) format'
   print*, '       spectrumfile is output of anicake (cake_out)'
   print*, '       bandphase is SKS (pulse width 6s), P (band 10-100s), or SI for splitting intensity'
   print*, '       save_rT: if >0, will save radial and transverse component into file'
   stop
end if
call getarg(1,infile)
call getarg(2,bandphase)
if(iargc().gt.2)then
   call getarg(3,save_rt_string)
   read(save_rt_string,*)i
   if(i.gt.0)then 
      save_rt = .true. 
   else 
      save_rt = .false.
   endif
else
   save_rt = .false.
endif
!       input 
specfile = 10
open(specfile,file = infile)
!       time domain output
timefile = 11
open(timefile,file = 'comb_out')
!       filtered output
filtfile = 12
open(filtfile,file = 'aniso_filt_out')
print_orig_seis = .true.
#else
if(iargc().lt.1) then
   print*, 'usage: spectoseis_stream bandphase [save_rt]'
   print*, '       calculates qP, qSV and qSH synthetics in stick (comb_out)'
   print*, '       and filtered (aniso_filt_out) format'
   print*, '       output of anicake (cake_out) is read from stdin'
   print*, '       bandphase is SI (center 12.2s), SKS (pulse width 6s) or P (band 10-100s)'
   print*, '       save_rT: if >0, will save radial and transverse component into file'
   stop
end if
call getarg(1,bandphase)
if(iargc().gt.1)then
   call getarg(2,save_rt_string)
   read(save_rt_string,*)i
   if(i.gt.0)then 
      save_rt = .true. 
   else 
      save_rt = .false.
   endif
else
   save_rt = .false.
endif

!       defaults:
!       get input from stdin
specfile = FTRN_STDIN		! stdin
!       timefile will be omitted,
!       filtered time file is stdout
filtfile = FTRN_STDOUT

print_orig_seis = .false.

#endif
if(save_rt)then
   sifile=20		! splitting intensity
   open(sifile,file = 'sirt.out')
endif
!
!read input        
!
call spec_read_file(specfile,qp,qsv,qsh,qpi,qsvi,qshi,&
     ppmat,psvmat,pshmat,svpmat,svsvmat,svshmat,shpmat,shsvmat,shshmat,&
     imaxfreq,freq,nfreq)

!
! compute filtered seismograms
!

call  spec_compute_seis(bandphase,qp,qsv,qsh,qpi,qsvi,qshi,&
     ppmat,psvmat,pshmat,svpmat,svsvmat,svshmat,shpmat,shsvmat,shshmat,&
     pxtrace,pytrace,pztrace,svxtrace,svytrace,svztrace,shxtrace,shytrace,shztrace,&
     imaxfreq,freq,nfreq,ntime,dt,&
     timefile,print_orig_seis)
!
! print
!
call spec_print_seis(dt,imaxfreq,pxtrace,pytrace,pztrace,svxtrace,svytrace,svztrace,&
     shxtrace,shytrace,shztrace,ntime,filtfile,sifile,save_rt)

close(filtfile)
if(save_rt)then	
   close(sifile)
endif

end program main

