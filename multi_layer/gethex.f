c
c aniosotropic reflectivity package from Vera Schulte-Pelkum in March 2005
c based on code by Peter Shearer and various others, see individual READMEs 
c and source code comments for details and description
c  
c minor modifications by Thorsten Becker
c
c $Id: gethex.f,v 1.7 2011/04/12 06:19:51 becker Exp $
c
c
cc subroutine GETHEX prompts and obtains an anisotropic tensor
c
c can obtain a  hexagonally symmetric
c anisotropy model.
c
c or a full tensor (see below)
c
c
c Requires:  vera_zeroout,vera_fillout,vera_rotate,VERA_SYM_ROT,HUDSON,CAL_TENSOR
* changed to add read in of arbitrary full tensor c(i,j,k,l)
* 30.9.2000

!
! additional input by twb:
!
! rot_tensor: will be returned true if azimuth is input as -999
! azi_start,azi_step: will be used if loop_azi is set to true
! corig: original tensor before azimuth rotation
!
!
!
!     
      subroutine vera_gethex(rho,c,corig,
     &     acfln,sym,rot_tensor,azi_start,azi)
      
      implicit double precision (a-h,o-z)
      real*8 c(3,3,3,3),corig(3,3,3,3),
     &     cr(3,3,3,3),acfln(5),sym(3,3),symr(3,3),azi_ref
      logical rot_tensor,first_call,verbose
      data first_call /.true./
      save first_call
      character*256 tensorfile
      verbose = .false.
      
      rot_tensor = .false.

      itensorfile = 11
      if(verbose)then
         print *,' '
         print *,'1)  Keith & Crampin upper mantle'
         print *,'2)  Stephen JGR 1985'
         print *,'3)  Any hexagonal model'
         print *,'4)  Simple % anisotropy hexagonal (poisson solid)'
         print *,'5)  Backus coefficients'
         print *,'6)  Hudson crack model'
         print *,'7)  Standard crack model'
         print *,'8)  Simple % anisotropy hexagonal (non-poisson solid)'
         print *,'9)  Input full elastic tensor c(i,j,k,l)'
         print *,' '
      endif
10    continue
      if(verbose)then
         print *,'Enter desired model: '
      endif
      read *,menu
      if (menu.lt.1.or.menu.gt.9) go to 10
      call vera_zeroout(c)
      if (menu.eq.1) then
*        print *,'Keith & Crampin upper mantle model'
         rho=3.3d0
         c(1,1,1,1)=260.78d0/rho
         c(1,1,2,2)=80.0d0/rho
         c(1,2,1,2)=72.9d0/rho
         c(2,2,2,2)=200.77d0/rho
         c(2,2,3,3)=72.99d0/rho
         c(2,3,2,3)=63.89d0/rho
         c(3,3,3,3)=200.77d0/rho
         c(3,3,1,1)=80.0d0/rho
         c(1,3,1,3)=72.9d0/rho
         isym=1
      else if (menu.eq.2) then
*        print *,'Stephen JGR 1985 model'
         c(1,1,1,1)=25.0d0
         c(2,2,2,2)=16.0d0
         c(1,2,1,2)=8.333d0
         c(1,1,2,2)=8.333d0
         c(1,3,1,3)=5.336d0
         c(3,3,3,3)=c(1,1,1,1)               !additional dependent terms for
         c(3,3,2,2)=c(1,1,2,2)               !hexagonal symmetry, see Musgrave,
         c(3,2,3,2)=c(1,2,1,2)               !Crystal Acoustics, for details
         c(1,1,3,3)=c(1,1,1,1)-2.0d0*c(1,3,1,3)
         rho=1.96d0
         isym=2
      else if (menu.eq.3) then
*        print *,'General hexagonal model'
*         print *,'Enter c1111,c2222,c1212,c1122,c1313'
         read *,c(1,1,1,1),c(2,2,2,2),c(1,2,1,2),c(1,1,2,2),c(1,3,1,3)
         c(3,3,3,3)=c(1,1,1,1)               !additional dependent terms for
         c(3,3,2,2)=c(1,1,2,2)               !hexagonal symmetry, see Musgrave,
         c(3,2,3,2)=c(1,2,1,2)               !Crystal Acoustics, for details
         c(1,1,3,3)=c(1,1,1,1)-2.0d0*c(1,3,1,3)
         print *,'Enter density'
         read *,rho
         isym=2
      else if (menu.eq.4) then
         rho=1.                  !assume normalized density
*        print *,'Simple % anisotropy hexagonal model'
*         print *,'Enter average Vp, % anisotropy (e.g. 4, 10)'
         read *,vpavg,fract
         fract=fract/100.0d0
*         print *,'P-wave 4-theta term will be set'
*         print *,'to 2-theta term times factor'
*         print *,'Enter factor for 4-theta term'
         read *,ccfact
*         print *,'(100) symmetry axis is:  (1) slow  or  (2) fast?'
         read *, isym
         vpmax=vpavg+(vpavg*fract)/2.0d0
         vpmin=vpavg-(vpavg*fract)/2.0d0
         if (isym.eq.1) ccfact=-ccfact
         scr=(1.0d0+ccfact)/(1.0d0-ccfact)
         a=(vpmax**2+vpmin**2*scr)/(1.0d0+scr)
         b=(vpmax**2-a)/(1.0d0+ccfact)
         if (isym.eq.1) ccfact=-ccfact
         cc=b*ccfact
         vsavg=vpavg/sqrt(3.)    !assume 0.25 Poisson's ratio
         vsmax=vsavg+(vsavg*fract)/2.0d0
         vsmin=vsavg-(vsavg*fract)/2.0d0
         d=0.5d0*(vsmax**2+vsmin**2)
         e=vsmax**2-d
         if (isym.eq.1) then
            b=-b
            cc=-cc
            e=-e
         end if
*c         print*, vpmax,vpmin,vsmax,vsmin
*         print *,'a,b,c,d,e=',a,b,cc,d,e
         call CAL_TENSOR(c,a,b,cc,d,e)
         isym=1
*         print *,'Enter density'         !****added for anianza
         read *,rho
      else if (menu.eq.8) then
         rho=1.                  !assume normalized density
*        print *,'Simple % anisotropy hexagonal model'
*         print *,'Enter average Vp, Vs,% anisotropy (e.g. 4, 10)'
         read *,vpavg,vsavg,fract
         fract=fract/100.0d0
*         print *,'P-wave 4-theta term will be set'
*         print *,'to 2-theta term times factor'
*         print *,'Enter factor for 4-theta term'
         read *,ccfact
*         print *,'(100) symmetry axis is:  (1) slow  or  (2) fast?'
         read *, isym
         vpmax=vpavg+(vpavg*fract)/2.0d0
         vpmin=vpavg-(vpavg*fract)/2.0d0
         if (isym.eq.1) ccfact=-ccfact
         scr=(1.0d0+ccfact)/(1.0d0-ccfact)
         a=(vpmax**2+vpmin**2*scr)/(1.0d0+scr)
         b=(vpmax**2-a)/(1.0d0+ccfact)
         if (isym.eq.1) ccfact=-ccfact
         cc=b*ccfact
*         vsavg=vpavg/sqrt(3.0d0)    !assume 0.25 Poisson's ratio
         vsmax=vsavg+(vsavg*fract)/2.0d0
         vsmin=vsavg-(vsavg*fract)/2.0d0
         d=0.5d0*(vsmax**2+vsmin**2)
         e=vsmax**2-d
         if (isym.eq.1) then
            b=-b
            cc=-cc
            e=-e
         end if
*c         print*, vpmax,vpmin,vsmax,vsmin
*         print *,'a,b,c,d,e=',a,b,cc,d,e
         call CAL_TENSOR(c,a,b,cc,d,e)
         isym=1
*         print *,'Enter density'         !****added for anianza
         read *,rho
      else if (menu.eq.5) then
*        print *,'Enter P terms:  A, B (2-theta), C (4-theta)'
         read *,a,b,cc
*         print *,'Enter SV terms:  D, E (2-theta)'
         read *,d,e
         call CAL_TENSOR(c,a,b,cc,d,e)
         isym=1
      else if (menu.eq.6) then
*        print *,'Enter pvel,svel,rho for isotropic matrix'
         read *,vp,vs,rho
         if (vs.lt.0.0d0) then
            poi=-vs
            if (poi.ne.0.5d0) then
               vs=vp/sqrt((2.0d0-2.0d0*poi)/(1.0d0-2.0d0*poi))
            else
               vs=0.0d0
            end if
         end if
*         print *,'Enter pvel,svel,rho for crack material'
         read *,vp1,vs1,rho1
*         print *,'Enter crack aspect ratio d'
         read *,d
*         print *,'Enter crack density parameter epsilon'
         read *,epsilon
         call HUDSON(vp,vs,rho,vp1,vs1,rho1,d,epsilon,c,vpiso,vsiso)
         isym=1
      else if (menu.eq.7) then
         vp=4.5d0
         poi=.27d0
         vs=vp/sqrt((2.0d0-2.0d0*poi)/(1.0d0-2.0d0*poi))
         rho=2.8d0
         vp1=1.5d0
         vs1=0.0
         rho1=1.0d0
*         print *,'Enter crack aspect ratio d'
         read *,d
*         print *,'Enter crack density parameter epsilon'
         read *,epsilon
         call HUDSON(vp,vs,rho,vp1,vs1,rho1,d,epsilon,c,vpiso,vsiso)
         isym=1
      else if (menu.eq.9) then
!
! full tensor read
!
!
         if(verbose)then
            print*, 'File to read tensor from: '
         endif
         read(*,'(a)') tensorfile
      	 open(itensorfile,file=tensorfile)
         call vera_read_cijkl(c,itensorfile)
         close(itensorfile)
         if(verbose)then
            print*, 'factor to multiply coefficients by'
            print*, '(to adjust for units, density:'
            print*, 'need coeffs in GPa, density normalized)'
         endif
     	 read*, cijklfact
     	 if(cijklfact.ne.1.) then
     	   do i = 1,3
     	     do j = 1,3
     	       do k = 1,3
     	         do l = 1,3
     	           c(i,j,k,l) = c(i,j,k,l)*cijklfact
     	         end do
     	       end do
     	     end do
     	   end do
         endif
         if(verbose)then
            print*, 'density'
         endif
     	 read*, rho
     	 isym = 1
      end if
      if(menu.lt.9) call vera_fillout(c)

*     dump full tensor and hardwire menu back to 8 for test
*     (result: menu not used from here on)
c      print*, 'dumping tensor from gethex, line 189'
c      open(398,file='tensor.dump')
c      write(398,'(a)') '-----------------layer------------'
c      write(398,'(1p,9e14.6)')((((c(i,j,k,l),
c     +              i=1,3),j=1,3),k=1,3),l=1,3)
c      print*, 'density ', rho
c      print*, ' scale factor', cijklfact


c
c Find 5 independent constants for transverse isotropy and set up sym(3,3)

      call vera_get_isym(sym,isym)

      acfln(1)=c(1,1,1,1) 
      acfln(2)=c(3,3,3,3)
      acfln(3)=c(1,1,3,3)
      acfln(4)=c(2,3,2,3)
      acfln(5)=c(1,2,1,2)
c
      iaxis=1
      theta=0.0d0
c
c
c     rotations of tensors
c
*      print *,'Coordinates can be rotated clockwise'
*      print *,'Enter desired rotation axis, angle'
      read *,iaxis,theta
      call VERA_ROTATE(c,cr,iaxis,theta)
      call VERA_SYM_ROT(sym,symr,iaxis,theta)
      if(iaxis.eq.3)then
         print *,'gethex:  expecting ax 3 only for second rotation'
         stop
      endif
      
!     
!     azimuth part
!     
!     
*     print *,'Coordinates can be rotated clockwise'
*     print *,'Enter desired rotation axis, angle'
      read *,iaxis,theta
      if(iaxis.eq.3)then
         if(theta.eq.-999)then
!     
!     angle of -999 switches on looping
!     
            theta = azi_start
            rot_tensor = .true.
         endif
         azi = theta
      endif
!
! save original tensor before last, azimuthal rotation
!
      call vera_cijkl_assign(cr,corig)
! rotate
      call VERA_ROTATE(cr,c,iaxis,theta)
      call VERA_SYM_ROT(symr,sym,iaxis,theta)

      return
      end
