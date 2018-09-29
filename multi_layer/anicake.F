c     
c     aniosotropic reflectivity package from Vera Schulte-Pelkum in March 2005
c     based on code by Peter Shearer and various others, see individual READMEs 
c     and source code comments for details and description
c     
c     this version can read full anisotropic tensors at up to 500 layers and compute
c     seismograms from them. produce_splitting_layers.awk can be used to generate an input file
c     depending on the idisp setting, spectral input for spectoseis, full seismograms, or splits can 
c     be computed. 
c     for input file examples, see the examples/ subdirectory and produce_splitting_layers.awk
c     
c     $Id: anicake.F,v 1.1 2010/12/31 19:00:09 becker Exp becker $
c     
c     
c     c Program ANICAKE is modified from ANIPLANE to permit each layer to be
c     defined separately.
c     
c     NOTE:  All variables (except input parameters) are complex!
c     Maximum number of layers = 501
c     
c     REQUIRES: ANIEIGEN, ANILAYER, ANIFACE, CINVERT, UTILCMAT, UTILCIJKL, GETHEX
c     
c     
c     minor modifications by Thorsten Becker (twb@usc.edu)
c     
c     - if output file is set to stdout, stdout is used to print
c     the harmonic response
c     - delaytime file output was switched off
c     - modified to allow looping over rotation azimuth     
c     - direct seismogram output (idisp == 3,4,5,6,7 for no, P, SKS, 
c     SI (will also write RT components), SK2 filters) add 10 or 20 for splitting computation (see below)
c     
c     the input file has a flag idisp at the end, which determines the operational mode:
c     
c     idisp = 1 regular, short output which can be read by spectoseis or spectoseis_stream
c     = 2 long output
c     = 3 unfiltered, straight seismogram output to output file
c     = 4 P   filtered seismogram output to output file
c     = 5 SKS filtered seismogram output to output file
c     = 6 SI filtered seismogram output to output file, also write radial and transverse traces of SV
c     = 7 SK2 filtered seismogram output 
c     = 8 SK3 filtered seismogram output 
c     = 15, 16, 17, 18, as above, but 
c     compute shear wave splitting with Levin & Park 1999's correlation method and write splitting output
c     = 25, 26, 27, 28, as above, but 
c     compute shear wave splitting with Menke's correlation method and write splitting output
c     examples: 
c     
c     - compute spectral response of layer stack, generate seismograms, filter to SKS band, 
c     and compute splitting with Levin & Park method:
c     
c     cat anicake.input | anicake | spectoseis_stream SKS | skssplit_xcorr_stream 0.0189 180
c     
c     where the last two parameters (0.0189 180) are just echoed. anicake input files are in examples/
c     and we've assumed all binaries are in your path. the above is equivalent to
c     
c     cat anicake.vsplit.input | ~/progs/src/fstrack/multi_layer/bin/i686/anicake
c     
c     where anicake.vsplit.input differs from anicake.input simply by a changed idisp, from 1 to 15
c     
c     - similarly, for the Menke method, use
c     
c     
c     cat anicake.input | anicake | spectoseis_stream SKS | ah_cross_conv_spectoseis  0.0189 180
c     
c     or 
c     
c     cat anicake.msplit.input | anicake
c     
c     - if you want to loop over backazimuths, use
c
c
c     cat anicake.vsplitrot.input | anicake
c
c
      program main

      implicit none

!     
!     max number of layers
!     
c      parameter max_layer=500
#define max_layer 500
!     
!     for spectoseis output
c      parameter imaxfreq=20000  ! max number of samples
!
#define imaxfreq 20000  

!     
!     
      character*80 vel_file,out_file,coup_file,suffix,filename
      integer laydep(5),outs
      logical debug,first
      real*8 c0(3,3,3,3),cr(3,3,3,3),c1(3,3,3,3),rho
      real*8 sym33(3,3),sym33r(3,3),acfln(5)
      real*8 cael(3,3,3,3),ciso(3,3,3,3),dcdz(3,3,3,3)
      real*8 pvel0,pgrad,dz,zmin,zmax,z,zmid(0:max_layer)
      real*8 rp,cond_num,rh(0:max_layer)
      double complex h(0:max_layer)
      real*8 freq,freqmin,freqmax,freqspace
      real*8 vtt(3)
c     
c     added for anicake to save info for each layer (max_layer layer maximum)
      real*8 clay(3,3,3,3,0:max_layer),rholay(0:max_layer),
     &     symlay(3,0:max_layer),clay_orig(3,3,3,3,0:max_layer)
c     
c     ANIEIGEN parameters
      real*8 sym(3)
      double complex q(3,6),pol(3,6,0:max_layer),
     &     e(6,6,0:max_layer),dscr(6,6)
      double complex d(6,6,0:max_layer),dinv(6,6,0:max_layer),
     &     flux(6,6,0:max_layer)
c     
c     ANIFACE parameters
      double complex tran_d(3,3,0:max_layer),tran_u(3,3,0:max_layer)
      double complex ref_d(3,3,0:max_layer),ref_u(3,3,0:max_layer)
      double complex lamda(6,6,0:max_layer),lam_d(3,3,max_layer),
     &     lam_u(3,3,0:max_layer)
      double complex vslow(6,0:max_layer),rt(6,6,0:max_layer)
c     
c     generalized reflection/transmission coefficients
      double complex gtran_d(3,3,0:max_layer),gtran_u(3,3,0:max_layer)
      double complex gref_d(3,3,0:max_layer),gref_u(3,3,0:max_layer)
      double complex sref_d(3,3),stran_d(3,3),stran_u(3,3),sref_u(3,3)
c     
c     coupling variables
      double complex dddz(6,6),couple(6,6)
c     
c     scratch and convenience variables
      double complex ci33(3,3),cs1
      double complex c2scr1(3,3),c2scr2(3,3),c2scr3(3,3),c2scr4(3,3)
      real*8 rscr1,rscr2,rscr3,theta,symt(3,3)

      double complex czero,cone,cpi,p,omega
      integer i,jj,j,k,icoefmat,num_dep,m,n,num_lay,idisp
      
c     
c     output straight to spectoseis
c     
      integer ifreqc,ntime,sifile,timefile,ifreqc_max
      logical save_rt,print_orig_seis
      character*256 bandphase
      real pxtraces(1:imaxfreq),pytraces(1:imaxfreq),
     &     pztraces(1:imaxfreq)
      real svxtraces(1:imaxfreq),svytraces(1:imaxfreq),
     &     svztraces(1:imaxfreq)
      real shxtraces(1:imaxfreq),shytraces(1:imaxfreq),
     &     shztraces(1:imaxfreq)

      real freqs(1:imaxfreq),dt
      complex ppmats(1:imaxfreq),svpmats(1:imaxfreq),
     &     shpmats(1:imaxfreq),
     &     psvmats(1:imaxfreq),svsvmats(1:imaxfreq),
     &     shsvmats(1:imaxfreq),
     &     pshmats(1:imaxfreq),svshmats(1:imaxfreq),
     &     shshmats(1:imaxfreq)
      complex qps(3),qpis(3),qsvs(3),qsvis(3),qshs(3),qshis(3)

c     
c     splitting computations
c     
      integer :: split_mode
      real :: time(1:imaxfreq),misfit,fastaz,delta
      integer :: null
c     
c     loop through azimuth, if azimuth is read in as -999
c     
      integer pass
      logical loop_azi,rot_tensor(0:max_layer)
      real*8 start_azi,stop_azi,dazi,azi_ref,azi_read
c     for back-azimuth loops

      do i=0,max_layer
         rot_tensor(i) = .false.
      enddo
      loop_azi = .true.
      start_azi = 180d0
      stop_azi = -178d0
      dazi = -2.0d0             ! will be changed to 10 for "SI"
      azi_ref = start_azi
c     
c     
c     init
c     

      czero=dcmplx(0.0d0,0.0d0)
      cone=dcmplx(1.0d0,0.0d0)
      cpi=dcmplx(3.1415926535897932384626433832795d0,0.0d0)

      call czeroout(ci33,3)
      do  i=1,3
         ci33(i,i)=cone         !ci33 is 3x3 identity matrix
      enddo
      
c     
c     main loop
c
c      print *,'OK start'
      pass=0
      do while(loop_azi.and.(azi_ref+dazi.ge.stop_azi)) ! main loop for azimuth scans

         pass=pass+1

         if(pass.eq.1)then
                                !
                                ! first initialization pass through routine
                                !
                                ! switch looping off, any azimuth == -999 will switch on
                                ! 
            loop_azi = .false.
                                !
                                ! read from anicake.input file, see examples/ subdirectory for 
                                ! input file examples
                                !
                                ! READ LAYER NUMBER
                                !
            read *,num_lay
            num_lay=num_lay-1   !shift to include layer # 0
            if (num_lay.gt.max_layer) then
               print *,'anicake: error: too many layers you have',
     &              num_lay
               stop
            end if
            z=0.0d0
            do k=0,num_lay      ! loop through layers
                                !
                                !
                                ! layer thickness, TOP TO BOTTOM
!     
               read *,dz
               h(k)=dcmplx(dz,0.0d0)
               rh(k)=dz
               zmid(k)=z        !used for coupling matrix, probably not needed now
               z=z+dz
!     
!     
!     read anisotropic  parameters for this layer
!     
               call vera_gethex(rholay(k),clay(1,1,1,1,k),
     &              clay_orig(1,1,1,1,k),acfln,sym33r,
     &              rot_tensor(k),start_azi,azi_read)
               if(.not.loop_azi)then
                  if(rot_tensor(k))then
                     loop_azi = .true.
                  endif
               endif
!               print *,k,rot_tensor(k),loop_azi
               if((k.eq.1).and.(.not.rot_tensor(k)))then
                  azi_ref = azi_read ! formally, use the read azi of the second layer 
               endif
               do  i=1,3
                  symlay(i,k)=sym33r(3,i)
               enddo
*     print *,'Symmetry axis =',(symlay(i,k),i=1,3)
c     call HOWWEAK(clay(1,1,1,1,k),ciso,cael,rscr1,rscr2,rscr3)
c     pvel0=dsqrt(ciso(1,1,1,1))
*     print *,'Isotropic P-wave velocity=',pvel0
               
c     print *,'rho=',rholay(k)
c     *     dump full tensor for test
c     print*, 'dumping tensor from main, line 80'
c     open(398,file='tensor.dump')
c     write(398,'(a)') '-----------------layer------------'
c     write(398,'(1p,9e14.6)')((((clay(idump,jdump,kdump,ldump,k),
c     + 	     idump=1,3),jdump=1,3),kdump=1,3),ldump=1,3)
               
            enddo               ! end layer loop

!     
                                !
                                ! RAY PARAMETER 
            read *,rp
                                !
            p=dcmplx(rp,0.0d0)

                                !
                                ! frequency range
                                !
            read *,freqmin,freqmax,freqspace

            if(((freqmax-freqmin)/freqspace+1) .ge.imaxfreq)then
               print *,'anicake: freqs:',freqmin,freqmax,freqspace
               print *,'anicake: increase imaxfreq'
               stop
            endif


*     print *,'Enter grt/c matrix desired:  (1) tran_u   or (2) ref_d'
            read *,icoefmat
*     print *,'Enter gr/tc output file name'
c     
c     output file (can be stdout)
c     
            read (*,('(a)')) out_file
            if(out_file.eq.'stdout')then
c     
C     output to stdout 
c     
               outs = 6
            else                ! output to file
               outs = 12
            endif
*     print *,'Enter number of output depths (maximum of 5)'
            read *,num_dep
            do  i=1,num_dep
*     print *,'Enter layer number for output layer #',i
               read *,laydep(i)
            enddo
!     
*     print *,'Enter  (1) short (2) long display, or (3) ONLY seismogram output'
!     
            read *,idisp
!     
!     
!     

            if(outs.eq.12)then ! write to file
               filename = out_file
               open (outs,file=filename)
               print *,'anicake: writing to ',trim(filename)
            endif
            if(idisp.lt.3)then

               write (outs,'(1x,a20)') filename
               write (outs,'(a)') 'spectrum'
               write (outs,*) '(1:gtran-u  2:gref_d) (nlay) (ndep)'
               write (outs,36) icoefmat,num_lay,num_dep
 36            format (3i4)
               write (outs,*) '(freqmin) (freqmax) (freqspace)'
               write (outs,37) freqmin,freqmax,freqspace
               write (outs,*) '(pgrad) (dz) (zmax)'
               write (outs,37) pgrad,dz,zmax
 37            format (3f10.4)
            else                ! only seismogram output
               if(num_dep.ne.1)then
                  print *,'anicake: error: idisp >=3 and num_dep != 1'
                  stop
               endif
               if(icoefmat.ne.1)then
                  print *,'anicake: error: idisp >=3 and icoefmat != 1'
                  stop
               endif
            endif
         else

!     
!     second pass
!     
!     
!     rotate the tensors further around the azimuth, axis 3
!     
!     
            if(outs.eq.12)then ! write to file
               filename = trim(out_file) //  suffix
               open (outs,file=filename)
               print *,'anicake: looping: writing to ',trim(filename)
            endif

            azi_ref = azi_ref + dazi
c
c     only rotate a selection of tensors 
c            
            do k=0,num_lay
               if(rot_tensor(k))then
!     rotate original tensor
                  call vera_rotate(clay_orig(1,1,1,1,k),clay(1,1,1,1,k),
     &                 3,azi_ref)
!     symmetry matrix
                  call vera_get_isym(symt,1)
                                ! rotate symmetry matrix
                  call vera_sym_rot(symt,sym33r,3,azi_ref)
                  do  i=1,3
                     symlay(i,k)=sym33r(3,i)
                  enddo
               endif
            enddo
                                ! end the second pass loop
         endif

c     
c     Get e and vslow for each layer
         vtt(1)=0.0d0
         vtt(2)=0.0d0
         vtt(3)=0.0d0
         do  i=0,num_lay        ! loop through layers


            call vera_anieigen(rholay(i),clay(1,1,1,1,i),p,symlay(1,i),
     &           q,pol(1,1,i),e(1,1,i),dscr,
     &           d(1,1,i),dinv(1,1,i),flux(1,1,i))
            

            do  j=1,6
               vslow(j,i)=q(3,j)
            enddo
c     
c     Remove eigenvector sign ambiguity by comparing with previous layer.
*     if (i.gt.1) then
            if (i.gt.0) then
               do j=1,6         !loop over eigenvectors
                  cs1=czero
                  do  m=1,6
                     cs1=cs1+e(m,j,i)*e(m,j,i-1)
                  enddo
                  if (dreal(cs1).lt.0.0d0) then
                     do  m=1,6
                        e(m,j,i)=-e(m,j,i)
                        d(m,j,i)=-d(m,j,i)
                        dinv(j,m,i)=-dinv(j,m,i)
                        if (m.le.3) pol(m,j,i)=-pol(m,j,i)
                     enddo
                  end if
               enddo
            end if
c     
            if (idisp.eq.2) then ! long output
*     print *,'LAYER',i,' polarization, vslow follow'
               print 43, (pol(jj,1,i),jj=1,3),vslow(1,i)
               print 43, (pol(jj,2,i),jj=1,3),vslow(2,i)
               print 43, (pol(jj,3,i),jj=1,3),vslow(3,i)
               print 43, (pol(jj,4,i),jj=1,3),vslow(4,i)
               print 43, (pol(jj,5,i),jj=1,3),vslow(5,i)
               print 43, (pol(jj,6,i),jj=1,3),vslow(6,i)
            end if
            vtt(1)=vtt(1)+rh(i)*abs(dreal(vslow(4,i)))
            vtt(2)=vtt(2)+rh(i)*abs(dreal(vslow(5,i)))
            vtt(3)=vtt(3)+rh(i)*abs(dreal(vslow(6,i)))
            
         enddo                  ! end layer loop
                                !
                                ! those three are needed for spectoseis, idisp>=3
                                !
         do jj=1,3              
            qps(jj) = pol(jj,4,0)
            qsvs(jj) = pol(jj,5,0)
            qshs(jj) = pol(jj,6,0)
            qpis(jj) = pol(jj,4,num_lay)
            qsvis(jj) = pol(jj,5,num_lay)
            qshis(jj) = pol(jj,6,num_lay)
         enddo
         
*     print *,'eigenvector matrix calculated for all layers'
         if(idisp.lt.3)then
                                !
                                ! regular output
                                ! 
            write (outs,*) ' '
            write (outs,*) 
     &           'Topmost layer complex eigenvectors, vslow follow'
            write (outs,43) (pol(jj,1,0),jj=1,3),vslow(1,0)
            write (outs,43) (pol(jj,2,0),jj=1,3),vslow(2,0)
            write (outs,43) (pol(jj,3,0),jj=1,3),vslow(3,0)
            write (outs,43) (pol(jj,4,0),jj=1,3),vslow(4,0)
            write (outs,43) (pol(jj,5,0),jj=1,3),vslow(5,0)
            write (outs,43) (pol(jj,6,0),jj=1,3),vslow(6,0)
            write (outs,*) 
     &           'Next down layer complex eigenvectors, vslow follow'
            write (outs,43) (pol(jj,1,1),jj=1,3),vslow(1,1)
            write (outs,43) (pol(jj,2,1),jj=1,3),vslow(2,1)
            write (outs,43) (pol(jj,3,1),jj=1,3),vslow(3,1)
            write (outs,43) (pol(jj,4,1),jj=1,3),vslow(4,1)
            write (outs,43) (pol(jj,5,1),jj=1,3),vslow(5,1)
            write (outs,43) (pol(jj,6,1),jj=1,3),vslow(6,1)
            write (outs,*) 
     &           'Bottom layer complex eigenvectors, vslow follow'
            write (outs,43) (pol(jj,1,num_lay),jj=1,3),
     &           vslow(1,num_lay)
            write (outs,43) (pol(jj,2,num_lay),jj=1,3),
     &           vslow(2,num_lay)
            write (outs,43) (pol(jj,3,num_lay),jj=1,3),
     &           vslow(3,num_lay)
            write (outs,43) (pol(jj,4,num_lay),jj=1,3),
     &           vslow(4,num_lay)
            write (outs,43) (pol(jj,5,num_lay),jj=1,3),
     &           vslow(5,num_lay)
            write (outs,43) (pol(jj,6,num_lay),jj=1,3),
     &           vslow(6,num_lay)

 43         format (4(2e12.5,2x))
            write (outs,'(a30,3(e12.5,1x))')
     &           'Upgoing tt through stack=',
     &           vtt(1),vtt(2),vtt(3)
c     
*     c      print*, 'shear wave splitting delay time ', vtt(3)-vtt(2),' s'
c     open(95,file='delaytime.out')
c     write(95,*) vtt(3)-vtt(2)
c     close(95)
         endif                  ! end idisp<3 branch
         
c     Get interface reflection/transmission coefficients
         do i=1,num_lay

c     ***NOTE:  I changed inputs to ANIFACE to e matrices from d matrices
c     to get this to work properly.  I don't know how it could have
c     been working in aniplane!
            

            
            call vera_aniface(e(:,:,i-1),e(:,:,i),p,
     &           rt(:,:,i),tran_d(:,:,i),tran_u(:,:,i),
     &           ref_d(:,:,i),ref_u(:,:,i))

         enddo

         if(idisp.lt.3)then
*     print *,'rt/c calculated for all layers'
            
*     print *,'FREE SURFACE ref_u follows'
c     call CDISPMAT(ref_u(1,1,1),3) 
*     print *,'FREE SURFACE tran_d follows'
c     call CDISPMAT(tran_d(1,1,1),3)
            write (outs,*) 'free surface ref_u follows'
            write (outs,311) (ref_u(1,i,1),i=1,3)
            write (outs,311) (ref_u(2,i,1),i=1,3)
            write (outs,311) (ref_u(3,i,1),i=1,3)
            write (outs,*) ' '
         endif
         
c     Start loop over frequency
         first=.true.
c         do  freq=freqmin,freqmax,freqspace
         ifreqc_max = (freqmax-freqmin)/freqspace
         freq=freqmin
c     print *,ifreqc_max,freqmax,freqmin,freqspace
         do ifreqc=1,ifreqc_max
!
                                ! single prec frequency
                                !
            freqs(ifreqc) = freq
!
            if (freq.ne.0.0d0) then
               omega=dcmplx(freq*2.0d0*3.1415927d0,0.0d0)
            else
               omega=dcmplx(0.001d0*2.0d0*3.1415927d0,0.0d0) !otherwise will bomb
            end if
            
c     Get all the layer matrices for each layer
            do  i=1,num_lay
               call vera_anilayer(vslow(1,i),h(i),omega,
     &              lamda(1,1,i),lam_d(1,1,i),lam_u(1,1,i))
            enddo
            
c     Now calculate generalized reflection/transmission coefficients
c     for the entire stack following Kennett section 6.2.1
c     Note that final generalized matrices do not include effect of
c     the thickness of layer 1.
            do i=1,3
               do j=1,3
                  gtran_d(i,j,num_lay)=tran_d(i,j,num_lay) !bottom interface
                  gtran_u(i,j,num_lay)=tran_u(i,j,num_lay) !is simply RT/C
                  gref_d(i,j,num_lay)=ref_d(i,j,num_lay) !for interface
                  gref_u(i,j,num_lay)=ref_u(i,j,num_lay)
               enddo
            enddo
            
*     open (13,file='xydump') !******
            
*     print *,'freq = ',freq
            do  n=num_lay-1,1,-1 !**last interface is 1 (include free-surface)
               
c     
c     lam_u used instead of lam_d since lam_u seems to match superplane.f lam_d
               call CMMULT(gref_d(1,1,n+1),lam_u(1,1,n),c2scr1,3) !Kennett
               call CMMULT(lam_u(1,1,n),c2scr1,sref_d,3) !Equa. 6.16
               call CMMULT(gtran_d(1,1,n+1),lam_u(1,1,n),stran_d,3)
               call CMMULT(lam_u(:,:,n),gtran_u(:,:,n+1),stran_u,3)
               
               call CMMULT(ref_u(1,1,n),sref_d,c2scr1,3) !Equa. 6.17
               call CSUB(ci33,c2scr1,c2scr2,3)
               call CINVERT(c2scr2,c2scr1,3,cond_num)
               call CMMULT(c2scr1,tran_d(1,1,n),c2scr2,3)
               call CMMULT(sref_d,c2scr2,c2scr1,3)
               call CMMULT(tran_u(1,1,n),c2scr1,c2scr3,3)
               call CADD(ref_d(1,1,n),c2scr3,gref_d(1,1,n),3)
               call CMMULT(stran_d,c2scr2,gtran_d(1,1,n),3)
               
               call CMMULT(sref_d,ref_u(1,1,n),c2scr1,3) !Equa. 6.18
               call CSUB(ci33,c2scr1,c2scr2,3)
               call CINVERT(c2scr2,c2scr1,3,cond_num)
               call CMMULT(c2scr1,stran_u,c2scr2,3)
               call CMMULT(ref_u(1,1,n),c2scr2,c2scr1,3)
               call CMMULT(stran_d,c2scr1,c2scr3,3)
               call CADD(gref_u(1,1,n+1),c2scr3,gref_u(1,1,n),3)
               call CMMULT(tran_u(:,:,n),c2scr2,gtran_u(:,:,n),3)
               
            enddo
            
c     ##############################################################################

            if(idisp.ge.3)then
!
! for spectoseis
!
               k=1
               ppmats(ifreqc) = gtran_u(1,1,laydep(k))
               svpmats(ifreqc) = gtran_u(1,2,laydep(k))
               shpmats(ifreqc) = gtran_u(1,3,laydep(k))
               psvmats(ifreqc) = gtran_u(2,1,laydep(k))
               svsvmats(ifreqc) = gtran_u(2,2,laydep(k))
               shsvmats(ifreqc) = gtran_u(2,3,laydep(k))
               pshmats(ifreqc) = gtran_u(3,1,laydep(k))
               svshmats(ifreqc) = gtran_u(3,2,laydep(k))
               shshmats(ifreqc) = gtran_u(3,3,laydep(k))
            else
!
! regular output
!     
               do k=1,num_dep   !loop over depths
                  if (icoefmat.eq.1) then
*     write (outs,310) freq
                     write (outs,'(a5,f26.16)') 'freq=',freq
                     write (outs,311) ((gtran_u(j,i,laydep(k)),
     &                    i=1,3),j=1,3)
*     write (13,*) freq,real(gtran_u(1,1,laydep(k))),
*     &              aimag(gtran_u(1,1,laydep(k))),
*     &              cabs(gtran_u(1,1,laydep(k)))
                     
                     
                     
                  else
                     write (outs,310) freq
                     write (outs,311) ((gref_d(j,i,laydep(k)),
     &                    i=1,3),j=1,3)
*     write (13,*) freq,real(gref_d(1,1,laydep(k))),
*     &              aimag(gref_d(1,1,laydep(k))),
*     &              cabs(gref_d(1,1,laydep(k)))
                  end if
               enddo
               
 310           format (f10.4)
 311           format (2e12.4,2x,2e12.4,2x,2e12.4)
            endif

            freq = freq + freqspace
         enddo                  ! end loop on frequencies
c         print *,freq,freqmax,ifreqc,imaxfreq
         if(idisp.ge.3)then
                                !
                                ! SEISMOGRAM OR SPLIT OUTPUT
                                !
            timefile=7          ! dummy for unfiltered, else will also be stdout
            print_orig_seis = .false.

            if(pass.eq.1)then
                                !
                                ! decide on the splitting mode
                                !
               save_rt=.false.
               if(idisp.gt.20)then
                  split_mode=2  ! menke
                  idisp=idisp-20
               else if(idisp.gt.10)then
                  split_mode=1  ! vera
                  idisp=idisp-10
               else
                  split_mode=0
               endif
               if(idisp.eq.3)then ! original seismogram to stdout
                  print_orig_seis = .true.
                  bandphase = ''
                  timefile=6
               else if(idisp.eq.4)then ! P wave filtered to stdout
                  bandphase = 'P'
               else if(idisp.eq.5)then ! SKS filtered to stdout
                  bandphase = 'SKS'
               else if(idisp.eq.6)then ! SI to stdout and also print out files for splitting intensity
                  bandphase = 'SI'
                  save_rt=.true.
                  dazi = -10d0
               else if(idisp.eq.7)then ! SKS, second filter to stdout
                  bandphase = 'SK2'
               else if(idisp.eq.8)then ! SKS, third filter to stdout
                  bandphase = 'SK3'
               else
                  print *,'anicake: error: idisp ',idisp,
     &                 'undefined'
                  stop
               endif
            endif
c            print *,'input OK'
!     
!     compute seismogram in time domain
!     

            call  spec_compute_seis(bandphase,qps,qsvs,qshs,qpis,qsvis,
     &           qshis,ppmats,psvmats,pshmats,svpmats,svsvmats,svshmats,
     &           shpmats,shsvmats,shshmats,
     &           pxtraces,pytraces,pztraces,svxtraces,svytraces,
     &           svztraces,shxtraces,shytraces,shztraces,imaxfreq,freqs,
     &           ifreqc,ntime,dt,timefile,print_orig_seis)
            if(save_rt)then
!
!     print rt components to file
!
               if(pass.eq.1)then
                  filename = 'sirt.out'
               else
                  filename = 'sirt.out' // trim(suffix)
               endif
c               print *,'writing rt to ',trim(filename)
               open(sifile,file = filename)
               call spec_print_rt_seis(dt,ntime,imaxfreq,
     &              svxtraces,svytraces,svztraces,sifile)
               close(sifile)
            endif
            if(split_mode.eq.0)then
!     
!     simple seismogram output
!     
               call spec_print_seis(dt,imaxfreq,pxtraces,pytraces,
     &              pztraces,svxtraces,svytraces,svztraces,
     &              shxtraces,shytraces,shztraces,ntime,outs)
            else                ! splitting branch
               do i = 1,ntime
                  time(i) = (i-1)*dt
               enddo
               if(split_mode.eq.1)then 
!     
!     Levin & Park splitting
!     
                  call vera_split(time,svxtraces,svytraces,
     &                 ntime,fastaz,delta,misfit,null)
               else if(split_mode.eq.2)then 
!     
!     menke splitting
!     
                  call ah_cross_conv_spectoseis(time,svxtraces,
     &                 svytraces,ntime,fastaz,delta,misfit)
                  if(misfit.gt.0.5)then
                     null=1
                  else
                     null=0
                  endif
               endif
!
!     output
!     

               if(fastaz.lt.0)then ! make fast azimuth be within 0 ... 180 
                  fastaz = fastaz + 180d0
               endif
               if(fastaz.gt.180)fastaz=fastaz-180
               if (null.eq.1) then 
                  write(outs,'(f7.4,1x,f6.1,1x,a7,1x,a5,1x,e14.7)')
     &                 rp, azi_ref,' nan ',' 0.0 ',misfit
              else
                 write(outs,'(f7.4,1x,f6.1,1x,f6.1,1x,f6.2,1x,e14.7)')
     &                 rp,azi_ref,fastaz,delta,misfit
               endif
            endif
         endif
c     
         if(outs.ne.6)then
            close (outs)
         endif
         
                                ! string to append to filenames
         write(suffix, '(i3.3)') pass
         suffix='.'//trim(suffix)
      enddo                     ! end global iteration loop
      
c     stop
      end
