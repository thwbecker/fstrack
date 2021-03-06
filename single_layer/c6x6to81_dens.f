      program tensorformat
c
c    based on code anisotropy written by Dan Raymer, Leeds University.
c    portion to put Cij into Cijkl format
* 30.10.2000

* this version normalises by density (fixed to 3.353 g/cm^3 - olivine)
*
* input: Cij in GPa, not normalised by density; order 11,22,..,12,13,..,23,..
*
* either in old format from file, or new format from stdin
*
* output: cijkl in GPa/density[g/cm^3] for input in christoffel equation solver
*         tensor_phase
* From Vera in Feb 2005
c
c $Id: c6x6to81_dens.f,v 1.2 2005/07/01 21:54:36 becker Exp $
c

      implicit none 
      real*8  c(6,6),cc(3,3,3,3),rho,xlon,xlat,depth
      character*30  fname,fout
      integer ip,op,i,j,iv,jv,inum
      parameter (rho=3.353)
      inum = 21
#ifndef STREAM
      print*,'input Cij filename: '
      read(*,9) fname
      print*,'output Cijkl filename: '
      read(*,9) fout
 9    format(a30)
      open(7,file=fname)
      open(8,file=fout,status='unknown')
      ip = 7
      op = 8
#else
      ip = 5                    ! stdin
      op = 6
#endif
      do i=1,6
         do j=1,6
            c(i,j)=0d0
         enddo
      enddo
#ifndef STREAM
c     read in upper right triangle part
      do i=1,inum
         read(ip,*)iv,jv,c(iv,jv) 
      enddo
#else                           ! new format
      read(ip,*)xlon,xlat,depth,((c(iv,jv),jv=iv,6),iv=1,6)

#endif
c     fill in and normalize
      call vera_fill_and_normalize_c2(c,rho)
c     convert
      call vera_c2c4(c,cc)
c     print
      call vera_print_cijkl_ftrn(cc,op)
#ifndef STREAM
      close(ip)
      close(op)
#endif

      stop
      end
      

      
