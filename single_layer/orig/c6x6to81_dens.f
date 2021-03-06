      program tensorformat
c
c    based on code anisotropy written by Dan Raymer, Leeds University.
c    portion to put Cij into Cijkl format
* 30.10.2000

* this version normalises by density (fixed to 3.353 g/cm^3 - olivine)
* input: Cij in GPa, not normalised by density; order 11,22,..,12,13,..,23,..
* output: cijkl in GPa/density[g/cm^3] for input in christoffel equation solver
*         tensor_phase
* Feb 2005

      implicit real*8 (a-h,n-z)
      real*8  c(6,6),cc(3,3,3,3),ec(21), rho
      character*30  fname,fout
      parameter (pi=3.1415926535897932)

	print*,'input Cij filename: '
      read(*,9) fname
	print*,'output Cijkl filename: '
	read(*,9) fout
9	format(a30)
	inum = 21

      open(7,file=fname)
      open(8,file=fout,status='unknown')
      
      do i=1,6
        do j=1,6
           c(i,j)=0d0
        enddo
      enddo
 
      do i=1,inum
         read(7,*)iv,jv,ec(i) 
         c(iv,jv)=ec(i)/rho
         c(jv,iv)=c(iv,jv)
      enddo
      call c2c4(c,cc)


      write(8,10) ((((cc(i,j,k,l),
     +              i=1,3),j=1,3),k=1,3),l=1,3)
 10   format(1p,9e14.6)
      
      stop
      end
      
      SUBROUTINE c2c4(c,cc)   
      REAL*8 c(6,6),cc(3,3,3,3)
      INTEGER i,j,k,l 
c
         do 20 i=1,3
            do 15 j=1,3
               do 10 k=1,3
                  do 5 l=1,3
                     cc(i,j,k,l) = 0.d0
 5                continue
 10            continue
 15         continue
 20      continue
c
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
c
c        do 10 i=1,6
c                do 5 j=i+1,6
c                        c(j,i)=c(i,j)
c5               continue  
c10      continue
c
c
      return
      end

      
