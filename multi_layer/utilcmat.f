c
c aniosotropic reflectivity package from Vera Schulte-Pelkum in March 2005
c based on code by Peter Shearer and various others, see individual READMEs 
c and source code comments for details and description
c  
c minor modifications by Thorsten Becker
c
c $Id: utilcmat.f,v 1.2 2005/03/24 23:38:11 becker Exp $
c
c
cc subs is set of utility subroutines for manipulating
c complex matrices.
c
c
c subroutine CZEROOUT sets complex matrix to zero
      SUBROUTINE CZEROOUT(c,n)
      implicit double precision (a-h,o-z)
      double complex c(n,n)
      do 20 i=1,n
         do 10 j=1,n
            c(i,j)=cmplx(0.0d0,0.0d0)
10       continue
20    continue
      return
      end
c
c
c subroutine CMMULT multiplies two complex square matrices
      SUBROUTINE CMMULT(a,b,c,n)
      implicit double precision (a-h,o-z)
      double complex a(n,n),b(n,n),c(n,n)
      do 40 i=1,n
      do 30 j=1,n
         c(i,j)=(0.0d0,0.0d0)
         do 20 k=1,n
20          c(i,j)=c(i,j)+a(i,k)*b(k,j)
30    continue
40    continue
      return
      end
c
c
c subroutine CADD adds two complex square matrices
      SUBROUTINE CADD(a,b,c,n)
      implicit double precision (a-h,o-z)
      double complex a(n,n),b(n,n),c(n,n)
      do 40 i=1,n
      do 30 j=1,n
         c(i,j)=a(i,j)+b(i,j)
30    continue
40    continue
      return
      end
c
c
c subroutine CSUB subtracts two complex square matrices
      SUBROUTINE CSUB(a,b,c,n)
      implicit double precision (a-h,o-z)
      double complex a(n,n),b(n,n),c(n,n)
      do 40 i=1,n
      do 30 j=1,n
         c(i,j)=a(i,j)-b(i,j)
30    continue
40    continue
      return
      end
c
c
c subroutine CTRANS finds transpose of complex square matrix
      SUBROUTINE CTRANS(a,b,n)
      implicit double precision (a-h,o-z)
      double complex a(n,n),b(n,n)
      do 20 i=1,n
      do 10 j=1,n
         b(i,j)=a(j,i)
10    continue
20    continue
      return
      end
c
c
c subroutine CMULT multiplies complex square matrix by complex vector
      SUBROUTINE CMULT(a,b,c,n)
      implicit double precision (a-h,o-z)
      double complex a(n,n),b(n),c(n)
      do 40 i=1,n
         c(i)=cmplx(0.0d0,0.0d0)
         do 20 k=1,n
20          c(i)=c(i)+a(i,k)*b(k)
40    continue
      return
      end
c
c
c subroutine CSCALE multiplies complex square matrix by complex scalar
      SUBROUTINE CSCALE(a,scale,b,n)
      implicit double precision (a-h,o-z)
      double complex a(n,n),b(n,n),scale
      do 20 i=1,n
      do 10 j=1,n
         b(i,j)=a(i,j)*scale
10    continue
20    continue
      return
      end
c
c
c CDISPMAT displays a complex square matrix
      SUBROUTINE CDISPMAT(c,n)
      implicit double precision (a-h,o-z)
      double complex c(n,n)
      do 10 i=1,n
10       print 20,(c(i,j),j=1,n)
20    format (6(2f10.4,1x))
      return
      end
c
c
c CVECDISP displays a complex vector
      SUBROUTINE CVECDISP(c,n)
      implicit double precision (a-h,o-z)
      double complex c(n)
      do 10 i=1,n
10       print 20, c(i)
20    format (2f10.4)
      return
      end
c
c
c CINVERT finds inverse of complex square matrix (max size 10x10)
      SUBROUTINE CINVERT(c,cinv,n,cond_num)
      implicit double precision (a-h,o-z)
      double complex c(n,n),cinv(n,n),z(10),det(2)
      real*8 cond_num
      integer ipvt(10),job
      do 20 i=1,n
         do 10 j=1,n
         cinv(i,j)=c(i,j)
10       continue
20    continue
c                {slatec subroutines follow
      call CGECO(cinv,n,n,ipvt,cond_num,z)
      job=1
      call CGEDI(cinv,n,n,ipvt,det,z,job)
      return
      end
