!  nrmod.f90
!  Independent "Numerical Recipes" methods
!
!  ARRAY_COPY     (subroutine)     : array copy
!  NRERROR        (subroutine)     : message error output
!  ARTH_I         (function)       : integer arth
!  ARTH_R         (function)       : simple precision arth
!  ARTH_D         (function)       : double precision arth
!  REALLOCATE_IV  (function)       : reallocate integer vector
!  REALLOCATE_IM  (function)       : reallocate integer matrix
!  OUTERPROD      (function)       : spreading product
!  OUTERDIFF      (function)       : spreading difference
!  GET_DIAG       (function)       : get diagonal coefficients
!  UNIT_MATRIX    (subroutine)     : form identity matrix
!  UPPER_TRIANGLE (function)       : form upper triangle logical mask
!  JACOBI         (subroutine)     : Jacobi diagonalization method
!  JROTATE        (int.subroutine) : Jacobi rotation
!  IMAXLOC        (function)       : maximum location of array
!  SWAP_D         (subroutine)     : swapping of dp values
!  SWAP_DV        (subroutine)     : swapping of dp arrays
!  EIGSRT         (subroutine)     : sorts eigenvalues & vectors
!  LUDCMP         (subroutine)     : LU decomposition
!  LUBKSB         (subroutine)     : LU back substitution
!  LUINV          (subroutine)     : LU inversion
!  GAULEG         (subroutine)     : Gaussian quadrature
!  LOCATE         (function)       : location value in array
!  SPLINT         (function)       : spline interpolation
!  TRIDAG_SER     (subroutine)     : tridiagonal serial algorithm
!  TRIDAG_PAR     (rec.subroutine) : tridiagonal parallel algorithm
!  SPLINE         (subroutine)     : spline second derivatives

      module NRMOD

      implicit none

      interface ARTH
       module procedure ARTH_I,ARTH_R,ARTH_D
      end interface

      interface REALLOCATE
       module procedure REALLOCATE_IV,REALLOCATE_IM
      end interface

      contains

      subroutine ARRAY_COPY(src,dest,n_copied,n_not_copied)
       !Copy of array "src" into "dest"
       real, dimension(:), intent(in) :: src
       real, dimension(:), intent(out) :: dest
       integer, intent(out) :: n_copied,n_not_copied
       n_copied=min(size(src),size(dest))
       n_not_copied=size(src)-n_copied
       dest(1:n_copied)=src(1:n_copied)
      end subroutine ARRAY_COPY

      subroutine NRERROR(string)
       !Message error output
       character(len=*), intent(in) :: string
       write (*,*) 'nrerror: ',string
       STOP 'program terminated by nrerror'
      end subroutine NRERROR

      function ARTH_I(first,increment,n)
       integer :: first,increment,n
       integer, dimension(n) :: arth_i
       integer :: k,k2,temp
       if (n > 0) arth_i(1)=first
       if (n<=16) then
        do k=2,n
         arth_i(k)=arth_i(k-1)+increment
        end do
       else
        do k=2,8
         arth_i(k)=arth_i(k-1)+increment
        end do
        temp=increment*8
        k=8
        do
         if (k >= n) exit
         k2=k+k
         arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
         temp=temp+temp
         k=k2
        end do
       end if
      end function ARTH_I

      function ARTH_R(first,increment,n)
       real, intent(in) :: first,increment
       integer, intent(in) :: n
       real, dimension(n) :: arth_r
       integer :: k,k2
       real :: temp
       if (n>0) arth_r(1)=first
       if (n<=16) then
        do k=2,n
         arth_r(k)=arth_r(k-1)+increment
        end do
       else
        do k=2,8
         arth_r(k)=arth_r(k-1)+increment
        end do
        temp=increment*8
        k=8
        do
         if (k>=n) exit
         k2=k+k
         arth_r(k+1:min(k2,n))=temp+arth_r(1:min(k,n-k))
         temp=temp+temp
         k=k2
        end do
       end if
      end function ARTH_R

      function ARTH_D(first,increment,n)
       double precision, intent(in) :: first,increment
       integer, intent(in) :: n
       double precision, dimension(n) :: arth_d
       integer :: k,k2
       double precision :: temp
       if (n>0) arth_d(1)=first
       if (n<=16) then
        do k=2,n
         arth_d(k)=arth_d(k-1)+increment
        end do
       else
        do k=2,8
         arth_d(k)=arth_d(k-1)+increment
        end do
        temp=increment*8
        k=8
        do
         if (k>=n) exit
         k2=k+k
         arth_d(k+1:min(k2,n))=temp+arth_d(1:min(k,n-k))
         temp=temp+temp
         k=k2
        end do
       endif
      end function ARTH_D

      function REALLOCATE_IV(p,n)
       integer, dimension(:), pointer :: p,reallocate_iv
       integer, intent(in) :: n
       integer :: nold,ierr
       allocate(reallocate_iv(n),stat=ierr)
       if (ierr/=0) call &
        NRERROR('reallocate_iv: problem in attempt to allocate memory')
       if (.not.associated(p)) RETURN
       nold=size(p)
       reallocate_iv(1:min(nold,n))=p(1:min(nold,n))
       deallocate(p)
      end function REALLOCATE_IV

      function REALLOCATE_IM(p,n,m)
       integer, dimension(:,:), pointer :: p,reallocate_im
       integer, intent(in) :: n,m
       integer :: nold,mold,ierr
       allocate(reallocate_im(n,m),stat=ierr)
       if (ierr/=0) call &
        NRERROR('reallocate_im: problem in attempt to allocate memory')
       if (.not.associated(p)) RETURN
       nold=size(p,1)
       mold=size(p,2)
       reallocate_im(1:min(nold,n),1:min(mold,m))= &
       p(1:min(nold,n),1:min(mold,m))
       deallocate(p)
      end function REALLOCATE_IM

      function OUTERPROD(a,b)
       double precision, dimension(:), intent(in) :: a,b
       double precision, dimension(size(a),size(b)) :: outerprod
       outerprod=spread(a,dim=2,ncopies=size(b))* &
       spread(b,dim=1,ncopies=size(a))
      end function OUTERPROD

      function OUTERDIFF(a,b)
       integer, dimension(:) :: a,b
       integer, dimension(size(a),size(b)) :: outerdiff
       outerdiff = spread(a,dim=2,ncopies=size(b)) - &
       spread(b,dim=1,ncopies=size(a))
      end function OUTERDIFF

      function GET_DIAG(mat)
       !Get diagonal coefficients
       double precision, dimension(:,:) :: mat
       double precision, dimension(size(mat,1)) :: get_diag
       integer :: j
       do j=1,size(mat,1)
        get_diag(j)=mat(j,j)
       end do
      end function GET_DIAG

      subroutine UNIT_MATRIX(mat)
       !Form identity matrix
       double precision, dimension(:,:) :: mat
       integer :: i,n
       n=min(size(mat,1),size(mat,2))
       mat(:,:)=0D0
       do i=1,n
        mat(i,i)=1D0
       end do
      end subroutine UNIT_MATRIX

      function UPPER_TRIANGLE(j,k,extra)
       !Form upper triangle logical mask
       integer :: j,k
       integer, optional :: extra
       logical, dimension(j,k) :: upper_triangle
       integer :: n
       n=0
       if (present(extra)) n=extra
       upper_triangle=(outerdiff(arth_i(1,1,j),arth_i(1,1,k)) < n)
      end function UPPER_TRIANGLE

      subroutine JACOBI(a,d,v,nrot)
       !Diagonalization Jacobi method
       !Computes eigenvalues and eigenvectors
       !of a real symmetric NxN matrix 'a'.
       !Elements of 'a' above diagonal are destroyed.
       !Vector 'd' of length N returns the eigenvalues of 'a'.
       !Matrix NxN 'v' columns contain normalized eigenvectors of 'a'.
       !'nrot' returns the number of Jacobi rotations required.
       !Underflows must be set to zero.
       integer, intent(out) :: nrot
       double precision, dimension(:) :: d
       double precision, dimension(:,:) :: a,v
       integer :: i,ip,iq,n
       double precision :: c,g,h,s,sm,t,tau,theta,tresh
       double precision, dimension(size(d)) :: b,z
       n=size(a,1)
       call UNIT_MATRIX(v(:,:))
       b(:)=GET_DIAG(a(:,:))
       d(:)=b(:)
       z(:)=0D0
       nrot=0
       do i=1,50
        sm=sum(abs(a),mask=UPPER_TRIANGLE(n,n))
        if (sm == 0D0) RETURN
        tresh=merge(0.2D0*sm/n**2,0D0, i < 4 )
        do ip=1,n-1
        do iq=ip+1,n
         g=100D0*abs(a(ip,iq))
         if ((i > 4) .and. (abs(d(ip))+g == abs(d(ip))) &
         .and. (abs(d(iq))+g == abs(d(iq)))) then
          a(ip,iq)=0D0
         else if (abs(a(ip,iq)) > tresh) then
          h=d(iq)-d(ip)
          if (abs(h)+g == abs(h)) then
           t=a(ip,iq)/h
          else
           theta=0.5D0*h/a(ip,iq)
           t=1D0/(abs(theta)+sqrt(1D0+theta**2))
           if (theta < 0D0) t=-t
          end if
          c=1D0/sqrt(1+t**2)
          s=t*c
          tau=s/(1D0+c)
          h=t*a(ip,iq)
          z(ip)=z(ip)-h
          z(iq)=z(iq)+h
          d(ip)=d(ip)-h
          d(iq)=d(iq)+h
          a(ip,iq)=0D0
          call jrotate(a(1:ip-1,ip),a(1:ip-1,iq))
          call jrotate(a(ip,ip+1:iq-1),a(ip+1:iq-1,iq))
          call jrotate(a(ip,iq+1:n),a(iq,iq+1:n))
          call jrotate(v(:,ip),v(:,iq))
          nrot=nrot+1
         end if
        end do
        end do
        b(:)=b(:)+z(:)
        d(:)=b(:)
        z(:)=0D0
       end do
       call nrerror('too many iterations in jacobi')
      contains
       subroutine JROTATE(a1,a2)
        !Jacobi rotation
        double precision, dimension(:), intent(inout) :: a1,a2
        double precision, dimension(size(a1)) :: wk1
        wk1(:)=a1(:)
        a1(:)=a1(:)-s*(a2(:)+a1(:)*tau)
        a2(:)=a2(:)+s*(wk1(:)-a2(:)*tau)
       end subroutine JROTATE
      end subroutine JACOBI

      function IMAXLOC(arr)
       !Maximum location of array
       double precision, dimension(:) :: arr
       integer :: imaxloc
       integer, dimension(1) :: imax
       imax=maxloc(arr(:))
       imaxloc=imax(1)
      end function IMAXLOC

      subroutine SWAP_D(a,b)
       !Swapping of double precision a,b
       double precision :: a,b,dum
       dum=a
       a=b
       b=dum
      end subroutine SWAP_D

      subroutine SWAP_DV(a,b)
       !Swapping of double precision array a,b
       double precision, dimension(:) :: a,b
       double precision, dimension(size(a)) :: dum
       dum=a
       a=b
       b=dum
      end subroutine SWAP_DV

      subroutine EIGSRT(d,v)
       !Given the eigenvalues 'd' and eigenvectors 'v' as
       !output from Jacobi, this routine sorts the
       !eigenvalues into descending order, and rearranges
       !the columns of 'v' correspondingly.
       !For example 1=max, 2=mid, 3=min
       double precision, dimension(:) :: d
       double precision, dimension(:,:) :: v
       integer :: i,j,n
       n=size(d)
       do i=1,n-1
        j=imaxloc(d(i:n))+i-1
        if (j/=i) then
         call swap_d(d(i),d(j))
         call swap_dv(v(:,i),v(:,j))
        end if
       end do
      end subroutine EIGSRT

      subroutine LUDCMP(a,n,indx,d)
       !Given a matrix a(1:n,1:n), this routine
       !replaces it by the LU decompostion of a
       !rowwise permutation of itself. a and n are
       !input. a is output, arranged. indx(1:n) is an
       !output vector that records the row permutation
       !effected by the partial pivoting. d is output as
       !+-1 depending on whether the number of row
       !interchanges was even or odd, respectively.
       !This routine is used with LUBKSB to solve
       !linear equations or invert a matrix.
       integer, intent(in) :: n
       double precision, dimension(n,n), intent(inout) :: a
       integer, dimension(n) , intent(inout) :: indx
       double precision, intent(out) :: d
       double precision, dimension(size(a,1)) :: vv
       double precision, parameter :: TINY=1D-20
       integer :: j,imax
       d=1D0
       vv=maxval(abs(a),dim=2)
       if (any(vv==0D0)) call nrerror('singular matrix in ludcmp')
       vv=1D0/vv
       do j=1,n
        imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j)))
        if (j/=imax) then
         call swap_dv(a(imax,:),a(j,:))
         d=-d
         vv(imax)=vv(j)
        end if
        indx(j)=imax
        if (a(j,j)==0D0) a(j,j)=TINY
        a(j+1:n,j)=a(j+1:n,j)/a(j,j)
        a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
       end do
      end subroutine LUDCMP

      subroutine LUBKSB(a,n,indx,b)
       !Solves the set of n linear equations AX=B. Here
       !a is input, not as matrix A, but as its LU
       !decomposition, determined by subroutine LUDCMP.
       !indx is input as the permutation vector returned
       !by LUDCMP. b(1:n) is input as the right-hand side
       !vector B, and returns with the solution vector X.
       !a,n and indx are not modified by this subroutine.
       integer, intent(in) :: n
       double precision, dimension(n,n) , intent(in) :: a
       integer, dimension(n), intent(in) :: indx
       double precision, dimension(n), intent(inout) :: b
       integer :: i,ii,ll
       double precision :: summ
       ii=0
       do i=1,n
        ll=indx(i)
        summ=b(ll)
        b(ll)=b(i)
        if (ii/=0) then
         summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
        else if (summ/=0D0) then
         ii=i
        end if
        b(i)=summ
       end do
       do i=n,1,-1
        b(i)=(b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
       end do
      end subroutine LUBKSB

      subroutine LUINV(a,n,y)
       !Inverse of matrix a of dimension n.
       !a is destroyed and y returns the inverse of a.
       !This routine uses the LU decomposition.
       integer :: n,i,j
       double precision :: d
       integer, dimension(n) :: indx
       double precision, dimension(n,n) :: a,y
       !Set up identity matrix
       do i=1,n
        do j=1,n
         y(i,j)=0D0
        end do
        y(i,i)=1D0
       end do
       !Decompose the matrix just once
       call LUDCMP(a,n,indx,d)
       !Find inverse by columns
       do j=1,n
        call LUBKSB(a,n,indx,y(1,j))
       end do
      end subroutine LUINV

      subroutine GAULEG(x1,x2,x,w,n)
       !Gaussian quadrature
       !Given the lower and upper limits of integration
       !x1 and x2, and given n, this routine returns
       !arrays x and w of length n, containing the
       !abscissas and weights of the Gauss-Legendre
       !n-point quadrature formula.
       integer :: n
       double precision, intent(in) :: x1,x2
       double precision, dimension(n), intent(out) :: x,w
       double precision, parameter :: EPS=3D-14
       integer :: its,j,m
       integer, parameter :: MAXIT=10
       double precision :: xl,xm
       double precision, dimension((n+1)/2) :: p1,p2,p3,pp,z,z1
       logical, dimension((n+1)/2) :: unfinished
       m=(n+1)/2
       xm=0.5D0*(x2+x1)
       xl=0.5D0*(x2-x1)
       z=cos(acos(-1D0)*(arth_i(1,1,m)-0.25D0)/(n+0.5D0))
       unfinished=.true.
       do its=1,MAXIT
        where (unfinished)
         p1=1D0
         p2=0D0
        end where
        do j=1,n
         where (unfinished)
          p3=p2
          p2=p1
          p1=((2D0*j-1D0)*z*p2-(j-1D0)*p3)/j
         end where
        end do
        where (unfinished)
         pp=n*(z*p1-p2)/(z*z-1D0)
         z1=z
         z=z1-p1/pp
         unfinished=(abs(z-z1)>EPS)
        end where
        if (.not.any(unfinished)) exit
       enddo
       if (its==MAXIT+1) call nrerror('too many iterations in gauleg')
       x(1:m)=xm-xl*z
       x(n:n-m+1:-1)=xm+xl*z
       w(1:m)=2D0*xl/((1D0-z**2)*pp**2)
       w(n:n-m+1:-1)=w(1:m)
      end subroutine GAULEG

      function LOCATE(xx,x)
       !Given a monotonic array xx(1:N) and value x, it
       !returns index j such that xx(j) < x < xx(j+1).
       !If x is out of range, j=0 or j=N.
       double precision, dimension(:), intent(in) :: xx
       double precision, intent(in) :: x
       integer :: locate
       integer :: n,jl,jm,ju
       logical :: ascnd
       n=size(xx)
       ascnd=(xx(n)>=xx(1))
       jl=0
       ju=n+1
       do
        if (ju-jl<=1) exit
        jm=(ju+jl)/2
        if (ascnd.EQV.(x>=xx(jm))) then
         jl=jm
        else
         ju=jm
        end if
       end do
       if (x==xx(1)) then
        locate=1
       else if (x==xx(n)) then
        locate=n-1
       else
        locate=jl
       end if
      end function LOCATE

      function SPLINT(xa,ya,y2a,x)
       !Spline interpolation
       !Given the arrays xa(1:n) and ya(1:n) tabulating
       !a function with the xa's in order and given y2a(1:n)
       !which the output from SPLINE and given a value of x,
       !this returns a cubic-spline interpolation y
       double precision, dimension(:), intent(in) :: xa,ya,y2a
       double precision, intent(in) :: x
       double precision :: splint
       integer :: khi,klo,n
       double precision :: a,b,h
       n=size(xa)
       klo=max(min(LOCATE(xa,x),n-1),1)
       khi=klo+1
       h=xa(khi)-xa(klo)
       if (h==0D0) call NRERROR('bad xa input in splint')
       a=(xa(khi)-x)/h
       b=(x-xa(klo))/h
       splint=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo) &
       +(b**3-b)*y2a(khi))*(h**2)/6D0
      end function SPLINT

      subroutine TRIDAG_SER(a,b,c,r,u)
       !Tridiagonal serial algorithm
       !Solves for a vector u(1:N) the tridiagonal system
       !Input vector b(1:N) is diagonal elements and r(1:N)
       !is right hand side.
       !a(1:N-1) and c(1:N-1) are off-diagonal elements
       double precision, dimension(:), intent(in) :: a,b,c,r
       double precision, dimension(:), intent(out) :: u
       double precision, dimension(size(b)) :: gam
       integer :: n,j
       double precision :: bet
       n=size(a)+1
       bet=b(1)
       if (bet==0D0) call NRERROR('tridag_ser: Error at code stage 1')
       u(1)=r(1)/bet
       !Decomposition and forward substitution
       do j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j-1)*gam(j)
        if (bet==0D0) call NRERROR('tridag_ser: Error at code stage 2')
        u(j)=(r(j)-a(j-1)*u(j-1))/bet
       end do
       !Backsubstitution
       do j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
       end do
      end subroutine TRIDAG_SER

      recursive subroutine TRIDAG_PAR(a,b,c,r,u)
       !Tridiagonal parallel algorithm
       !Solves for a vector u(1:N) the tridiagonal system
       !Input vector b(1:N) is diagonal elements and r(1:N)
       !is right hand side.
       !a(1:N-1) and c(1:N-1) are off-diagonal elements
       double precision, dimension(:), intent(in) :: a,b,c,r
       double precision, dimension(:), intent(out) :: u
       integer, parameter :: NPAR_TRIDAG=4
       integer :: n,n2,nm,nx
       double precision, dimension(size(b)/2) :: y,q,piva
       double precision, dimension(size(b)/2-1) :: x,z
       double precision, dimension(size(a)/2) :: pivc
       n=size(a)+1
       if (n<NPAR_TRIDAG) then
        call TRIDAG_SER(a,b,c,r,u)
       else
        if (maxval(abs(b(1:n)))==0D0) then
         call NRERROR('tridag_par: possible singular matrix')
        end if
        n2=size(y)
        nm=size(pivc)
        nx=size(x)
        piva=a(1:n-1:2)/b(1:n-1:2)
        pivc=c(2:n-1:2)/b(3:n:2)
        y(1:nm)=b(2:n-1:2)-piva(1:nm)*c(1:n-2:2)-pivc*a(2:n-1:2)
        q(1:nm)=r(2:n-1:2)-piva(1:nm)*r(1:n-2:2)-pivc*r(3:n:2)
        if (nm<n2) then
         y(n2)=b(n)-piva(n2)*c(n-1)
         q(n2)=r(n)-piva(n2)*r(n-1)
        end if
        x=-piva(2:n2)*a(2:n-2:2)
        z=-pivc(1:nx)*c(3:n-1:2)
        call TRIDAG_PAR(x,y,z,q,u(2:n:2))
        u(1)=(r(1)-c(1)*u(2))/b(1)
        u(3:n-1:2)=(r(3:n-1:2)-a(2:n-2:2)*u(2:n-2:2) &
         -c(3:n-1:2)*u(4:n:2))/b(3:n-1:2)
        if (nm==n2) u(n)=(r(n)-a(n-1)*u(n-1))/b(n)
       end if
      end subroutine TRIDAG_PAR

      subroutine SPLINE(x,y,n,yp1,ypn,y2)
       !Spline interpolation second derivatives
       !Given arrays x(1:n) and y(1:n) containing
       !a tabulated function, i.e., yi=f(xi), with
       !x1<x2<...<xn, and given values yp1 and ypn
       !for the first derivative of the interpolating
       !function at points 1 and n, this returns an
       !array y2(1:n) of length n containing the second
       !derivatives of the interpolating function at the
       !tabulated points xi
       double precision, dimension(:), intent(in) :: x,y
       double precision, intent(in) :: yp1,ypn
       double precision, dimension(:), intent(out) :: y2
       integer :: n
       double precision, dimension(size(x)) :: a,b,c,r
       n=size(x)
       c(1:n-1)=x(2:n)-x(1:n-1)
       r(1:n-1)=6D0*((y(2:n)-y(1:n-1))/c(1:n-1))
       r(2:n-1)=r(2:n-1)-r(1:n-2)
       a(2:n-1)=c(1:n-2)
       b(2:n-1)=2D0*(c(2:n-1)+a(2:n-1))
       b(1)=1D0
       b(n)=1D0
       if (yp1>0.99D30) then
        r(1)=0D0
        c(1)=0D0
       else
        r(1)=(3D0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
        c(1)=0.5D0
       end if
       if (ypn>0.99D30) then
        r(n)=0D0
        a(n)=0D0
       else
        r(n)=(-3D0/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
        a(n)=0.5D0
       end if
       call TRIDAG_PAR(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
      end subroutine SPLINE

      end module NRMOD
