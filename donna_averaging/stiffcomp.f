
      SUBROUTINE stiffcomp(c,s)
c
c  Subroutine to convert c(i,j) to s(i,j)
c  or vice-versa, tensor to be converted first

c  From Dan Raymer 3/00; minor mods by dkb to
c  clean up double precision lingo
c
      PARAMETER (np=6)
      IMPLICIT real*8(a-h,o-z) 
      INTEGER n,np,indx(np)
      DIMENSION c(np,np),s(np,np),rmat(np,np)
      n=6
      do i=1,n 
        do j=1,n 
           rmat(i,j)=c(i,j) 
           s(i,j)=0.d0
        enddo 
        s(i,i)=1.d0 
      enddo
c
c  Actual inversion routine
c
      call ludcmp(rmat,n,np,indx,d)
      do j=1,n            
        call lubksb(rmat,n,np,indx,s(1,j))
      enddo
c
      return
      END
c
c=========================================================================
c
      SUBROUTINE ludcmp(a,n,np,indx,d)
      IMPLICIT real*8(a-h,o-z) 
      INTEGER n,np,indx(n),NMAX
      REAL*8 a(np,np)
      PARAMETER (NMAX=500,TINY=1.0d-20)
      INTEGER i,imax,j,k
      REAL*8 vv(NMAX)

      d=1.d0

      do 12 i=1,n
        aamax=0.d0
        do 11 j=1,n
          if (dabs(a(i,j)).gt.aamax) aamax=dabs(a(i,j))
11      continue
        if (aamax.eq.0.d0) pause 'singular matrix in ludcmp'
        vv(i)=1.d0/aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.d0
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*dabs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.d0)a(j,j)=TINY
        if(j.ne.n)then
          dum=1.d0/a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software *%&&,1{.


      SUBROUTINE lubksb(a,n,np,indx,b)
      IMPLICIT real*8(a-h,o-z) 
      INTEGER n,np,indx(n)
      REAL*8 a(np,np),b(n)
      INTEGER i,ii,j,ll

      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.d0) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END
c  (C) Copr. 1986-92 Numerical Recipes Software *%&&,1{.

