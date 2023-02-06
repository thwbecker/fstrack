
! Module library for the Lambert projection (equal area
! projection) of the orientation distribution function
! Normalized by mean random deviation (i.e. random ODF)
!
! LAMBERT_GRID  (subroutine) : forms spherical net
! LAMBERT_PROJ  (subroutine) : Lambert projection
! LAMBERT_GMT   (subroutine) : ouput for GMT
! LAMBERT_MRD   (subroutine) : mean random deviation
!
! by Jules Browaeys, as of August 2005
!
! modifications by TWB
!
! $Id: drex_lambert.f90,v 1.2 2005/08/11 22:28:58 becker Exp becker $
!
#include "drex_fconst.h"
  
subroutine drex_lambert_proj(ngr,lacs,lodf,rotmat,rotate,nlat,ldist)
  !
  ! Lambert projection on lower hemisphere
  ! Each point has a gaussian amplitude
  ! a=[100],b=[010],c=[001]
  ! Longitude = netpsi [0,2*pi]
  ! Colatitude = netheta [0,pi]
  !
  !
  integer, intent(in) :: ngr  ! number of grains
  integer, intent(in) :: nlat ! number of points in latitude
  double precision, intent(inout), dimension(3,nlat,nlat*2) :: ldist ! distribution grid
  double precision, intent(in), dimension(3,3) :: rotmat ! rotation matrix to apply before summation
  logical, intent(in) :: rotate ! should we rotate?
  double precision, intent(in), dimension(ngr) :: lodf ! orientation density function for each grain
  double precision, intent(in), dimension(ngr,3,3):: lacs ! direction cosines, lacs(igrain,iaxis,idir)
  !
  ! ldist(:,nlat,1) will hold the polar value
  !

  !
  ! normalization for random distribution
  double precision, parameter :: invmrd=1D0/2.35218456950443D-2
  !
  ! local
  !
  integer :: n,ia,i,j,is,js,nlon,nlatm1,nlathm1
  logical :: init
  double precision :: xp,xd,sxd,ldistn,df,netpsi,sin_lat,fac1,fac2
  double precision, dimension(3) :: vec1,vec2
  double precision, dimension(3,0:1) :: ldistp
  double precision, dimension(:,:), allocatable :: lat,long
  double precision, dimension(:), allocatable :: sin_netheta,cos_netheta,cos_netpsi,sin_netpsi
  save sin_netheta,cos_netheta,cos_netpsi,sin_netpsi,init
  data init /.false./
  
  !
  ! bounds
  !
  nlatm1 = nlat - 1
  nlathm1 = nlat/2 - 1
  nlon = nlat * 2
  if(.not.init)then
     !
     ! init trig factors
     !
     print *,'drex_lambert_proj: initializing for nlat:',nlat
     allocate(sin_netpsi(nlon));allocate(cos_netpsi(nlon));
     allocate(sin_netheta(nlat));allocate(cos_netheta(nlat));
     df = DREX_PI/nlat
     netpsi=0
     do i=1, nlon
        netpsi = netpsi + df
        sin_netpsi(i) = sin(netpsi)
        cos_netpsi(i) = cos(netpsi)
        if (i < nlat)then       ! address only up to nlat-1
           cos_netheta(i) = cos_netpsi(i)
           sin_netheta(i) = sin_netpsi(i)
        endif
     end do
     init = .true.
  endif
  ! 
  ! local coordinate arrays
  !
  allocate(lat(3,ngr))
  allocate(long(3,ngr))
  !
  ! Gaussian characteristics
  !
  ! this was 4*cp180, but cp180 was never defined? XXX
  sxd=4D0*(DREX_PI/180.)/sqrt(log(2D0))
  ldistn=1D0/sqrt(DREX_PI)/sxd
  sxd=1D0/sxd**2
  !
  !
  
  !
  ! Grains position on sphere
  !
  do n=1,ngr
     !
     ! Axes 1,2,3
     !
     do ia=1,3
        !
        ! Latitude and longitude on sphere
        !
        vec1 = lacs(n,ia,:)
        if(rotate)then          ! need to rotate this grain
           call drex_abase_vec2bbase_vec_ftrn(vec1,vec2,rotmat)
        else
           vec2 = vec1
        endif
        if ((vec2(2) == 0D0).and.(vec2(1) == 0D0)) then
           long(ia,n) = 0D0
           if (vec2(3) >= 0D0) lat(ia,n)=0D0
           if (vec2(3) <= 0D0) lat(ia,n)=DREX_PI
        elseif (vec2(3) >= 1D0) then
           long(ia,n)=0D0
           lat(ia,n)=0D0
        elseif (vec2(3) <= -1D0) then
           long(ia,n)= 0D0
           lat(ia,n) = DREX_PI
        else
           lat(ia,n) =acos(vec2(3)) ! theta
           long(ia,n)=atan2(vec2(2),vec2(1)) !azimuth
           if (long(ia,n) <= 0D0) long(ia,n) = 2D0*DREX_PI + long(ia,n)
        endif
        !
        ! Restrict to upper hemisphere
        !
        if (lat(ia,n) > DREX_HALF_PI) then
           lat(ia,n) =  DREX_PI - lat(ia,n) 
           long(ia,n) = modulo(DREX_PI + long(ia,n),DREX_TWO_PI)
        endif
     end do
  end do                        ! end computation of lon/lat
  !
  ! Lambert projection amplitude
  !
  ldist=0D0                     !init
  do j=1, nlon              ! longitude loop
     do i=1, nlatm1           ! latitude loop
        fac1 = sin_netheta(i) * cos_netpsi(j)
        fac2 = sin_netheta(i) * sin_netpsi(j)
        do n=1, ngr           ! grain loop 
           do ia=1, 3         ! axis loop
              !
              ! Gaussian distribution on sphere
              !
              sin_lat = sin(lat(ia,n))
              xp =      fac1 * sin_lat * cos(long(ia,n))  
              xp = xp + fac2 * sin_lat * sin(long(ia,n))
              xp = xp + cos_netheta(i) * cos(lat(ia,n))
              ! compute distance squared
              if (xp >= 1D0) then
                 xd = 0D0
              else if (xp<=-1D0) then
                 xd = DREX_PI**2
              else
                 xd= (acos(xp))**2
              endif
              ldist(ia,i,j) = ldist(ia,i,j) + lodf(n) * exp(-xd*sxd)
           end do
        end do
     end do
  end do
  !
  ! normalize
  !
  ldist = ldist * ldistn
  !
  ! Poles
  !
  ldistp = 0D0
  ! Theta = 0
  do n=1,ngr
     do ia=1,3
        ! Gaussian distribution on sphere
        xp=cos(lat(ia,n))
        if (xp>=1D0) then
           xd=0D0
        else if (xp<=-1D0) then
           xd=DREX_PI**2
        else
           xd=(acos(xp))**2
        endif
        ldistp(ia,0)=ldistp(ia,0)+lodf(n)*exp(-xd*sxd)
     end do
  end do
  !
  ! Theta = pi
  !
  do n=1,ngr
     do ia=1,3
        ! Gaussian distribution on sphere
        xp = -cos(lat(ia,n))
        if (xp >= 1D0) then
           xd = 0D0
        else if (xp <= -1D0) then
           xd=DREX_PI**2
        else
           xd=(acos(xp))**2
        endif
        ldistp(ia,1)=ldistp(ia,1)+lodf(n)*exp(-xd*sxd)
     end do
  end do
  ! normalize
  ldistp = ldistp * ldistn
  deallocate(long)
  deallocate(lat)
  !
  ! Restrict Lambert projection to lower hemisphere
  !
  ldistp(:,1) =ldistp(:,1) + ldistp(:,0)
  !
  ! store pole in ldist
  !
  ldist(:,nlat,1) =  ldistp(:,1)
  
  do i=1,nlathm1
     is=nlat-i
     do j=1,nlon
        js=j+nlat
        if (j > nlat) js = j - nlat
        ldist(:,is,js) = ldist(:,is,js) + ldist(:,i,j)
     end do
  end do
  do j=1, nlon
     ldist(:,nlathm1,j) = ldist(:,nlat/2+1,j)
  end do
  !
  ! Normalization by mean random deviation (i.e. mrd)
  ldist = ldist * invmrd
end subroutine DREX_LAMBERT_PROJ



subroutine drex_lambert_gmt(nframe,inmx,nlat,ldist)
  double precision, intent(in), dimension(3,nlat,nlat*2) :: ldist ! distribution grid
  integer, intent(in) :: nlat
  double precision, parameter :: ldmin=1D-6
  integer :: ia,i,j,nlatm1
  integer :: nframe,inmx
  integer :: ilength,nval,i0,isl,il
  double precision :: ldgmt,ldgmtp,df
  character(len=200) :: file1,filegrd

  nlatm1 = nlat - 1
  nlon = nlat * 2 
  file1='gmtvisual/L'

  if (inmx==1) file1='gmtvisual/L.M1.'
  if (inmx==2) file1='gmtvisual/L.M2.'
  if (nframe==0) file1=file1(1:len_trim(file1))//'0'
  if (nframe>0) then
     ilength=1
     nval=10
     do while (nframe>=nval)
        ilength=ilength+1
        nval=10*nval
     end do
     i0=iachar('0')
     isl=nframe
     do il=1,ilength
        nval=nval/10
        i=isl/nval
        file1=file1(1:len_trim(file1))//achar(i0+i)
        isl=mod(isl,nval)
     end do
  endif
  df = 90.0/nlat 

  do ia=1,3
     if (ia==1) filegrd=file1(1:len_trim(file1))//'.100.dat'
     if (ia==2) filegrd=file1(1:len_trim(file1))//'.010.dat'
     if (ia==3) filegrd=file1(1:len_trim(file1))//'.001.dat'
     open(21,file=filegrd,status='replace')

     ! pole
     ldgmtp=ldist(ia,nlat,1)
     
     if (ldgmtp<ldmin) ldgmtp=0D0
     ! Longitude=0 
     do i=nlat/2-1,nlatm1
        ldgmt=ldist(ia,i,nlon)
        if (ldgmt < ldmin) ldgmt=0D0
        write(21,'(1x,e15.7,2x,e15.7,2x,e15.7)') &
             0D0,90D0 - i * df,ldgmt
     end do
     write(21,'(1x,e15.7,2x,e15.7,2x,e15.7)') &
          0D0,-90D0,ldgmtp
     ! Longitude>0
     do j=1,nlon
        do i=nlat/2-1,nlatm1
           ldgmt=ldist(ia,i,j)
           if (ldgmt < ldmin) ldgmt=0D0
           write(21,'(1x,e15.7,2x,e15.7,2x,e15.7)') &
                j*df,90D0-i*df,ldgmt
        end do
        write(21,'(1x,e15.7,2x,e15.7,2x,e15.7)') &
             j*df, -90D0, ldgmtp
     end do
     close(21)
  end do
end subroutine DREX_LAMBERT_GMT
  
!
! compute mean random deviation from ldist
!
  subroutine drex_lambert_mrd(mrd,ldist,nlat)
    integer, intent(in) :: nlat
    double precision, intent(in), dimension(3,nlat,nlat*2) :: ldist
    double precision, intent(out) :: mrd


    integer :: ia,i,j,nlon,nlath

    nlon = nlat * 2
    nlath = nlat / 2
    mrd = 0D0
    do j=1,nlon
       do i=nlath,nlat
          do ia=1, 3
             mrd = mrd + ldist(ia,i,j)
          end do
       end do
    end do
    mrd=mrd/(3D0*(2.0d0*nlat)*(nlat/2+1))
    !write(*,*) ' mean random deviation =',mrd
  end subroutine DREX_LAMBERT_MRD
  

