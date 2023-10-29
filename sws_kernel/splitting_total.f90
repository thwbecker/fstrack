
!----------------------------------------------------------------------------
! This program calculates the splitting intensita at the surface produced by
! the hexagonal anisotropy distribution described by any orientation of the
! symmetry axis and by the two anisotropic parameters
! 
!
! from Sebastien Chevrot as of June 2005
!
! slight modifications by TWB

!
! $Id: splitting_total.f90,v 1.5 2005/11/10 23:37:27 becker Exp becker $
!

!----------------------------------------------------------------------------

!
!
!
!
program splitting_total

  implicit none

  ! For loop on backazimuths

  integer, parameter :: pasx = 2, pasy = 2, pasz = 2
  integer, parameter :: nptmax = 8000000
  !
  double precision, parameter ::  pas_alpha = 2.d0
  integer, parameter :: length_alpha = 360 
  integer, parameter :: ind_alpha=length_alpha/pas_alpha+1
  !

  integer :: i,ind,ib,npt,Nx,Ny,Nz,idepth,k,off,j,ixs,iys,indref,irc,izskip

  double precision :: x(nptmax),y(nptmax),z(nptmax),r(nptmax),ro(nptmax)
  double precision :: azi(nptmax),inc(nptmax),cos_inc(nptmax),sin_inc(nptmax),&
       cos_azi(nptmax),sin_azi(nptmax),cos_phi(nptmax),sin_phi(nptmax)
  double precision :: gamma(nptmax),epsilon(nptmax),delta(nptmax)

  double precision :: A_FF,A_MF,A_NF,A_LF
  double precision :: herm_FF,herm_MF,herm_NF,herm_LF
  double precision :: F_gamma_FF,F_gamma_MF,F_gamma_NF,F_gamma_LF
  double precision :: F_epsilon_FF,F_epsilon_MF,F_epsilon_NF,F_epsilon_LF
  double precision :: F_delta_FF,F_delta_MF,F_delta_NF,F_delta_LF

  double precision :: deg2rad,beta,tau,lambda,aa
  double precision :: kr,zmin,zmax
  double precision :: phi,sinth,costh,alpha,inci
  double precision :: xsta,ysta,slon,slat,error
  double precision :: p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12
  double precision :: k1(3),k2(3),t(3),p(3),s(3)
  double precision :: IS(ind_alpha),Kz1(100),Kz2(100),Kz3(100),Kz4(100)
  double precision :: dxvol,zold,epsmean,dx,dy,dz,xloc,yloc,zloc,xmin,xmax,ymin,ymax

  character (len=200)    :: mode_string

  double precision, parameter :: pi = 3.141592653589793d0,two_pi = 2.0d0 * pi,&
       sqrt_two = 1.4142135623730950d0

  double precision, parameter :: dlmin = 0.0d0,dlmax = 50.0d0
  !
  logical, parameter :: print_kernel=.false.
  !
  integer, parameter :: icas=1  ! 1: Gaussian 2: Ricker
  !
  logical, parameter :: large_region = .true. ! large: -200...200 not: -100...100
  !logical, parameter :: large_region = .false. 

  ! if true, will make all parameters only vary with depth according to the station
  ! in the middle
  integer :: fake_even

  fake_even = 0                 ! 0: use all data 
  ! 1: distribute center data evenly for each layer 
  ! 2: use hexagonal
  ! 3: use hexgonal for dlmin< z <=50km, isotropic else
  ! 4: use 45 dipping hexgonal for dlmin < z<=dlmax, isotropic else
  ! 5: use 30 dipping hexgonal for dlmin < z<=dlmax, isotropic else

  if(iargc().gt.0)then
     call getarg(1,mode_string)
     read(mode_string,*)fake_even
  endif

  print *,'splitting_total: using mode',fake_even



  deg2rad = pi/180.d0

  ! station location 

  slon = 242.d0
  slat = 34.d0

  error = 0.1d0

  !
  ! output
  if(fake_even.eq.1)then
     print *,'testing mode: using center station for all nodes in layer'
     print *,'writing to splitting.even.dat'
     open(200,file='splitting.even.dat') !splitting intensity
     open(300,file='kernel.even.dat',form='unformatted',status='unknown') ! 3-D kernel
     open(400,file='profil.even.dat')   ! kernel profile
  else if(fake_even.eq.2)then
     print *,'testing mode: using single hexagonal tensor '
     print *,'writing to splitting.hex.dat'
     open(200,file='splitting.hex.dat') !splitting intensity
     open(300,file='kernel.hex.dat',form='unformatted',status='unknown') ! 3-D kernel
     open(400,file='profil.hex.dat')   ! kernel profile
  else if(fake_even.eq.3)then
     print *,'testing mode: using single hexagonal tensor for z<=50'
     print *,'writing to splitting.tophex.dat'
     open(200,file='splitting.tophex.dat') !splitting intensity
     open(300,file='kernel.tophex.dat',form='unformatted',status='unknown') ! 3-D kernel
     open(400,file='profil.tophex.dat')   ! kernel profile
  else if(fake_even.eq.4)then
     print *,'testing mode: using single 45 dipping hexagonal tensor for z<=50'
     print *,'writing to splitting.tophexdip.dat'
     open(200,file='splitting.tophexdip.dat') !splitting intensity
     open(300,file='kernel.tophexdip.dat',form='unformatted',status='unknown') ! 3-D kernel
     open(400,file='profil.tophexdip.dat')   ! kernel profile
  else if(fake_even.eq.5)then
     print *,'testing mode: using single 30 dipping hexagonal tensor for z<=50'
     print *,'writing to splitting.tophexdip2.dat'
     open(200,file='splitting.tophexdip2.dat') !splitting intensity
     open(300,file='kernel.tophexdip2.dat',form='unformatted',status='unknown') ! 3-D kernel
     open(400,file='profil.tophexdip2.dat')   ! kernel profile
  else
     print *,'regular mode, using all data as read in'
     print *,'writing to splitting.dat'
     open(200,file='splitting.dat') !splitting intensity
     open(300,file='kernel.dat',form='unformatted',status='unknown') ! 3-D kernel
     open(400,file='profil.dat')   ! kernel profile
  endif

  ! Reading the reference model
  beta = 4.47d0

  !tau = 10.0d0                   ! for SKS
  tau = 12.5d0
  !tau = 15.d0
  !tau = 30.d0                   

  lambda = beta*tau
  !
  ! incidence of wave
  !inci = 15.d0*deg2rad
  inci = 5.d0*deg2rad

  ! Factors related to the choice of teh wave form
  if (icas.eq.1) then
     ! Gaussian wavelet
     A_FF = 1.d0/(2.d0*beta**2*tau*sqrt_two)
     A_MF = 1.d0/(4.d0*pi*beta)
     A_NF = tau/(4.d0*sqrt_two*pi**2)
     A_LF = beta*tau**2/(8.d0*pi**3)
  elseif (icas.eq.2) then
     ! Ricker wavelet, 2nd derivative of Gaussian
     A_FF = 1.d0/(120.d0*beta**2*tau*sqrt_two)
     A_MF = 1.d0/(240.d0*pi*beta)
     A_NF = tau/(240.d0*sqrt_two*pi**2)
     A_LF = beta*tau**2/(480.d0*pi**3)
  endif
  ! Definition of station position
  !  xsta = 200.d0
  !  ysta = 200.d0
  xsta = 0.d0
  ysta = 0.d0

  ! Reading the anisotropic model

  print *,'reading the model...'

  !  read(100) dx,dy,dz
  !  read(100) nx,ny,nz
  !  read(100) zmin
  !  read(100) zmax
  dx = pasx;
  dy = pasy;
  dz = pasz;
  if(large_region)then
     nx=201                   !-200:2:200
     ny=201                   !-200:2:200
     xmin=-200d0;xmax=200d0
     ymin=-200d0;ymax=200d0
  else
     nx=101                   !-100:2:100
     ny=101                   !-100:2:100
     xmin=-100d0;xmax=100d0
     ymin=-100d0;ymax=100d0
  endif
  
  nz=175;                  ! 2:2:350
  zmin=2d0;zmax=350d0
  ! center
  ixs=(nx-1)/2;iys=(ny-1)/2
  
  ! volumen
  dxvol = dx*dy*dz


  write(300) dx,dy,dz
  write(300) nx,ny,nz
  write(300) zmin
  write(300) zmax
  print*,dx,dy,dz
  print*,nx,ny,nz
  print*,zmin
  print*,zmax

  npt = nx * ny * nz      
  if(fake_even.le.1)then
     print *,'reading model'
     ! Open files containing the anisotropic model and
     ! the splitting intensity depending on the backazimuth
     !open(100,file='aniso.mod',form='unformatted',status='old')
     open(100,file='inp.ref',form='binary')
     irc=0
     epsmean=0.0d0
     do ind =1,npt
        ! this now expects radians for azi and inc
        read(100) x(ind),y(ind),z(ind),azi(ind),inc(ind),  &
             gamma(ind),epsilon(ind),delta(ind)
        epsmean=epsmean+epsilon(ind)
        irc=irc+1
        !if(mod(ind,nx*ny).eq.1)print *,z(ind)
        if(z(ind).ne.zold)then
           zold=z(ind)
           print *,zold,epsmean/(nx*ny)
           epsmean=0d0
        endif
     enddo
1    continue
     close(100)

     if(irc.ne.npt)then
        print *,'read error: read :',irc
        print *,'expected ntotal:nx,ny,nz ',npt,nx,ny,nz
        stop
     end if
     if(irc.eq.nptmax)then
        print *,'irc ',irc
        print *,'reached nptmax, increase!',nptmax
        stop
        
     end if

     print *,'done reading model'
  else
     print *,'faking structure'
     do i=1,nz
        off = (i-1)*nx*ny
        zloc = zmin+dble(i-1)*dz
        do j=1,ny
           yloc = ymin+dble(j-1)*dy
           do k = 1,nx
              ind = off + (j-1)*nx + k
              x(ind) = xmin+dble(k-1)*dx
              y(ind) = yloc
              z(ind) = zloc
           enddo
        enddo
     enddo
     azi=0.0d0
     inc=0.0d0
     gamma=0.0d0
     epsilon=0.0d0
     delta=0.0d0
  endif
  !open(500,file='cdata.dat')   ! converted data
  !do ind=1,npt
  !   write(500,113) x(ind),y(ind),z(ind),azi(ind),inc(ind),  &
  !        gamma(ind),epsilon(ind),delta(ind)
  !enddo
  !close(500)
  !print *,'converted layer input in cdata.dat'
  if(fake_even.eq.1)then
     !
     ! use the central data values at each layer to reassign
     !
     print *,'faking layer structure!!!'
     do i=1,nz
        off = (i-1)*nx*ny
        indref = off + (iys-1)*nx + ixs
        do j=1,ny
           do k = 1,nx
              ind = off + (j-1)*nx + k
              if(ind.ne.indref)then !reassign
                 azi(ind) = azi(indref)
                 inc(ind) = inc(indref)
                 gamma(ind) = gamma(indref)
                 epsilon(ind) = epsilon(indref)
                 delta(ind) = delta(indref)
              endif
           enddo
        enddo
     enddo
     !open(500,file='cdata.dat')   ! converted data
     !do ind=1,npt
     !   write(500,113) x(ind),y(ind),z(ind),azi(ind),inc(ind),  &
     !        gamma(ind),epsilon(ind),delta(ind)
     !enddo
     !close(500)
     !print *,'converted layer input in cdata.dat'

  else if(fake_even.eq.2)then
     !
     print *,'using single hex tensor everywhere!!!'
     !
     ! use NS oriented hexagonal anisotropy without dip
     do i=1,npt
        azi(i) = 0.0d0
        inc(i) = 0.0d0
        epsilon(i)= -0.096541d0          ! 10% anisootrpy
        gamma(i) =  -0.103181d0
        delta(i) =  -0.103182d0
     enddo
  else if(fake_even.eq.3)then
     print *,'using single hex tensor in top layer'
     !
     ! single hex for z<=50km
     !
     do i=1,npt
        azi(i) = 0.0d0
        inc(i) = 0.0d0
        if((z(i).le.dlmax).and.(z(i).gt.dlmin))then
           epsilon(i)= -0.096541d0          ! 10% hex anisootrpy
           gamma(i) =  -0.103181d0
           delta(i) =  -0.103182d0
        else
           epsilon(i)= 0.0d0    !no anisotropy deeper
           gamma(i) =  0.0d0
           delta(i) =  0.0d0
        endif
     enddo
  else if(fake_even.eq.4)then
     print *,'using single 45 dipping hex tensor in top layer'
     !
     ! single hex for z<=50km
     !
     do i=1,npt
        azi(i) = 0.0d0
        inc(i) = 45.d0*deg2rad
        if((z(i).le.dlmax).and.(z(i).gt.dlmin))then
           epsilon(i)= -0.096541d0          ! 10% hex anisootrpy
           gamma(i) =  -0.103181d0
           delta(i) =  -0.103182d0
        else
           epsilon(i)= 0.0d0    !no anisotropy deeper
           gamma(i) =  0.0d0
           delta(i) =  0.0d0
        endif
     enddo

  else if(fake_even.eq.5)then
     print *,'using single 30 dipping hex tensor in top layer'
     !
     ! single hex for z<=50km
     !
     do i=1,npt
        azi(i) = 0.0d0
        inc(i) = 30.d0*deg2rad
        if((z(i).le.dlmax).and.(z(i).gt.dlmin))then
           epsilon(i)= -0.096541d0          ! 10% hex anisootrpy
           gamma(i) =  -0.103181d0
           delta(i) =  -0.103182d0
        else
           epsilon(i)= 0.0d0    !no anisotropy deeper
           gamma(i) =  0.0d0
           delta(i) =  0.0d0
        endif
     enddo

  endif
  !
  ! precompute a few things
  !
  do i=1,npt
     r(i) = sqrt((x(i)-xsta)**2+(y(i)-ysta)**2+z(i)**2)
     if(abs(r(i)).lt.1e-7)then
        print *,'error: r=',r(i)
        print *,'station: ',xsta,ysta
        print *,'obs: ',x(i),y(i),z(i)
        print *,'i = ',i
        stop
     endif
     ro(i) = sqrt((x(i)-xsta)**2+(y(i)-ysta)**2)
     !
     cos_inc(i) = cos(inc(i))
     sin_inc(i) = sin(inc(i))
     cos_azi(i) = cos(azi(i))
     sin_azi(i) = sin(azi(i))

     if (abs(x(i)-xsta) > 0.d0) then
        !           phi = atan2((y(i)-ysta),(x(i)-xsta))
        phi = atan2((x(i)-xsta),(y(i)-ysta))
     else
        phi = 0.d0
     endif
     cos_phi(i) = cos(phi)
     sin_phi(i) = sin(phi)
  enddo

  ! Initialization of Kz
  do i = 1,Nz
     Kz1(i) = 0.d0
     Kz2(i) = 0.d0
     Kz3(i) = 0.d0
     Kz4(i) = 0.d0
  end do


  ! Loop on backazimuths
  do ib = 1,ind_alpha

     print *,ib,' out of ',ind_alpha, ' backazimuths'

     alpha = (dble(ib-1)*pas_alpha)*deg2rad
     !alpha = 45.d0*deg2rad
     k1(1) = -sin(inci)*cos(alpha)
     k1(2) = -sin(inci)*sin(alpha)
     k1(3) = -cos(inci)
     p(1) = cos(inci)*cos(alpha)
     p(2) = cos(inci)*sin(alpha)
     p(3) = -sin(inci)
     t(1) = -sin(alpha)
     t(2) = cos(alpha)

     ! init splitting intensity
     IS(ib) = 0.d0

     ! Loop on all grid points
     do i = 1,npt
        sinth = z(i)/r(i)
        costh = ro(i)/r(i)
     
        !        kr = abs(k1(1)*(x(i)-xsta)+k1(2)*(y(i)-ysta)+k1(3)*z(i))
        kr = abs(k1(1)*(y(i)-ysta)+k1(2)*(x(i)-xsta)+k1(3)*z(i))
        aa = sqrt_two*(r(i)-kr)*pi/lambda

        if (aa < two_pi) then
           ! changed, as azi and inc are now in radians
           s(1) = cos_inc(i)*cos_azi(i)
           s(2) = cos_inc(i)*sin_azi(i)
           s(3) = sin_inc(i)
           !
           k2(1) = -costh*cos_phi(i)
           k2(2) = -costh*sin_phi(i)
           k2(3) = -sinth
           !
           p1 = k1(1)*k2(1)+k1(2)*k2(2)+k1(3)*k2(3)
           p2 = k2(1)*p(1)+k2(2)*p(2)+k2(3)*p(3)
           p3 = k2(1)*t(1)+k2(2)*t(2)
           p4 = p(1)*s(1)+p(2)*s(2)+p(3)*s(3)
           p5 = k1(1)*s(1)+k1(2)*s(2)+k1(3)*s(3)
           p6 = k2(1)*s(1)+k2(2)*s(2)+k2(3)*s(3)
           p7 = t(1)*s(1)+t(2)*s(2)

           ! Les Polynomes d'Hermite
           if (icas.eq.1) then
              herm_FF = (8.d0*aa**3-12.d0*aa)*exp(-(aa**2))
              herm_MF = (4.d0*aa**2-2.d0)*exp(-(aa**2))
              herm_NF = 2.d0*aa*exp(-(aa**2))
              herm_LF = exp(-(aa**2))
           elseif (icas.eq.2) then
              herm_FF = (128.d0*aa**7-1344.d0*aa**5+3360.d0*aa**3 -1680.d0*aa)*exp(-(aa**2))
              herm_MF = (64.d0*aa**6-480.d0*aa**4+720.d0*aa**2-120.d0)*exp(-(aa**2))
              herm_NF = (32.d0*aa**5-160.d0*aa**3+120.d0*aa)*exp(-(aa**2))
              herm_LF = (16.d0*aa**4-48.d0*aa**2+12.d0)*exp(-(aa**2))
           endif

          ! Termes gamma
           !      F_gamma_FF = ((p4*p1+p5*p2)*(4.d0*p3*p6-2.d0*p7)-2.d0*p3*p2*p1)/r
           !      F_gamma_FF = ((p4*p1+p5*p2)*(4.d0*p3*p6-2.d0*p7)-4.d0*p3*p2*p1)/r
           F_gamma_FF = (p4*p1+p5*p2)*(4.d0*p3*p6-2.d0*p7)/r(i)
           !      F_gamma_MF = ((p4*p1+p5*p2)*(24.d0*p3*p6-6.d0*p7)-p3*(8.d0*p5*p4+12.d0*p2*p1))/r**2
           !      F_gamma_MF = ((p4*p1+p5*p2)*(24.d0*p3*p6-6.d0*p7)-p3*(8.d0*p5*p4+24.d0*p2*p1))/r**2
           F_gamma_MF = ((p4*p1+p5*p2)*(24.d0*p3*p6-6.d0*p7)-8.d0*p3*p5*p4)/(r(i)**2)
           !      F_gamma_NF = ((p4*p1+p5*p2)*(60.d0*p3*p6-12.d0*p7)-p3*(24.d0*p5*p4+30.d0*p2*p1))/r**3
           !      F_gamma_NF = ((p4*p1+p5*p2)*(60.d0*p3*p6-12.d0*p7)-p3*(24.d0*p5*p4+60.d0*p2*p1))/r**3
           F_gamma_NF = ((p4*p1+p5*p2)*(60.d0*p3*p6-12.d0*p7)-24.d0*p3*p5*p4)/(r(i)**3)
           F_gamma_LF = F_gamma_NF/r(i)
           ! Termes epsilon
           F_epsilon_FF = -6.d0*p4*p5*p6*(p6*p3-p7)/r(i)
           F_epsilon_MF = 3.d0*p4*p5*(2.d0*p3-12.d0*p3*p6**2+6.d0*p6*p7)/(r(i)**2)
           F_epsilon_NF = 3.d0*p4*p5*(6.d0*p3-30.d0*p3*p6**2+12.d0*p6*p7)/(r(i)**3)
           F_epsilon_LF = F_epsilon_NF/r(i)
           ! Termes delta
           F_delta_FF = -F_epsilon_FF
           F_delta_MF = -F_epsilon_MF
           F_delta_NF = -F_epsilon_NF
           F_delta_LF = F_delta_NF/r(i)

           ! Splitting total
           IS(ib) = IS(ib)+gamma(i)*F_gamma_FF*herm_FF*A_FF &
                + gamma(i)*F_gamma_MF*herm_MF*A_MF        &
                + gamma(i)*F_gamma_NF*herm_NF*A_NF        &
                + gamma(i)*F_gamma_LF*herm_LF*A_LF        &
                + epsilon(i)*F_epsilon_FF*herm_FF*A_FF    &
                + epsilon(i)*F_epsilon_MF*herm_MF*A_MF    &
                + epsilon(i)*F_epsilon_NF*herm_NF*A_NF    &
                + epsilon(i)*F_epsilon_LF*herm_LF*A_LF    &
                + delta(i)*F_delta_FF*herm_FF*A_FF        &
                + delta(i)*F_delta_MF*herm_MF*A_MF        &
                + delta(i)*F_delta_NF*herm_NF*A_NF        &
                + delta(i)*F_delta_LF*herm_LF*A_LF

           !      print*,ib,IS(ib),i,gamma(i),F_gamma_FF,herm_FF,A_FF

           idepth = int((z(i)-zmin)/dz)+1
           Kz1(idepth) = Kz1(idepth)+F_gamma_FF*herm_FF*A_FF*dxvol
           Kz2(idepth) = Kz2(idepth)+F_gamma_MF*herm_MF*A_MF*dxvol
           Kz3(idepth) = Kz3(idepth)+F_gamma_NF*herm_NF*A_NF*dxvol
           Kz4(idepth) = Kz4(idepth)+F_gamma_LF*herm_LF*A_LF*dxvol
           !      Kz1(idepth) = Kz1(idepth)+F_epsilon_FF*herm_FF*A_FF*dxvol
           !      Kz2(idepth) = Kz2(idepth)+F_epsilon_MF*herm_MF*A_MF*dxvol
           !      Kz3(idepth) = Kz3(idepth)+F_epsilon_NF*herm_NF*A_NF*dxvol
           !      Kz4(idepth) = Kz4(idepth)+F_epsilon_LF*herm_LF*A_LF*dxvol
           !      Kz1(idepth) = Kz1(idepth)+F_delta_FF*herm_FF*A_FF*dxvol
           !      Kz2(idepth) = Kz2(idepth)+F_delta_MF*herm_MF*A_MF*dxvol
           !      Kz3(idepth) = Kz3(idepth)+F_delta_NF*herm_NF*A_NF*dxvol
           !      Kz4(idepth) = Kz4(idepth)+F_delta_LF*herm_LF*A_LF*dxvol
           !      write(300) x(i),y(i),z(i),F_gamma_FF*herm_FF*A_FF*dxvol
           if(print_kernel)then
              if(ib.eq.1)then
                 write(300) x(i),y(i),z(i),&
                      F_gamma_FF*herm_FF*A_FF*dxvol
              endif
           endif
        endif
     enddo

     !write(200,112)slon,alpha/deg2rad,sin(inci)*111.1949266d0/beta,tau,IS(ib)*dxvol,error
     write(200,112) alpha/deg2rad,IS(ib)*dxvol
112  format(f5.1,1x,f6.3)

  enddo
  !
  ! profile
  !
  do i = 1,Nz
     write(400,111)zmin+float(i-1)*pasz+pasz/2.,&
          Kz1(i),Kz2(i),Kz3(i),Kz4(i),Kz1(i)+Kz2(i)+Kz3(i)+Kz4(i)
  end do
111 format(f5.1,5(1x,f8.4))
113 format(8(1x,f12.4))


  close(200)
  close(300)
  close(400)

end program splitting_total

