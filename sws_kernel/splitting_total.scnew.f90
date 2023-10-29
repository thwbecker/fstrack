
!----------------------------------------------------------------------------
! Ce programme calcule l'intensite du splitting en surface produit par
! anisotropie hexagonale distribuee decrite par une orientation quelconque
! de l'axe de symetrie et par les deux parametres anisotropes
!----------------------------------------------------------------------------

program splitting_total

  implicit none

  ! pour boucle sur les backazimuths
  integer, parameter :: nbtotal = 180

  integer, parameter :: nptmax = 6000000

  integer :: i,ind,ib,npt,Nx,Ny,Nz,idepth,icas

  double precision :: x(nptmax),y(nptmax),z(nptmax)
  double precision :: azi(nptmax),inc(nptmax)
  double precision :: gamma(nptmax),epsilon(nptmax),delta(nptmax)

  double precision :: A_FF,A_MF,A_NF,A_LF
  double precision :: herm_FF,herm_MF,herm_NF,herm_LF
  double precision :: F_gamma_FF,F_gamma_MF,F_gamma_NF,F_gamma_LF
  double precision :: F_epsilon_FF,F_epsilon_MF,F_epsilon_NF,F_epsilon_LF
  double precision :: F_delta_FF,F_delta_MF,F_delta_NF,F_delta_LF

  double precision :: deg2rad,beta,tau,lambda,aa,r,ro,dx,dy,dz
  double precision :: kr,zmin,zmax,pas_alpha
  double precision :: phi,sinth,costh,cosp,sinp,alpha,inci
  double precision :: xsta,ysta,slon,slat,error
  double precision :: p1,p2,p3,p4,p5,p6,p7
  double precision :: k1(3),k2(3),t(3),p(3),s(3)
  double precision :: IS(100),Kz1(100),Kz2(100),Kz3(100),Kz4(100)

  ! pi
  double precision, parameter :: pi = 3.141592653589793d0

  deg2rad = pi/180.d0
  slon = 242.d0
  slat = 34.d0
  error = 0.1d0
  pas_alpha = 360.d0/nbtotal
  !  pas_alpha = 2.d0
  icas = 1

  ! Ouverture des fichiers contenant le modele anisotrope et
  ! l'intensite du splitting en fonction du backazimuth
  open(100,file='aniso.mod',form='unformatted',status='old')
  open(200,file='splitting.dat')
  open(300,file='kernel.dat',form='unformatted',status='unknown')
  open(400,file='profil.dat')

  ! Lecture du modele de reference
  beta = 4.47d0
  tau = 10.d0
  lambda = beta*tau
  inci = 5.0d0*deg2rad

  ! facteurs lies au choix de la forme d'onde
  if (icas.eq.1) then
     A_FF = 1.d0/(2.d0*beta**2*tau*sqrt(2.d0))
     A_MF = 1.d0/(4.d0*pi*beta)
     A_NF = tau/(4.d0*sqrt(2.d0)*pi**2)
     A_LF = beta*tau**2/(8.d0*pi**3)
  elseif (icas.eq.2) then
     A_FF = 1.d0/(120.d0*beta**2*tau*sqrt(2.d0))
     A_MF = 1.d0/(240.d0*pi*beta)
     A_NF = tau/(240.d0*sqrt(2.d0)*pi**2)
     A_LF = beta*tau**2/(480.d0*pi**3)
  endif

  ! Definition de la position de la station
  xsta = 0.d0
  ysta = 0.d0

  ! Lecture du modele anisotrope
  ind = 0
  print *,'Lecture du modele...'

  read(100) dx,dy,dz
  read(100) nx,ny,nz
  read(100) zmin
  read(100) zmax

  write(300) dx,dy,dz
  write(300) nx,ny,nz
  write(300) zmin
  write(300) zmax
  print*,dx,dy,dz
  print*,nx,ny,nz
  print*,zmin
  print*,zmax

  do i = 1,nptmax
     ind = ind+1
     read(100,end=1) x(ind),y(ind),z(ind),azi(ind),inc(ind),  &
          gamma(ind),epsilon(ind),delta(ind)
  enddo

1 continue

  print *,'Fin de lecture...'
  print*,ind-1,' lignes lues ...'

  npt = ind-1

  ! Initialisation de Kz
  do i = 1,Nz
     Kz1(i) = 0.d0
     Kz2(i) = 0.d0
     Kz3(i) = 0.d0
     Kz4(i) = 0.d0
  end do

  ! Boucle sur les backazimuths
  do ib = 1,nbtotal

     print *,ib,' sur un total de ',nbtotal

     alpha = (dble(ib-1)*pas_alpha)*deg2rad
     !  alpha = 45.d0*deg2rad
     k1(1) = -sin(inci)*cos(alpha)
     k1(2) = -sin(inci)*sin(alpha)
     k1(3) = -cos(inci)
     p(1) = cos(inci)*cos(alpha)
     p(2) = cos(inci)*sin(alpha)
     p(3) = -sin(inci)
     t(1) = -sin(alpha)
     t(2) = cos(alpha)
     IS(ib) = 0.d0

     ! Boucle sur tous les points de la grille
     do i = 1,npt
        r = sqrt((x(i)-xsta)**2+(y(i)-ysta)**2+z(i)**2)
        ro = sqrt((x(i)-xsta)**2+(y(i)-ysta)**2)
        sinth = z(i)/r
        costh = ro/r
        if (abs(x(i)-xsta) > 0.d0) then
           !      phi = atan2((y(i)-ysta),(x(i)-xsta))
           phi = atan2((x(i)-xsta),(y(i)-ysta))
        else
           phi = 0.d0
        endif
        cosp = cos(phi)
        sinp = sin(phi)
        !    kr = abs(k1(1)*(x(i)-xsta)+k1(2)*(y(i)-ysta)+k1(3)*z(i))
        kr = abs(k1(1)*(y(i)-ysta)+k1(2)*(x(i)-xsta)+k1(3)*z(i))
        aa = sqrt(2.d0)*(r-kr)*pi/lambda
        if (aa < 2.d0*pi) then
           s(1) = cos(inc(i)*deg2rad)*cos(azi(i)*deg2rad)
           s(2) = cos(inc(i)*deg2rad)*sin(azi(i)*deg2rad)
           s(3) = sin(inc(i)*deg2rad)
           k2(1) = -costh*cosp
           k2(2) = -costh*sinp
           k2(3) = -sinth
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
           F_gamma_FF = (p4*p1+p5*p2)*(4.d0*p3*p6-2.d0*p7)/r
           !      F_gamma_MF = ((p4*p1+p5*p2)*(24.d0*p3*p6-6.d0*p7)-p3*(8.d0*p5*p4+12.d0*p2*p1))/r**2
           !      F_gamma_MF = ((p4*p1+p5*p2)*(24.d0*p3*p6-6.d0*p7)-p3*(8.d0*p5*p4+24.d0*p2*p1))/r**2
           F_gamma_MF = ((p4*p1+p5*p2)*(24.d0*p3*p6-6.d0*p7)-8.d0*p3*p5*p4)/r**2
           !      F_gamma_NF = ((p4*p1+p5*p2)*(60.d0*p3*p6-12.d0*p7)-p3*(24.d0*p5*p4+30.d0*p2*p1))/r**3
           !      F_gamma_NF = ((p4*p1+p5*p2)*(60.d0*p3*p6-12.d0*p7)-p3*(24.d0*p5*p4+60.d0*p2*p1))/r**3
           F_gamma_NF = ((p4*p1+p5*p2)*(60.d0*p3*p6-12.d0*p7)-24.d0*p3*p5*p4)/r**3
           F_gamma_LF = F_gamma_NF/r
           ! Termes epsilon
           F_epsilon_FF = -6.d0*p4*p5*p6*(p6*p3-p7)/r
           F_epsilon_MF = 3.d0*p4*p5*(2.d0*p3-12.d0*p3*p6**2+6.d0*p6*p7)/r**2
           F_epsilon_NF = 3.d0*p4*p5*(6.d0*p3-30.d0*p3*p6**2+12.d0*p6*p7)/r**3
           F_epsilon_LF = F_epsilon_NF/r
           ! Termes delta
           F_delta_FF = -F_epsilon_FF
           F_delta_MF = -F_epsilon_MF
           F_delta_NF = -F_epsilon_NF
           F_delta_LF = F_delta_NF/r

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
           Kz1(idepth) = Kz1(idepth)+F_gamma_FF*herm_FF*A_FF*dx*dy*dz
           Kz2(idepth) = Kz2(idepth)+F_gamma_MF*herm_MF*A_MF*dx*dy*dz
           Kz3(idepth) = Kz3(idepth)+F_gamma_NF*herm_NF*A_NF*dx*dy*dz
           Kz4(idepth) = Kz4(idepth)+F_gamma_LF*herm_LF*A_LF*dx*dy*dz
           !      Kz1(idepth) = Kz1(idepth)+F_epsilon_FF*herm_FF*A_FF*dx*dy*dz
           !      Kz2(idepth) = Kz2(idepth)+F_epsilon_MF*herm_MF*A_MF*dx*dy*dz
           !      Kz3(idepth) = Kz3(idepth)+F_epsilon_NF*herm_NF*A_NF*dx*dy*dz
           !      Kz4(idepth) = Kz4(idepth)+F_epsilon_LF*herm_LF*A_LF*dx*dy*dz
           !      Kz1(idepth) = Kz1(idepth)+F_delta_FF*herm_FF*A_FF*dx*dy*dz
           !      Kz2(idepth) = Kz2(idepth)+F_delta_MF*herm_MF*A_MF*dx*dy*dz
           !      Kz3(idepth) = Kz3(idepth)+F_delta_NF*herm_NF*A_NF*dx*dy*dz
           !      Kz4(idepth) = Kz4(idepth)+F_delta_LF*herm_LF*A_LF*dx*dy*dz
           !      write(300) x(i),y(i),z(i),F_gamma_FF*herm_FF*A_FF*dx*dy*dz

        endif
     enddo

     !  write(200) slat,slon,alpha/deg2rad,sin(inci)*111.1949266d0/beta,tau,IS(ib)*dx*dy*dz,error
     write(200,112) alpha/deg2rad,IS(ib)*dx*dy*dz
112  format(f5.1,1x,f6.3)
     !  print*,dx,dy,dz
     !  print*,alpha/deg2rad,IS(ib)*dx*dy*dz

  enddo

  do i = 1,Nz
     write(400,111)zmin+float(i-1)*dz+dz/2.,Kz1(i),Kz2(i),Kz3(i),Kz4(i),Kz1(i)+Kz2(i)+Kz3(i)+Kz4(i)
  end do
111 format(f5.1,5(1x,f8.4))

  close(100)
  close(200)
  close(300)
  close(400)

end program splitting_total

