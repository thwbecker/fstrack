c
c aniosotropic reflectivity package from Vera Schulte-Pelkum in March 2005
c based on code by Peter Shearer and various others, see individual READMEs 
c and source code comments for details and description
c  
c minor modifications by Thorsten Becker
c
c $Id: utilcijkl.f,v 1.4 2010/03/04 21:51:31 becker Exp $
c
c
c

c subroutine GETMODEL prompts and obtains an anisotropy model
      subroutine GETMODEL(rho,cr)
      implicit double precision (a-h,o-z)
      dimension c(3,3,3,3),cr(3,3,3,3)
      real*8 lamda,mu
c
*      print *,' '
*      print *,'1)  Keith & Crampin isotropic'
*      print *,'2)  Keith & Crampin olivine'
*      print *,'3)  Keith & Crampin upper mantle'
*      print *,'4)  Stephen JGR 1985'
*      print *,'5)  Any hexagonal model'
*      print *,'6)  Any isotropic model'
*      print *,'7)  Almost isotropic test'
*      print *,'8)  Simple % anisotropy hexagonal'
*      print *,'9)  Backus coefficients'
*      print *,'10) Hudson crack model'
*      print *,'11) Standard crack models'
*      print *,'12) Iron (cubic)'
*      print *,'13) alpha-quatrz (trigonal)'
*      print *,'14) titanium'
*      print *,' '
*      print *,'Enter desired model: '
10     continue
      read *,menu
      if (menu.lt.1.or.menu.gt.14) go to 10
      call VERA_ZEROOUT(c)
      if (menu.eq.1) then
         print *,'Keith and Crampin isotropic model'
         vp=10.0d0
         vs=5.77d0
         rho=1.0d0                     !assume normalized density
         lamda=rho*(vp**2-2.0d0*vs**2)
         mu=rho*vs**2
         c(1,1,1,1)=lamda+2.0d0*mu
         c(2,2,2,2)=lamda+2.0d0*mu
         c(3,3,3,3)=lamda+2.0d0*mu
         c(1,1,2,2)=lamda
         c(1,1,3,3)=lamda
         c(2,2,3,3)=lamda
         c(1,2,1,2)=mu
         c(1,3,1,3)=mu
         c(2,3,2,3)=mu
         rho=3.6
      else if (menu.eq.2) then
         print *,'Keith and Crampin olivine model'
         rho=3.324d0
         c(1,1,1,1)=324.0d0/rho
         c(1,1,2,2)=59.0d0/rho
         c(1,2,1,2)=79.0d0/rho
         c(2,2,2,2)=198.0d0/rho
         c(2,2,3,3)=78.0d0/rho
         c(2,3,2,3)=66.7d0/rho
         c(3,3,3,3)=249.0d0/rho
         c(3,3,1,1)=79.0d0/rho
         c(1,3,1,3)=81.0d0/rho
      else if (menu.eq.3) then
         print *,'Keith & Crampin upper mantle model'
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
      else if (menu.eq.4) then
         print *,'Stephen JGR 1985 model'
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
      else if (menu.eq.5) then
         print *,'General hexagonal model'
*         print *,'Enter c1111,c2222,c1212,c1122,c1313'
         read *,c(1,1,1,1),c(2,2,2,2),c(1,2,1,2),c(1,1,2,2),c(1,3,1,3)
         c(3,3,3,3)=c(1,1,1,1)               !additional dependent terms for
         c(3,3,2,2)=c(1,1,2,2)               !hexagonal symmetry, see Musgrave,
         c(3,2,3,2)=c(1,2,1,2)               !Crystal Acoustics, for details
         c(1,1,3,3)=c(1,1,1,1)-2.0d0*c(1,3,1,3)
*         print *,'Enter density'
         read *,rho
      else if (menu.eq.6) then
         print *,'General isotropic model'
*         print *,'Enter Vp, Vs, density'
         read *,vp,vs,den
         rho=1.0d0                     !assume normalized density
         lamda=rho*(vp**2-2.0d0*vs**2)
         mu=rho*vs**2
         c(1,1,1,1)=lamda+2.0d0*mu
         c(2,2,2,2)=lamda+2.0d0*mu
         c(3,3,3,3)=lamda+2.0d0*mu
         c(1,1,2,2)=lamda
         c(1,1,3,3)=lamda
         c(2,2,3,3)=lamda
         c(1,2,1,2)=mu
         c(1,3,1,3)=mu
         c(2,3,2,3)=mu
         rho=den
      else if (menu.eq.7) then
         print *,'Almost isotropic test model'
*         print *,'Enter Vp, Vs, density'
         read *,vp,vs,den
         rho=den
*         print *,'Enter 2theta term (example: .1)'
         read*,err
         a=vp**2
         b=err
         cc=0.0d0
         d=vs**2
         e=err
         c(1,1,1,1)=a+b+cc
         c(2,2,2,2)=a-b+cc
         c(1,3,1,3)=d+e
         c(2,3,2,3)=d-e
         c(3,3,3,3)=c(2,2,2,2)
         c(1,2,1,2)=c(1,3,1,3)
         c(1,1,3,3)=a-3.0d0*cc-2.0d0*(d+e)
         c(1,1,2,2)=c(1,1,3,3)
         c(2,2,3,3)=a-b+cc-2.0d0*(d-e)
      else if (menu.eq.8) then
         rho=1.                  !assume normalized density
         print *,'Simple % anisotropy hexagonal model'
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
         vsavg=vpavg/sqrt(3.0d0)    !assume 0.25 Poisson's ratio
         vsmax=vsavg+(vsavg*fract)/2.0d0
         vsmin=vsavg-(vsavg*fract)/2.0d0
         d=0.5d0*(vsmax**2+vsmin**2)
         e=vsmax**2-d
         if (isym.eq.1) then
            b=-b
            cc=-cc
            e=-e
         end if
*         print *,'a,b,c,d,e=',a,b,cc,d,e
         call CAL_TENSOR(c,a,b,cc,d,e)
      else if (menu.eq.9) then
*         print *,'Enter P terms:  A, B (2-theta), C (4-theta)'
         read *,a,b,cc
*         print *,'Enter SV terms:  D, E (2-theta)'
         read *,d,e
         call CAL_TENSOR(c,a,b,cc,d,e)
      else if (menu.eq.10) then
*         print *,'Enter pvel,svel,rho for isotropic matrix'
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
      else if (menu.eq.11) then
         vp=4.5d0
         poi=.27d0
         vs=vp/sqrt((2.0d0-2.0d0*poi)/(1.0d0-2.0d0*poi))
         rho=2.8d0
         vp1=1.5d0
         vs1=0.0d0
         rho1=1.0d0
*         print *,'Enter crack aspect ratio d'
         read *,d
*         print *,'Enter crack density parameter epsilon'
         read *,epsilon
         call HUDSON(vp,vs,rho,vp1,vs1,rho1,d,epsilon,c,vpiso,vsiso)
      else if (menu.eq.12) then
         rho=7.86d0
         c(1,1,1,1)=233.0d0/rho              !values from Musgrave
         c(2,2,2,2)=233.0d0/rho              !cubic symmetry
         c(3,3,3,3)=233.0d0/rho
         c(1,1,2,2)=139.2d0/rho
         c(1,1,3,3)=139.2d0/rho
         c(2,2,3,3)=139.2d0/rho
         c(2,3,2,3)=116.2d0/rho
         c(1,3,1,3)=116.2d0/rho
         c(1,2,1,2)=116.2d0/rho
      else if (menu.eq.13) then
         rho=2.65d0
         c(1,1,1,1)=86.74d0/rho                !values from Bechmann 1958
         c(3,3,3,3)=107.2d0/rho
         c(1,1,2,2)=6.99d0/rho
         c(1,1,3,3)=11.91d0/rho
         c(2,3,2,3)=57.94d0/rho
         c(1,1,2,3)=-17.91d0/rho
         c(2,2,2,2)=c(1,1,1,1)          !trigonal symmetry
         c(1,3,1,3)=c(2,3,2,3)
         c(1,2,1,2)=0.5d0*(c(1,1,1,1)-c(1,1,2,2))
         c(2,2,3,3)=c(1,1,3,3)
         c(2,2,2,3)=-c(1,1,2,3)
         c(1,3,1,2)=c(1,1,2,3)
      else if (menu.eq.14) then    !titanium values from Wenk et al. 1988
         c(1,1,1,1)=160.0d0
         c(2,2,2,2)=c(1,1,1,1)     !hexagonal symmetry (001) axis
         c(3,3,3,3)=181.0d0
         c(2,3,2,3)=46.5d0
         c(1,3,1,3)=c(2,3,2,3)
         c(1,1,2,2)=90.0d0
         c(1,1,3,3)=66.0d0
         c(2,2,3,3)=c(1,1,3,3)
         c(1,2,1,2)=0.5d0*(c(1,1,1,1)-c(1,1,2,2))
      end if
      call vera_fillout(c)
      iaxis=1
      theta=0.0d0
*      print *,'Coordinates can be rotated clockwise'
*      print *,'Enter desired rotation axis, angle'
      read *,iaxis,theta
      call vera_rotate(c,cr,iaxis,theta)
      return
      end
c
c
c
c
c subroutine GETCIJKL returns the elastic tensor and its depth
c derivative for a given depth.  Note that Cijkl is scaled by
c a simple factor so that its elements stay in proportion to
c each other.  The scaling is such that the P-wave velocity
c of the isotropic part of the tensor increases linearly with depth.
c
c    Inputs:    c(3,3,3,3)  =  elastic tensor at zero depth
c               pvel0       =  Isotropic P-wave velocity at zero depth
c               pgrad       =  P-wave velocity gradient (km/s/km)
c               z           =  depth
c    Returns:   a(3,3,3,3)  =  elastic tensor at depth z
c             dadz(3,3,3,3) =  derivative of a with respect to depth
c
      subroutine GETCIJKL(c,pvel0,pgrad,z,a,dadz)
      implicit double precision (a-h,o-z)
      real*8 c(3,3,3,3),a(3,3,3,3),dadz(3,3,3,3)
c
      afact=(pvel0+pgrad*z)**2/pvel0**2
      dafact=2.0d0*pgrad*(pvel0+pgrad*z)/pvel0**2
c
      do 40 i=1,3
      do 30 j=1,3
      do 20 k=1,3
      do 10 l=1,3
         a(i,j,k,l)=c(i,j,k,l)*afact
         dadz(i,j,k,l)=c(i,j,k,l)*dafact
10    continue
20    continue
30    continue
40    continue
c
      return
      end
c
c
c
c
c HUDSON calculates anisotropic elastic parameters for a material
c containing aligned cracks and also the isotropic velocities for
c randomly distributed cracks.  Based on papers by:
c
c           Hudson, J.A.  Overall properties of a cracked solid.
c        Math. Proc. Camb. Phil. Soc. (1980?), 88, 371.
c
c           Crampin, S.  Effective anisotropic elastic constants
c        for wave propagation through cracked solids.  Geophys.
c        J. R. astr. Soc. (1984),76,135-145.
c
c The routine assumes hexagonal symmetry with a horizontal symmetry
c axis corresponding to (100) which is the convention used in the
c Crampin paper.
c
c    Inputs:   vp          =  host P-wave velocity
c              vs          =  host S-wave velocity
c              rho         =  host density
c              vp1         =  crack material P-wave velocity
c              vs1         =  crack material S-wave velocity
c              rho1        =  crack material density
c              d           =  crack aspect ratio (<<1)
c              epsilon     =  crack density parameter (<<1)
c    Returns:  c(3,3,3,3)  =  elastic tensor of composite material
c              vpiso       =  P-wave velocity for random cracks
c              vsiso       =  S-wave velocity for random cracks
c
c    Requires:  VERA_ZEROOUT,VERA_FILLOUT
c
      subroutine HUDSON(vp,vs,rho,vp1,vs1,rho1,d,epsilon,c,vpiso,vsiso)
      implicit double precision (a-h,o-z)
      real*8 c(3,3,3,3),c0(3,3,3,3),c1(3,3,3,3),c2(3,3,3,3)
      real*8 lamda,mu,lamda1,mu1,k1,lamda2,mu2,lamda3,mu3
      pi=3.1415926d0
c
      lamda=rho*(vp**2-2.0d0*vs**2)         !Lame parameters
      mu=rho*vs**2                          !of isotropic matrix
      lamda1=rho1*(vp1**2-2.0d0*vs1**2)     !Lame parameters
      mu1=rho1*vs1**2                       !of crack material
      k1=lamda1+(2.0d0/3.0d0)*mu1           !bulk modulus of crack material
c
c Calculate isotropic tensor
      call VERA_ZEROOUT(c0)
      c0(1,1,1,1)=lamda+2.0d0*mu
      c0(2,2,2,2)=lamda+2.0d0*mu
      c0(3,3,3,3)=lamda+2.0d0*mu
      c0(1,1,2,2)=lamda
      c0(1,1,3,3)=lamda
      c0(2,2,3,3)=lamda
      c0(1,2,1,2)=mu
      c0(1,3,1,3)=mu
      c0(2,3,2,3)=mu
      call VERA_FILLOUT(c0)
c
c Calculate first order perturbation
      call VERA_ZEROOUT(c1)
      fk=((k1+(4.0d0/3.0d0)*mu1)/(pi*d*mu))*((lamda+2.0d0*mu)/
     &     (lamda+mu))
      fm=(4.0d0*mu1/(pi*d*mu))*((lamda+2.0d0*mu)/(3.0d0*lamda+4.0d0*mu))
      u11=(4.0d0/3.0d0)*((lamda+2.0d0*mu)/(lamda+mu))/(1.0d0+fk)
      u33=(16.0d0/3.0d0)*((lamda+2.0d0*mu)/(3.0d0*lamda+4.0d0*mu))/
     &     (1.0d0+fm)
      scr1=-u11*epsilon/mu
      scr3=-u33*epsilon/mu
      c1(1,1,1,1)=scr1*(lamda+2.0d0*mu)**2
      c1(2,2,2,2)=scr1*(lamda**2)
      c1(3,3,3,3)=c1(2,2,2,2)
      c1(1,1,2,2)=scr1*lamda*(lamda+2.0d0*mu)
      c1(1,1,3,3)=c1(1,1,2,2)
      c1(2,2,3,3)=c1(2,2,2,2)
      c1(1,3,1,3)=scr3*mu**2
      c1(1,2,1,2)=c1(1,3,1,3)
      call VERA_FILLOUT(c1)
c
c Calculate second order perturbation
      call VERA_ZEROOUT(c2)
      scr1=u11**2*epsilon**2/15.0d0
      scr3=u33**2*epsilon**2/15.0d0
      fq=15.0d0*(lamda/mu)**2+28.0d0*(lamda/mu)+28.0d0
      fx=2.0d0*mu*(3.0d0*lamda+8.0d0*mu)/(lamda+2.0d0*mu)
      c2(1,1,1,1)=scr1*(lamda+2.0d0*mu)*fq
      c2(2,2,2,2)=scr1*lamda**2*fq/(lamda+2.0d0*mu)
      c2(3,3,3,3)=c2(2,2,2,2)
      c2(1,1,2,2)=scr1*lamda*fq
      c2(1,1,3,3)=c2(1,1,2,2)
      c2(2,2,3,3)=c2(2,2,2,2)
      c2(1,3,1,3)=scr3*fx
      c2(1,2,1,2)=c2(1,3,1,3)
c      call VERA_FILLOUT(c1)  !***this error was in program until 19 Feb. 87***
      call VERA_FILLOUT(c2)
c
c Sum c0,c1,c2
      do 24 i=1,3
      do 23 j=1,3
      do 22 k=1,3
      do 21 l=1,3
         c(i,j,k,l)=(c0(i,j,k,l)+c1(i,j,k,l)+c2(i,j,k,l))/rho
21    continue
22    continue
23    continue
24    continue
c
c Now calculate isotropic velocities for random case in
c Hudson equations 59, 61, 62.  Note that we overwrite the
c old values for mu1 and lamda1 (sloppy, sloppy!)
c Note also that u11 and u33 are switched in these equations
c compared to the Hudson equations to account for the switch
c in symmetry axis to match the Crampin paper.
      mu1=-mu*epsilon*2.0d0*(3.0d0*u33+2.*u11)/15.0d0
      scr=-epsilon*(3.0d0*lamda+2.0d0*mu)*u11/(3.0d0*mu)
c note typo in Hudson Eqn. (59): extra lamda term!
      lamda1=(1.0d0/3.0d0)*(scr*(3.0d0*lamda+2.0d0*mu)-2.0d0*mu1)
      mu2=mu*(2.0d0/15.0d0)*(mu1/mu)**2*(3.0d0*lamda+8.0d0*mu)/
     &     (lamda+2.0d0*mu)
      scr=(3.0d0*lamda1+2.0d0*mu1)**2/(3.0d0*(lamda+2.0d0*mu)*
     &     (3.0d0*lamda+2.0d0*mu))
      lamda2=(1.0d0/3.0d0)*(scr*(3.0d0*lamda+2.0d0*mu)-2.0d0*mu2)
      lamda3=lamda+lamda1+lamda2
      mu3=mu+mu1+mu2
      poros=(4.0d0/3.0d0)*pi*epsilon*d
      rho3=rho*(1.0d0-poros)+rho1*poros
      vpiso=sqrt((lamda3+2.0d0*mu3)/rho3)
      vsiso=sqrt(mu3/rho3)
c
      return
      end
c
c
c
c
c VERA_ZEROOUT is designed to fill Cijkl with zeros
      SUBROUTINE VERA_ZEROOUT(c)
      implicit double precision (a-h,o-z)
      dimension c(3,3,3,3)
      do 80 i=1,3
      do 70 j=1,3
      do 60 k=1,3
      do 50 l=1,3
         c(i,j,k,l)=0.0d0
50    continue
60    continue
70    continue
80    continue
      return
      end
c
c
c
c
c SWITCH is designed to switch pairs of indices in Cijkl
      SUBROUTINE SWITCH(c,i1,i2)
      implicit double precision (a-h,o-z)
      dimension c(3,3,3,3),d(3,3,3,3)
      do 40 i=1,3
      do 30 j=1,3
      do 20 k=1,3
      do 10 l=1,3
         ii=i
         jj=j
         kk=k
         ll=l
         if (i.eq.i1) ii=i2
         if (i.eq.i2) ii=i1
         if (j.eq.i1) jj=i2
         if (j.eq.i2) jj=i1
         if (k.eq.i1) kk=i2
         if (k.eq.i2) kk=i1
         if (l.eq.i1) ll=i2
         if (l.eq.i2) ll=i1
         d(ii,jj,kk,ll)=c(i,j,k,l)
10    continue
20    continue
30    continue
40    continue
c
      do 80 i=1,3
      do 70 j=1,3
      do 60 k=1,3
      do 50 l=1,3
      c(i,j,k,l)=d(i,j,k,l)
50    continue
60    continue
70    continue
80    continue
c
      return
      end
c
c
c
c
c ROTATE rotates an elastic tensor about a specified angle, 
c angle in degrees
c
c
      SUBROUTINE vera_rotate(c,cc,iaxis,theta)
      implicit double precision (a-h,o-z)
      dimension c(3,3,3,3),cc(3,3,3,3),r(3,3)
      pifac=57.295779513082320876798154814105d0
c
      thetar = theta/pifac
      do i=1,3
         do j=1,3
            if (i.eq.iaxis.and.j.eq.iaxis) then
               r(i,j)=1.d0
            else if (i.eq.iaxis.or.j.eq.iaxis) then
               r(i,j)=0.d0
            else if (i.eq.j) then
               r(i,j)=cos(thetar)
            else if (i.gt.j) then
               r(i,j)=sin(thetar)
            else
               r(i,j)=-sin(thetar)
            end if
         enddo
      enddo
      
c
      do ii=1,3
         do jj=1,3
            do kk=1,3
               do ll=1,3
c     
                  cc(ii,jj,kk,ll)=0.0d0
                  do i=1,3
                     do  j=1,3
                        do k=1,3
                           do l=1,3
                              cc(ii,jj,kk,ll)=cc(ii,jj,kk,ll)+
     &                             c(i,j,k,l)*r(ii,i)*r(jj,j)*
     &                             r(kk,k)*r(ll,l)
                              
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

c
900   return
      end
c
c
c
c
c VERA_FILLOUT is designed to fill out Cijkl to its proper 81 values by using
c appropriate symmetry relationships for an elastic tensor.  Thus, by
c using this subroutine, the user need only specify the non-redundant values
c in Cijkl.  This routine makes no symmetry assumptions other than those
c that are always true for an elastic tensor (i.e. 21 independent values).
      SUBROUTINE VERA_FILLOUT(c)
      implicit double precision (a-h,o-z)
      dimension c(3,3,3,3)
      do 80 i=1,3              !use symmetry relationships to fill
      do 70 j=1,3              !in Cijkl where appropriate
      do 60 k=1,3
      do 50 l=1,3
         if (c(i,j,k,l).eq.0.) c(i,j,k,l)=c(j,i,k,l)
         if (c(i,j,k,l).eq.0.) c(i,j,k,l)=c(i,j,l,k)
         if (c(i,j,k,l).eq.0.) c(i,j,k,l)=c(j,i,l,k)
         if (c(i,j,k,l).eq.0.) c(i,j,k,l)=c(k,l,i,j)
         if (c(i,j,k,l).eq.0.) c(i,j,k,l)=c(l,k,i,j)
         if (c(i,j,k,l).eq.0.) c(i,j,k,l)=c(k,l,j,i)
         if (c(i,j,k,l).eq.0.) c(i,j,k,l)=c(l,k,j,i)
50    continue
60    continue
70    continue
80    continue

      return
      end
c
c
c
c
c SHOW66 displays an elastic tensor in a convenient symmetric 6x6 form.
c See Musgrave for details.
c Key:  1=11, 2=22, 3=33, 4=23=32, 5=31=13, 6=12=21
      SUBROUTINE vera_show66(c)
      implicit double precision (a-h,o-z)
      dimension c(3,3,3,3)
      print 10,c(1,1,1,1),c(1,1,2,2),c(1,1,3,3),c(1,1,2,3),c(1,1,3,1)
     .        ,c(1,1,1,2)
10    format (6f8.3)
      print 20,c(2,2,2,2),c(2,2,3,3),c(2,2,2,3),c(2,2,3,1),c(2,2,1,2)
20    format (8x,5f8.3)
      print 30,c(3,3,3,3),c(3,3,2,3),c(3,3,3,1),c(3,3,1,2)
30    format (2(8x),4f8.3)
      print 40,c(2,3,2,3),c(2,3,3,1),c(2,3,1,2)
40    format (3(8x),3f8.3)
      print 50,c(3,1,3,1),c(3,1,1,2)
50    format (4(8x),2f8.3)
      print 60,c(1,2,1,2)
60    format (5(8x),f8.3)
      return
      end
c
c
c
c
c HOWWEAK splits an elastic tensor into an isotropic tensor and an aelotropic
c tensor, calculates the Lame parameters for the isotropic tensor, and
c computes ae, a measure of the relative strength of the anisotropy.
c For details, see Backus, 1982, JGR v. 87, p. 4641-4644.
      SUBROUTINE HOWWEAK(c,ciso,cael,lamda,mu,ae)
      implicit double precision (a-h,o-z)
      real*8 c(3,3,3,3),ciso(3,3,3,3),cael(3,3,3,3)
      real*8 s(3,3,3,3),n(3,3,3,3),del(3,3)
      real*8 lamda, mu
c
      do 20 i=1,3                !set up delta function
         do 10 j=1,3
            del(i,j)=0.0d0
            if (i.eq.j) del(i,j)=1.0d0
10       continue
20    continue
c
      cdots=0.              !This section defines s and n, an orthonormal basis
      cdotn=0.              !in the two-dimensional space of isotropic elastic
      do 60 i=1,3           !tensors.  It also computes the inner product of
      do 50 j=1,3           !these tensors with our elastic tensor.
      do 40 k=1,3
      do 30 l=1,3
      s(i,j,k,l)=del(i,j)*del(k,l)+del(i,k)*del(j,l)+del(i,l)*del(j,k)
      s(i,j,k,l)=s(i,j,k,l)/sqrt(45.0d0)
      n(i,j,k,l)=2.0d0*del(i,j)*del(k,l)-del(i,k)*del(j,l)-
     &     del(i,l)*del(j,k)
      n(i,j,k,l)=n(i,j,k,l)/6.0d0
      cdots=cdots+c(i,j,k,l)*s(i,j,k,l)
      cdotn=cdotn+c(i,j,k,l)*n(i,j,k,l)
30    continue
40    continue
50    continue
60    continue
c
c      print *, 's tensor'
c      call SHOW66(s)
c      print *, 'n tensor'
c      call SHOW66(n)
c
      cisomag=0.0d0            !This section computes ciso as the orthogonal projec
      caelmag=0.0d0            !of c onto the space defined by s and n.  cael is th
      do 100 i=1,3          !the difference between c and ciso.  We also compute
      do 90 j=1,3           !magnitude of ciso and ceal in this section.
      do 80 k=1,3
      do 70 l=1,3
      ciso(i,j,k,l)=cdots*s(i,j,k,l)+cdotn*n(i,j,k,l)
      cael(i,j,k,l)=c(i,j,k,l)-ciso(i,j,k,l)
      cisomag=cisomag+ciso(i,j,k,l)**2
      caelmag=caelmag+cael(i,j,k,l)**2
70    continue
80    continue
90    continue
100   continue
      cisomag=sqrt(cisomag)
      caelmag=sqrt(caelmag)
c      print *, 'cisomag= ',cisomag
c      print *, 'caelmag= ',caelmag
c
      lamda=ciso(1,1,2,2)
      mu=ciso(1,2,1,2)
      ae=caelmag/(lamda+mu)
c
      esum1=0.0d0            !This is a test to make sure cael
      esum2=0.0d0            !is aelotropic.  Both sums should
      do 120 i=1,3        !be zero.  For details, see eq. 6
      do 110 j=1,3        !in the Backus paper.
         esum1=esum1+cael(i,i,j,j)
         esum2=esum2+cael(i,j,i,j)
110   continue
120   continue
c      print *, 'cael(iijj)= ',esum1
c      print *, 'cael(ijij)= ',esum2
c
      return
      end
c
c
c
c
c CAL_ABCDEFG calculates the Backus velocity coefficients
      SUBROUTINE CAL_ABCDEFG(x,a,b,c,d,e,f,g)
      implicit double precision (a-h,o-z)
      dimension x(3,3,3,3)
      a=(3.0d0*x(1,1,1,1)+3.*x(2,2,2,2)+2.0d0*x(1,1,2,2)+4.0d0*
     &     x(1,2,1,2))/8.0d0
      b=(x(1,1,1,1)-x(2,2,2,2))/2.0d0
      c=(x(1,1,1,1)+x(2,2,2,2)-2.0d0*x(1,1,2,2)-4.0d0*x(1,2,1,2))/8.0d0
      d=(x(1,3,1,3)+x(2,3,2,3))/2.0d0
      e=(x(1,3,1,3)-x(2,3,2,3))/2.0d0
      f=c+d+e
      g=-c
      return
      end
c
c
c
c
c subroutine SHOW_ABCDEFG prints useful information about
c the Backus coefficients.
      SUBROUTINE SHOW_ABCDEFG(a,b,c,d,e,f,g)
      implicit double precision (a-h,o-z)
*     print *,' '
*     print *,'ANISOTROPIC VELOCITY COEFFICIENTS'
*     print *,'for HEXAGONAL SYMMETRY WITH 100 AXIS'
*     print *,' P-wave terms A= ',a
*     print *,' 	     B= ',b
*     print *,' 	     C= ',c
*     print *,'SV-wave terms D= ',d
*     print *,' 	     E= ',e
*     print *,'SH-wave terms F= ',f
*     print *,' 	     G= ',g
      p0=sqrt(a+b+c)
      p45=sqrt(a-c)
      p90=sqrt(a-b+c)
      sv0=sqrt(d+e)
      sv90=sqrt(d-e)
      sh0=sqrt(f+g)
      sh45=sqrt(f-g)
*     print *,'Vp(0)=	',p0
*     print *,'Vp(45)=  ',p45
*     print *,'Vp(90)=  ',p90
*     print *,'Vsv(0)=  ',sv0
*     print *,'Vsv(90)= ',sv90
*     print *,'Vsh(0)=  ',sh0
*     print *,'Vsh(45)= ',sh45
      return
      end
c
c
c
c
c CAL_TENSOR calculates the elastic tensor from the coefficients
      SUBROUTINE CAL_TENSOR(x,a,b,c,d,e)
      implicit double precision (a-h,o-z)
      dimension x(3,3,3,3)
      call VERA_ZEROOUT(x)
      x(1,1,1,1)=a+b+c
      x(2,2,2,2)=a-b+c
      x(3,3,3,3)=x(2,2,2,2)
      x(1,3,1,3)=d+e
      x(1,2,1,2)=x(1,3,1,3)
      x(2,3,2,3)=d-e
      x(1,1,2,2)=a-3.0d0*c-2.0d0*(d+e)
      x(1,1,3,3)=x(1,1,2,2)
      x(2,2,3,3)=a-b+c-2.0d0*(d-e)
      call VERA_FILLOUT(x)
      return
      end

      subroutine vera_cijkl_assign(source,target)
      implicit none
      double precision source(3,3,3,3),target(3,3,3,3)
      integer i,j,k,l
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  target(i,j,k,l) = source(i,j,k,l)
               enddo
            enddo
         enddo
      enddo
      end 
               
      subroutine vera_get_isym(sym,isym)
      implicit none
      double precision sym(3,3)
      integer isym,i1,i2,i3,i,j
c      
c Find 5 independent constants for transverse isotropy and set up sym(3,3)
c     
      do  i=1,3
         do  j=1,3
            sym(i,j)=0.0d0
         enddo
      enddo
      if (isym.eq.3) then
         i1=1
         i2=2
         i3=3
         sym(1,1)=1.0d0
         sym(2,2)=1.0d0
         sym(3,3)=1.0d0
      else if (isym.eq.2) then
         i1=1
         i2=3
         i3=2
         sym(1,3)=1.0d0
         sym(2,1)=1.0d0
         sym(3,2)=1.0d0
      else if (isym.eq.1) then
         i1=3
         i2=2
         i3=1
         sym(1,2)=1.0d0
         sym(2,3)=1.0d0
         sym(3,1)=1.0d0
      end if
      end

c
c
c
c
c VERA_SYM_ROT rotates sym(3,3) matrix about a specified angle
c (given in degrees)
c
      SUBROUTINE vera_sym_rot(sym,symr,iaxis,theta)
      implicit double precision (a-h,o-z)
      dimension sym(3,3),symr(3,3),r(3,3)

      pif=57.295779513082320875d0
c
      thetar = theta/pif
      do  i=1,3
         do  j=1,3
            if (i.eq.iaxis.and.j.eq.iaxis) then
               r(i,j)=1.0d0
            else if (i.eq.iaxis.or.j.eq.iaxis) then
               r(i,j)=0.0d0
            else if (i.eq.j) then
               r(i,j)=cos(thetar)
            else if (i.gt.j) then
               r(i,j)=sin(thetar)
            else
               r(i,j)=-sin(thetar)
            end if
         enddo
      enddo
c
      do  i=1,3
         do  j=1,3
            symr(i,j)=sym(i,1)*r(j,1)+sym(i,2)*r(j,2)+sym(i,3)*r(j,3)
         enddo
      enddo
c
      return
      end

