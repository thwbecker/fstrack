c
c     calculate time derivatives of work array (output: workder) at time t
c     given the velocity arrays vr,vtheta,vphi (dimensions: nr*ntheta*nphi*nvtimes)
c     specified at times vtimes(nvtimes) and rlevels(nr)
c     
c     the time derivative of work(1...3), the location x, is the velocity at x
c     in a spherical coordinate system
c
c     
c
c
c
c     the time derivative of work(4..12), the strain F, is the velocity gradient
c            matrix, G,  dotted with the strain, ie. G . F
c
c     WARNING: when nwork<=3, the routine will not determine
c     the time derivate fo the strain but only the velocities at location x, ie.
c     the time derivatives of work(1...3)
c
c     ifremovestrain: flag to possibly remove the symmetric strain and only 
c     rotate F, for this we also need vhdivmax, the max horizontal divergence
c     used to determine the plate `age'
c
c     
c     returns VGM, the velocity gradient matrix, in FORTRAN style sorting
c     (this is different from the workder[4-12] derivative part if work[4-12] 
c     is the unity matrix. in this case, workder[4-12] is the VGM in C style 
c     sorting.)
c
c
      subroutine fse_derivs(t,work,workder,nwork,
     &     vr,vtheta,vphi,nr,ntheta,nphi,dtheta,dphi,rlevels,
     &     vtimes,nvtimes,ifremovestrain,vhdivmax,vgm) 
c_______________________________________________________________________
c
c     Note:  passed muti-dimensional arrays are in c-order
c
c_______________________________________________________________________

#include "fortran_precision.h"
      dimension vr(*), vtheta(*), vphi(*), rlevels(nr),
     &     vtimes(*),vhdivmax(*),vgm(3,3)
      dimension work(nwork),workder(nwork)
      dimension grid(0:3,3),weight(0:3,0:1,3)
      dimension igrid(0:3,3),ider(0:3), dvr(0:3),dvtheta(0:3),
     &     dvphi(0:3)

      logical icalc_vgm

      save ntnp,nrntnp,init
      data init /0/
      
      if(init.eq.0)then
         ntnp = ntheta * nphi
         nrntnp = ntnp * nr
c
c     do some checks if called for the first time
c     
         if(nr.lt.4)then
            print *,'ellderivs: error: need at least four depth levels'
            stop
         endif
         do i=2,nr
            if(rlevels(i).le.rlevels(i-1))then
               print *,'ellderivs: error: rlevels have to be ascending'
               stop
            endif
         enddo


         if(ifremovestrain.eq.1)then
            do i=1,nvtimes
               if(vhdivmax(i).le.0)then
                  print *,'ellderivs: error: vhdivmax should be != 0'
                  print *,'ellderivs: i:',i,'vhdivmax: ',vhdivmax(i)
                  stop
               endif
            enddo
         endif
         init=1
      endif

c     decide if we need to calculate the velocity gradient matrix
      if(nwork.gt.3)then 
                                ! for vel diff loop
         icalc_vgm=.true. 
         idindex=3 
      else 
                                ! only velocities
         icalc_vgm=.false.
         idindex=0 
      endif

      !pick cell for tracer

      ! adjust to physical limits
      call check_phys_lim_tracer(XTRACER(1))

      !first, radial component
      ixtracer=0
      i=0
      do while ((ixtracer.eq.0).and.(i.lt.nr))
         i=i+1
         if (XTRACER(1).le.rlevels(i)) ixtracer=i
      end do
      if(ixtracer .eq. 0)then  ! no depth levels found, tracer is above surface
         ixtracer = nr ! assign last layer, x_r should be corrected by the RK routines
      endif


      !pick grid points
      do i=0,3
         igrid(i,1)=ixtracer-1+i
      end do

      !make sure all grid points exist
      do while (igrid(0,1).le.0)
         do i=0,3
            igrid(i,1)=igrid(i,1)+1
         end do
      end do
      do while (igrid(3,1).gt.nr)
         do i=0,3
            igrid(i,1)=igrid(i,1)-1
         end do
      end do

      !find values of r for each grid point
      do i=0,3
         grid(i,1)=rlevels(igrid(i,1))
      end do


      !calculate 1-d weights for each grid point
      ngrid=3 ! ie. fourth order since this is n-1
      iorder=1
      call weights(XTRACER(1),grid(0,1),ngrid,iorder,weight(0,0,1))

      !now, do theta component
      ixtracer=nint(XTRACER(2)/dtheta)
      
      !pick grid points
      do i=0,3
         igrid(i,2)=ixtracer-1+i
      end do

      !adust grid points to avoid wrap-around
      do while (igrid(0,2).le.0)
         do i=0,3
            igrid(i,2)=igrid(i,2)+1
         end do
      end do
      do while (igrid(3,2).gt.ntheta)
         do i=0,3
            igrid(i,2)=igrid(i,2)-1
          end do
      end do
      !
      ! find values of theta: since given on dtheta/2 .... pi-dtheta/2
      ! theta_i = (i-0.5)*dtheta
      !
      do i=0,3
         grid(i,2)=(igrid(i,2)-.5d0)*dtheta
      end do

      !calculate 1-d weights for each grid point
      call weights(XTRACER(2),grid(0,2),ngrid,iorder,weight(0,0,2))
      
      ixtracer=int(XTRACER(3)/dphi)+1

      !pick grid points
      do i=0,3
         igrid(i,3)=ixtracer-1+i
      end do
      !
      ! find values of phi. phi_i = i * dphi
      !
      do i=0,3
         grid(i,3)=(igrid(i,3)-1)*dphi
      end do
      !calculate 1-d weights for each grid point
      call weights(XTRACER(3),grid(0,3),ngrid,iorder,weight(0,0,3))
      
      !now make sure grid points exist
      do i=0,3
         if (igrid(i,3).lt.1) igrid(i,3)=igrid(i,3)+nphi
         if (igrid(i,3).gt.nphi) igrid(i,3)=igrid(i,3)-nphi
      end do

#ifdef DEBUG 
c     check if all indices are ok
      do i=0,3
         if((igrid(i,1).lt.1).or.(igrid(i,1).gt.nr))then
            print *,'ellsphere: row ',i,'r index ',igrid(i,1), ' error'
            stop
         endif
         if((igrid(i,2).lt.1).or.(igrid(i,2).gt.ntheta))then
            print *,'ellsphere: row ',i,'th index ',igrid(i,2), ' error'
            stop
         endif
         if((igrid(i,3).lt.1).or.(igrid(i,3).gt.nphi))then
            print *,'ellsphere: row ',i,'ph index ',igrid(i,3), ' error'
            stop
         endif
      enddo
#endif
c
c
c     check if we are in polar region, ie. closer than one grid point from
c     pole
c
      if ((XTRACER(2).gt.dtheta).and.((pi-XTRACER(2)).gt.dtheta)) then
c
c     NOT in polar regions
c
c     first calculate velocities only (idindex=0) or vel and  derivatives 
c     of velocity (e.g. v_(r,r)) if needed (idindex=3)
c
         do m=0, idindex        !m=0 -> no derivative
            dvr(m)=0.d0         !m=1 -> derivative wrt r
            dvtheta(m)=0.d0     !m=2 -> derivative wrt theta
            dvphi(m)=0.d0       !m=3 -> derivative wrt phi
            do msub=0,3
               ider(msub)=0
            end do
            ider(m)=1 
            do i=0,3
               do j=0,3
                  do k=0,3
                     rnet=weight(i,ider(1),1)*
     &                    weight(j,ider(2),2)*
     &                    weight(k,ider(3),3)
                     index=(igrid(i,1)-1)*ntnp + 
     &                    (igrid(j,2)-1)*nphi +
     &                    igrid(k,3)
                     call get_velocities(vrloc,vthetaloc,vphiloc,
     &                    vr,vtheta,vphi,vtimes,nvtimes,index,
     &                    nrntnp,t)
                     dvr(m)=    dvr(m)     + rnet * vrloc
                     dvtheta(m)=dvtheta(m) + rnet * vthetaloc
                     dvphi(m)=  dvphi(m)   + rnet * vphiloc
                  end do
               end do
            end do         
         end do
         
      else    
c
c     polar regions
c     
c         print *,'polar! theta(for example):',XTRACER(2)
c
         do m=0,idindex
            dvr(m)=0.d0      
            dvtheta(m)=0.d0  
            dvphi(m)=0.d0    
            do msub=0,3
               ider(msub)=0
            end do
            ider(m)=1 
            do i=0,3            ! radial node
               do k=0,3         !phi node
                  if (m.ne.2) then ! just use values for theta=dtheta
                     j=1
                     !check if n pole or s pole
                     if (pi-XTRACER(2).lt.dtheta) j=ntheta
                     index=(igrid(i,1)-1)*ntnp+(j-1)*nphi + 
     &                    igrid(k,3)
                     rnet=0.5*weight(i,ider(1),1)*weight(i,ider(3),3)
                     call get_velocities(vrloc,vthetaloc,vphiloc,
     &                    vr,vtheta,vphi,vtimes,nvtimes,
     &                    index,nrntnp,t)
                     dvr(m)=dvr(m)+rnet*vrloc
                     dvtheta(m)=dvtheta(m)+rnet*vthetaloc
                     dvphi(m)=dvphi(m)+rnet*vphiloc
                  else
                  !for theta derivative, project derivative in from 
                  ! two nodes closest to pole
                     do j=0,1
                        jindex=j
                        ishift=1 !see rnet calculation.  this gets the
                                 !sign of the weight right for the deriv
                        if (pi-XTRACER(2).lt.dtheta) then !s pole
                           jindex=ntheta-j
                           shift=0
                        end if
                        index=(igrid(i,1)-1)*ntnp+(jindex-1)*nphi + 
     &                       igrid(k,3)
                        rnet=((-1)**(j+ishift)/dtheta)*
     &                       weight(i,ider(1),1)*
     &                       weight(i,ider(3),3)
                        call get_velocities(vrloc,vthetaloc,vphiloc,
     &                       vr,vtheta,vphi,
     &                       vtimes,nvtimes,index,nrntnp,t)
                        dvr(m)=dvr(m)+rnet*vrloc
                        dvtheta(m)=dvtheta(m)+rnet*vthetaloc
                        dvphi(m)=dvphi(m)+rnet*vphiloc
                     end do
                  end if
               end do
            end do
         end do
      end if
#ifdef DEBUG
      if(XTRACER(1).le.0.0)then
         print *,'ellderivs: out of range error'
         print *,'ellderivs: r: ',XTRACER(1)
         stop
      endif
#endif       
c  
c     factors
c
      finvr=1.d0/XTRACER(1)    ! 1/r
c
c     the first three derivatives of the location vector are the spherical
c     system velocities
c      
c     the workder(3) part had a 1/sin(theta) factor (wrongly) 
c     in front of it
c     
      workder(1)=dvr(0)
      workder(2)=finvr*dvtheta(0)
      workder(3)=finvr*dvphi(0)
c
c     calculate the velocity gradient matrix
c
      if(icalc_vgm)then
         sintheta=sin(XTRACER(2)) ! sin(theta)
         costheta=cos(XTRACER(2))
         cottheta=1.d0/tan(XTRACER(2))
c
c     build velocity gradient matrix. transpose of grad v
c
c     (see, e.g., Dahlen & Tromp, p. 836)
c     
c     vgm(i,j) = d_j v_i
         
         vgm(1,1)=dvr(1)                       ! d_r v_i
         vgm(2,1)=dvtheta(1)                   
         vgm(3,1)=dvphi(1)                     
         vgm(1,2)=finvr*(dvr(2)-dvtheta(0))    ! d_theta v_i
         vgm(2,2)=finvr*(dvtheta(2)+dvr(0))    
         vgm(3,2)=finvr*dvphi(2)               
         vgm(1,3)=finvr*(dvr(3)/sintheta-dvphi(0)) ! d_phi v_i
         vgm(2,3)=finvr*(dvtheta(3)/sintheta-cottheta*dvphi(0)) 
         vgm(3,3)=finvr*(dvphi(3)/sintheta+dvr(0)+cottheta*dvtheta(0)) 

c
c     calculate divergence and adjust vgm matrix so that the divergence
c     is zero (not ideal, but functional)
c
         div=(vgm(1,1)+vgm(2,2)+vgm(3,3))/3.0d0
         do i=1,3
            vgm(i,i)=vgm(i,i)-div
         end do
c         print *,vgm(1,1)+vgm(2,2)+vgm(3,3)
c     
c     remove the strain part of the velocity gradient tensor
c     
         if(ifremovestrain.eq.1)then 
c     horizontal component of the divergence
            vhdiv = vgm(2,2) + vgm(3,3) 
c
c     get the max horizontal divergence at time t
            call get_vhdivmax(t,vtimes,nvtimes,vhdivmax,vhdm_loc)
c
c     determine an empirical depth z of the T = 0.5 T_m isoline to check if 
c     the tracer is shallower than that depth. if so, freeze in the strain by
c     removing the symmetric part of the velocity gradient matrix so that 
c     there's only rotations affecting the strain matrix F
            if(XTRACER(1).gt.0.98901267310195d0)then 
c     tracer is shallower than ~70km
               vhdivnd = vhdiv/vhdm_loc ! normalized horizontal divergence
               if(vhdivnd.gt.0.91047619047619d0)then 
c     strongly positive divergence -> ridge
                  rlim = 1.0d0  ! ie., all strains matter, z_T0.5 = 0
               else if(vhdivnd.lt.-0.00527472527472533d0) then 
c     strongly negative, away from ridge -> z_T0.5 = 70km
                  rlim = 0.98901267310195d0
               else             
c     calculate depth of 0.5 T_m from an upper bound for z_T05 from 
c     `experimental' data. the approx is z(x) = -21+25/(x+0.28) [km], x is
c     relative horizontal divergence. this is converted into r=(1-z/6371)
                  rlim = 1.0032961858421d0 - 0.00392403076440119d0 / 
     &                 (vhdivnd + 0.28d0)
               endif
               if(XTRACER(1).gt.rlim)then ! if shallower than critical temperature
                                          ! for this `plate age', freeze in strain
                  call remove_symm_part(vgm)   
               endif
            endif 
c     end of remove-symmetrix-strain loop
         endif
c     
c
c     the time derivative of F is L dot F, repack into workder array
c
         index=3
         do i=1,3
            do j=1,3
               index=index+1
               workder(index)=vgm(i,1)*DEF1(j)+vgm(i,2)*DEF2(j)+
     &              vgm(i,3)*DEF3(j)
            end do
         end do
      endif

      return
      end
c
c
c     remove the symmetric part, D, 
c     of the velocity gradient matrix, G, by 
c     calculating the vorticity, W
c     
c     G = D + W ==> W = G - D 
c                     = G - 0.5(G + G^T) 
c                     = 0.5 ( G - G^T ) = W
c
c     returns the vorticity W
c     
      subroutine remove_symm_part(g)
#include "fortran_precision.h"
      dimension w(3,3),g(3,3)
      do i=1,3
         do j=1,3
            w(i,j) = 0.5d0 * (g(i,j) - g(j,i))
         enddo
      enddo
      do i=1,3
         do j=1,3
            g(i,j) = w(i,j)
         enddo
      enddo
      end
!
! check the physical bounds of a spherical system point
! in r,theta,phi space
!
      subroutine check_phys_lim_tracer(xloc)
#include "fortran_precision.h"
      parameter(r_cmb = 0.5493)
      COMP_PRECISION xloc(3)
      !
      ! adjust position values

      if(xloc(1).gt.1.0)then
         print *,'check_phys_lim_tracer: r surface adjust: ',xloc(1)
         xloc(1)=1.0d0
      endif
      if(xloc(1).lt.r_cmb)then
         print *,'check_phys_lim_tracer: r CMB adjust: ',xloc(1)
         xloc(1) = r_cmb
      endif
      !
      ! check the theta range
      !
      ! 0 <= theta <= pi
      !
      do while (xloc(2).gt.twopi)
         xloc(2)=xloc(2)-twopi
      end do
      do while (xloc(2).lt.0.d0)
         xloc(2)=xloc(2)+twopi
      end do
      if (xloc(2).gt.pi) then
         xloc(2)=twopi - xloc(2)
         xloc(3)=xloc(3) + pi
      end if
      !
      ! check the phi range
      !
      ! 0<= phi <= 2*pi
      do while (xloc(3).gt.twopi)
         xloc(3)=xloc(3)-twopi
      end do
      do while (xloc(3).lt.0.d0)
         xloc(3)=xloc(3)+twopi
      end do


      end 
!
!     find a physics a priori knowledge based error condition
!
      subroutine check_physics_based_error(y,dydx,n,eps,errmax,h)
#include "fortran_precision.h"
      COMP_PRECISION y(n),dydx(n),eps,errmax,h,tmp_err
      integer n
c
c     check for radial direction coordinate
c
      if(y(1)-1.0.gt.eps)then
c
c     check for tracer above surface at next step
c
         tmp_err = (y(1) - 1.0d0)/eps
         errmax = MAX(tmp_err,errmax)
#ifdef DEBUG
         print *,'check_phy_based error: dydx/h/r:',dydx(1),h,y(1)
#endif
      endif
      end
