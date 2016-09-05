c     
c     calculate time derivatives of work array (output: workder) at time t
c     given the velocity arrays vr,vtheta,vphi (dimensions: nr*ntheta*nphi*nvtimes)
c     specified at times vtimes(nvtimes*3) (left, mid, right end of interval) 
c     and rlevels(nr)
c     
c     
c     the time derivative of work(1...3), the location x, is the change
c     in the coordinates at x in a spherical coordinate system (not the velocities)
c
c     d_t r      =                  u_r
c     d_t theta  = 1/r              u_t
c     d_t phi    = 1/r 1/sin(theta) u_p
c
c     the velocities are output as vout(3). (see below for ref. frame)
c
c     the time derivative of work(4..12), the strain F, is the velocity gradient
c     matrix, G,  dotted with the strain, ie. G . F
c     
c     WARNING: when nwork <= 3, the routine will not determine
c     the time derivate fo the strain but only the velocities at location x, ie.
c     the time derivatives of work(1...3)
c     
c     ifremovestrain: 
c     
c     flag to remove the symmetric strain and only rotate F, for this we
c     also need vhdivmax, the max horizontal divergence used to
c     determine the plate `age'
c     
c     if nwork.gt.3: 

c     fse_derivs also returns vgm_c, the velocity gradient matrix, in C
c     style sorting. if ifremovestrain is set, the VGM matrix will
c     reflect this
c     
C     IF IDER_CARTESIAN IS UNITY, VGM, THE CORRESPONDING TIME
C     DERIVATIVES, AND THE OUTPUT VELOCITIES WILL BE IN A CARTESIAN
C     SYSTEM. ELSE, THEY WILL BE IN SPHERICAL.
c
c     ifraction_strain(input): is set to unity, will determine the degree to 
c     which the actual strain-rate matrix (as opposed to the rotational matrix)
c     will be used according to alpha
c
c
c ONLY FOR THE DEBUGGING VERSION:
c
c     if ivgm_constant is set, will use a constant velocity gradient
c     matrix for debugging purposes. if this is the case, then the [3]
c     vector cvgm_stw will hold the S, T, and W values from Mckenzie and
c     Jackson (1982)
c
c     dextral simple shear:  S=T=0, W>0
c     pure shear:            T=W=0, S>0
c     isostatic contraction: W=S=0, T>0
c     uniaxial contraction:  S=T>0, W=0
c
c     the coordinate system is rotated such that x is 0,0, y is 90, 0, 
c     and z is 0,90 (the north pole), in a lon,lat system. ie. x = -y' and 
c     y=x', where x',y' are the normal 2-D coordinates as used in McKenzie &
c     Jackson
c  END DEBUGGING VERSION
c
c
c     ifrozen (output) is 1 if the tracer strain will actually be removed (ifremovestrain),
c     else it's zero
c     
c     
      subroutine fse_derivs(time,work,workder,nwork,
     &     vr,vtheta,vphi,nr,ntheta,nphi,dtheta,dphi,rlevels,
     &     vtimes,nvtimes,ifremovestrain,vhdivmax,vout,vgm_c,
     &     ider_cartesian,ivgm_constant,cvgm_stw,ifrozen,
     &     ialpha_strain,alpha) 
c_______________________________________________________________________
c     
c     Note:  passed muti-dimensional arrays are in c-order
c     
c_______________________________________________________________________
      
#include "fortran_precision.h"

      VPRECF vr,vtheta,vphi

      dimension vr(*), vtheta(*), vphi(*), rlevels(nr),vout(3),voutr(3),
     &     vtimes(*),vhdivmax(*),vgm_f(3,3),xloc(3),vgm_c(9),vgm_cc(9),
     &      vgm_f2(3,3)
      dimension work(nwork),workder(nwork),dtlim(2)
      dimension dvr(0:3),dvtheta(0:3),dvphi(0:3),cvgm_stw(3)
      integer ifrozen
      logical icalc_vgm
      save dtlim
      save init
      data init /0/

      ifrozen = 0 
      pih = acos(0.0d0)         ! pi/2.0

      if(init.eq.0)then         
c
c     init some variables
c     
         dtlim(1) = dtheta / 20.0
         dtlim(2) = pi - dtlim(1)
c     
         if(ifremovestrain.eq.1)then
            print *,'fse_derivs: WARNING: removing shallow strain'
            do i=1,nvtimes
               if(vhdivmax(i).le.0)then
                  print *,'fse_derivs: vhdivmax should be not be 0'
                  print *,'fse_derivs: i:',i,'vhdivmax: ',vhdivmax(i)
                  stop
               endif
            enddo
         endif
         if(ialpha_strain.eq.1)then
            print *,'fse_derivs: WARNING: removing fraction of strain'
            print *,'fse_derivs: according to alpha'
         endif

         init=1
      endif
c     
c     decide if we need to calculate the velocity gradient matrix
      if(nwork.gt.3)then 
                                ! for vel diff loop
         icalc_vgm=.true. 
      else 
                                ! only velocities
         icalc_vgm=.false.
      endif
c     
c     deflect the polar paths
      do while (work(_T_).lt.dtlim(1))
         work(_T_) = dtlim(1)
         work(_P_) = work(_P_) + pi
      end do
      if (work(_T_).gt.dtlim(2)) then
         work(_T_)=dtlim(2)
         work(_P_)=work(_P_) + pi
      end if
c     
c     pick cells for tracer
c     
c     adjust tracer location 
c     to physical limits (0<=phi<=2pi, 0<=theta<=pi)
c     but leave original state untouched
c     
      call check_phys_lim_tracer(work(_R_),xloc)
c
c     call the interpolation routine and return the velocities dvr(0), 
c     dvtheta(0), dvphi(0), and their r,theta,phi derivatives dv_(1,2,3)
c
      call find_vel_and_der(xloc,time,nr,ntheta,nphi,
     &     vr,vtheta,vphi,rlevels,dtheta,dphi,vtimes,nvtimes,icalc_vgm,
     &     dvr,dvtheta,dvphi)
c     
c     factors
c     
      finvr=1.d0/xloc(_R_)      ! 1/r
      sintheta=sin(xloc(_T_))   ! sin(theta)
c     
c     the first three derivatives of the location vector are NOT the
c     spherical system velocities but related!
c     
      workder(1)=dvr(0)                            ! d_t r =         u_r
      workder(2)=finvr * dvtheta(0)                ! d_t theta = 1/r u_theta
      workder(3)=finvr * dvphi(0) / sintheta       ! d_t phi = 1/(r sin(theta)) u_phi
c
c     those are the spherical velocities
c
      vout(_R_) = dvr(0)
      vout(_T_) = dvtheta(0)
      vout(_P_) = dvphi(0)
      if(ider_cartesian.eq.1)then
c
c     convert velocities to Cartesian
c         
         call polar_vec2cart_vec_at_xp(vout(_R_),voutr,xloc(_R_))
         do i=1,3
            vout(i) = voutr(i)
         enddo
      endif
c     
c     calculate the velocity gradient matrix
c     
      if(icalc_vgm)then
         costheta=cos(xloc(_T_))
         cottheta=1.d0/tan(xloc(_T_))
c     
c     build velocity gradient matrix. transpose of grad v
c     
c     (see, e.g., Dahlen & Tromp, p. 836)
c     
c     vgm_f(i,j) = d_j v_i
c
c     | d_r(vr) d_t(vr) d_p(vr) |
c     | d_r(vt) d_t(vt) d_p(vt) |
c     | d_r(vp) d_t(vp) d_p(vp) |
c     
c     
         vgm_f(1,1)=dvr(1)      ! d_r v_i, first column
         vgm_f(2,1)=dvtheta(1)                   
         vgm_f(3,1)=dvphi(1)                     
         vgm_f(1,2)=finvr*(dvr(2)-dvtheta(0)) ! d_theta v_i, second column
         vgm_f(2,2)=finvr*(dvtheta(2)+dvr(0))    
         vgm_f(3,2)=finvr*dvphi(2)               
         vgm_f(1,3)=finvr*(dvr(3)/sintheta-dvphi(0)) ! d_phi v_i, third column
         vgm_f(2,3)=finvr*(dvtheta(3)/sintheta-cottheta*dvphi(0)) 
         vgm_f(3,3)=finvr*(dvphi(3)/sintheta+dvr(0)+
     &        cottheta*dvtheta(0)) 
c     
c     calculate divergence and adjust vgm_f matrix so that the divergence
c     is zero (not ideal, but functional)
c     
         div=(vgm_f(1,1)+vgm_f(2,2)+vgm_f(3,3))/3.0d0
#ifdef FSTRACK_DEBUG
      if(abs(div).gt.5e-2)then
         print *,'fse_deriv: t:',time,' div (before corr):',div
      endif
#endif
      do i=1,3
         vgm_f(i,i)=vgm_f(i,i)-div
      end do
c
c     ifremovestrain is set:
c
c     for shallow regions of the plate ("frozen") remove the strain part of the velocity gradient tensor 
c     and only leave the rotational components
c
      if(ifremovestrain.eq.1)then 
         call determine_rs_status(time,vtimes,nvtimes,vhdivmax,vgm_f,
     &        xloc,ifrozen)
         if(ifrozen.eq.1)then
            call remove_symm_part(vgm_f,0.0d0)
         endif
      endif
c
c
c     remove a fraction of the strain according to alpha
c
c 
      if(ialpha_strain.eq.1)then
c         write(*,678)xloc(_P_)*57.295779513082320,
c     &        90.-xloc(_T_)*57.295779513082320,
c     &        (1-xloc(_R_))*6371,alpha
 678     format(4(f12.3,1x))
         call remove_symm_part(vgm_f,alpha)
      endif


#ifdef FSTRACK_DEBUG
      if(ivgm_constant.eq.1)then
c
c
c set all velocities to zero, no transport
c
c
         workder(1) = 0.0d0
         workder(2) = 0.0d0
         workder(3) = 0.0d0
         vout(_R_) = 0.0d0
         vout(_T_) = 0.0d0
         vout(_P_) = 0.0d0

c
c     assign a constant velocity gradient matrix instead for debugging
c     purposes in Cartesian system. this assumes that x is east and 
c     y is north
c     
         vgm_f2(1,1)=cvgm_stw(1)-cvgm_stw(2) ! S - T     
         vgm_f2(2,1)=0.0d0
         vgm_f2(3,1)=0.0d0
         vgm_f2(1,2)=cvgm_stw(3) ! W
         vgm_f2(2,2)= -(cvgm_stw(1)+cvgm_stw(2)) ! -( S + T)
         vgm_f2(3,2)=0.0d0              
         vgm_f2(1,3)=0.0d0
         vgm_f2(2,3)=0.0d0
         vgm_f2(3,3)=2.0d0*cvgm_stw(2) ! 2T
c
c     rotate the matrix into our coordinate frame (leave this in here)
c     90 degrees CCW -> alpha = 90 = pi/2, beta = 0
c
         xalpha = pih
         xbeta =  0.0d0
         xgamma = 0.0d0
         call rotate_cart_mat_ftrn(vgm_f2,vgm_f,xalpha,xbeta,xgamma)

c     assign to C version
         k=0
         do i=1,3
            do j=1,3
               vgm_c(k + j) = vgm_f(i,j)
            enddo
            k=k+3
         enddo
         if(ider_cartesian.eq.0)then
c
c     need to convert to spherical, this routine operates on a C style matrix
c
            call cart_to_polar_mat_at_r(vgm_c,vgm_cc,xloc(_R_))
            k=0
            do i=1,3
               do j=1,3
                  l = k + j
                  vgm_f(i,j) = vgm_cc(l) 
                  vgm_c(l) = vgm_cc(l)
               enddo
               k=k+3
            enddo
         endif
      else
#endif
c
c     use the real velocity gradient matrix which is so far in spherical 
c     coordinates
c     
c     assign velocity gradient matrix to C style output
c     
         k=0
         do i=1,3
            do j=1,3
               vgm_c(k + j) = vgm_f(i,j)
            enddo
            k=k+3
         enddo
c     
c     
c     should we operate on the strain in Cartesian space?
c     
c
         if(ider_cartesian.eq.1)then
c     
c     convert the C-style matrix to cartesian
c     
            call polar_to_cart_mat_at_r(vgm_c,vgm_cc,xloc(_R_))
c
c     assign to the fortran style, transposed velocity gradient matrix
c     to Cartesian
c     
            k=0
            do i=1,3
               do j=1,3
                  l = k + j
c     assign to the fortran matrix
                  vgm_f(i,j) = vgm_cc(l) 
c     re-assign to the C style matrix
                  vgm_c(l) = vgm_cc(l)
               enddo
               k=k+3
            enddo
         endif                  ! end Cartesian branch
#ifdef FSTRACK_DEBUG
      endif                     ! end of ivmg not constant branch
#endif
c     
c     
c     the time derivative of F is G dot F, 
c     where G is the deformation rate tensor, vgm
c     
c     d_t F = G . F
c
c     repack into workder array
c     
      index=3
      do i=1,3
         do j=1,3
            index=index+1
c     c(i,j) = sum_k a(i,k) * b(k,j)
            workder(index)=
     &           vgm_f(i,1) * DEF1(j)+
     &           vgm_f(i,2) * DEF2(j)+
     &           vgm_f(i,3) * DEF3(j)
         end do
      end do
      endif
      
      return
      end
c     
c     
c     removes a fraction of the symmetric part, D, 
c     of the velocity gradient matrix, G
c     calculating the vorticity, W
c     
c     
c     if delta = 0, returns the vorticity W
c     if delta = 1, returns the full velocity gradient matrix
c     
      subroutine remove_symm_part(g,delta)
#include "fortran_precision.h"
      dimension w(3,3),g(3,3),d(3,3)
      do i=1,3 ! decompose the velocity gradient matrix
         do j=1,3
            w(i,j) = 0.5d0 * (g(i,j) - g(j,i)) ! vorticity
            d(i,j) = 0.5d0 * (g(i,j) + g(j,i)) ! strain-rate
         enddo
      enddo
      do i=1,3
         do j=1,3
            g(i,j) = w(i,j) + delta * d(i,j)
         enddo
      enddo
      end



      
