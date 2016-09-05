c
c
c     routines that deal with velocity interpolation in time and space
c
c     $Id: velinterpol.F,v 1.9 2011/04/12 06:19:09 becker Exp $
c
c
c     find the velocities (and their first derivatives, if icalc_der is set)
c     with respect to r,theta,phi
c
c     dvr(0,1,2,3):     v_r, d_r v_r, d_theta v_r, d_phi v_r
c     dvtheta(0,1,2,3): v_theta, d_r v_theta, d_theta v_theta, d_phi v_theta
c     dvphi(0,1,2,3):   v_phi, d_r v_phi, d_theta v_phi, d_phi v_phi
c     
c     if icalc_der is TRUE, will compute derivatives, else only velocities
c     
c     input are velocity fields vr,vtheta,vphi given laterally in dtheta/dphi
c     spaced ntheta*nphi layers at nr radial levels at rlevels[nr]
c
c
C     default order of interpolation - 1
c     ie. if 3, this means will do fourth order on 0,1,2,DEF_ORDER arrays
c#define DEF_ORDER 3
#define DEF_ORDER 2
c
c
c
      subroutine find_vel_and_der(xloc,time,nr,ntheta,nphi,
     &     vr,vtheta,vphi,rlevels,dtheta,dphi,vtimes,nvtimes,
     &     icalc_der,dvr,dvtheta,dvphi)
#include "fortran_precision.h"
c
c     output:
      dimension dvr(0:3),dvtheta(0:3),dvphi(0:3)
c
c     input:

c those have possibly different 
      VPRECF vr,vtheta,vphi
      
      dimension vr(*),vtheta(*),vphi(*),vtimes(*),xloc(3),
     &     rlevels(nr)
c     derivatives yes/no?
      logical icalc_der
c
c     order of derivatives needed in general for the weighting routine
      parameter (iorder=1)
c
c     internal arrays for interpolation
      dimension istencil(3),isshift(3),ider(0:3)
      dimension igrid(0:DEF_ORDER,3)
      dimension grid(0:DEF_ORDER,3),weight_r(0:DEF_ORDER,0:iorder),
     &     weight_t(0:DEF_ORDER,0:iorder),weight_p(0:DEF_ORDER,0:iorder)
c
      save ntnp,nrntnp,init,ider,istencil,isshift
      data init /0/
      if(init.eq.0)then
c     
c     do some checks if called for the first time
c     
         if(nr.lt.DEF_ORDER+1)then
            print *,'find_vel_and_der: need at least four depth levels'
            print *,'find_vel_and_der: using polynomial interpolation'
            stop
         endif
         if(nphi.lt.DEF_ORDER+1)then
            print *,'find_vel_and_der: need at least four lon levels'
            print *,'find_vel_and_der: using polynomial interpolation'
            stop
         endif
         if(ntheta.lt.DEF_ORDER+1)then
            print *,'find_vel_and_der: need at least four lat levels'
            print *,'find_vel_and_der: using polynomial interpolation'
            stop
         endif
         do i=2,nr
            if(rlevels(i).le.rlevels(i-1))then
               print *,'find_vel_and_der: error:'
               print *,'rlevels have to be ascending'
               print *,i,rlevels(i),rlevels(i-1)
               stop
            endif
         enddo
c
c     set some of the parameters
c     
         ntnp = ntheta * nphi
         nrntnp = ntnp * nr
c
c     these are offsets for the interpolation routine
c
         do i=1,3
            istencil(i) = DEF_ORDER
            isshift(i) = int(istencil(i)/2.0)
         enddo
c     
c     this is for derivatives, initialize once as zeroes
         do i=0,3
            ider(i)=0
         end do
         init = 1
      endif
c
c     
c     for the summing up routine of weights
c
      if(icalc_der)then         ! compute derivatives
         idindex=3
      else
         idindex=0
      endif
c     
c     first, find radial component
c     
c     default order-1 interpolation for radial direction
      ixtracer=0
      i=0
      do while ((ixtracer.eq.0).and.(i.lt.nr))
         i=i+1
         if (xloc(_R_).le.rlevels(i)) ixtracer=i
      end do
      if(ixtracer .eq. 0)then   ! no depth levels found, tracer is above surface
         ixtracer = nr          ! assign last layer, x_r should be corrected by the RK routines
      endif
c     
c     pick indices of grid points for the radial stencil
c     
      do i=0,istencil(_R_)
         igrid(i,_R_) = ixtracer - isshift(_R_) + i
      end do
c     
c     make sure all grid points exist 
c     
      ishift = igrid(0,_R_) - 1
      if(ishift.lt.0)then
         do i=0,istencil(_R_)
            igrid(i,_R_) = igrid(i,_R_) - ishift
         end do
      endif
c     same for upper limit
      ishift = igrid(istencil(_R_),_R_) - nr
      if(ishift.gt.0)then
         do i=0,istencil(_R_)
            igrid(i,_R_)=igrid(i,_R_) - ishift
         enddo
      endif
c     find values of r for each grid point
      do i=0,istencil(_R_)
         grid(i,_R_)=rlevels(igrid(i,_R_))
      end do
c     
c     theta component
c
      ixtracer = nint(xloc(_T_)/dtheta)
c     pick grid points
      do i=0,istencil(_T_)
         igrid(i,_T_)=ixtracer-isshift(_T_)+i
      end do
c     
c     adust grid points to avoid wrap-around
c     
c     
      ishift = igrid(0,_T_) - 1 
      if(ishift.lt.0)then
         do i=0,istencil(_T_)
            igrid(i,_T_) = igrid(i,_T_) - ishift
         end do
      endif
c     same for upper limit
      ishift = igrid(istencil(_T_),_T_) - ntheta 
      if(ishift.gt.0)then
         do i=0,istencil(_T_)
            igrid(i,_T_)=igrid(i,_T_) - ishift
         enddo
      endif
                                !
                                ! find values of theta: since given on dtheta/2 .... pi-dtheta/2
                                ! theta_i = (i-0.5)*dtheta
                                !
      do i=0,istencil(_T_)
         grid(i,_T_)=(igrid(i,_T_)-.5d0)*dtheta
      end do
c     
c     now for phi
c     
      ixtracer=int(xloc(_P_)/dphi)+1
                                !pick grid points
      do i=0,istencil(_P_)
         igrid(i,_P_)=ixtracer-isshift(_P_)+i
c     
c     wrap around 
c     
         if(igrid(i,_P_).gt.nphi)then
            igrid(i,_P_) = igrid(i,_P_) - nphi
         endif
         if(igrid(i,_P_).lt.1)then 
            igrid(i,_P_) = igrid(i,_P_) + nphi
         endif
      end do
                                !
                                ! find values of phi. phi_i = i * dphi
                                !
      do i=0,istencil(_P_)
         grid(i,_P_)=(igrid(i,_P_)-1)*dphi
      end do
                                !calculate 1-d weights for each grid point
#ifdef FSTRACK_DEBUG
c     
c     check if all indices are ok
c     
      do i=0,istencil(_R_)
         if((igrid(i,_R_).lt._R_).or.(igrid(i,_R_).gt.nr))then
            print *,'find_vel_and_der: row ',i,' r index ',
     &           igrid(i,_R_), 'error'
            stop
         endif
      enddo
      do i=0,istencil(_T_)
         if((igrid(i,_T_).lt.1).or.(igrid(i,_T_).gt.ntheta))then
            print *,'find_vel_and_der: row ',i,' theta index ',
     &           igrid(i,_T_), 'error'
            stop
         endif
      enddo
      do i=0,istencil(_P_)
         if((igrid(i,_P_).lt.1).or.(igrid(i,_P_).gt.nphi))then
            print *,'find_vel_and_der: row ',i,' phi index ',
     &           igrid(i,_P_), 'error'
            stop
         endif
      enddo
#endif

c     
c     POLYNOMIAL
c     
c     compute the weights for each stencil
      call weights(xloc(_R_),grid(0,_R_),istencil(_R_),
     &     iorder,weight_r(0,0))
      call weights(xloc(_T_),grid(0,_T_),istencil(_T_),
     &     iorder,weight_t(0,0))
      call weights(xloc(_P_),grid(0,_P_),istencil(_P_),
     &     iorder,weight_p(0,0))
      
c     
c     first calculate velocities only (idindex=0) or vel and  derivatives 
c     of velocity (e.g. v_(r,r)) if needed (idindex=3)
c      
      do m=0, idindex           !m=0 -> no derivative
         dvr(m)=0.d0            !m=_R_ -> derivative wrt r
         dvtheta(m)=0.d0        !m=_T_ -> derivative wrt theta
         dvphi(m)=0.d0          !m=_P_ -> derivative wrt phi
c     
c     this is the derivative yes/no array
c     set once, and switch off again below!
c     
         ider(m) = 1 
         do i=0,istencil(_R_)   ! radial 
            do j=0,istencil(_T_) ! theta 
               do k=0,istencil(_P_) ! phi
                  rnet=weight_r(i,ider(_R_))*
     &                 weight_t(j,ider(_T_))*
     &                 weight_p(k,ider(_P_))
                  index=(igrid(i,_R_)-1)*ntnp + 
     &                 (igrid(j,_T_)-1)*nphi +
     &                 igrid(k,_P_)
                  call get_velocities(vrloc,vthetaloc,vphiloc,
     &                 vr,vtheta,vphi,vtimes,nvtimes,index,
     &                 nrntnp,time)
                  dvr(m)=    dvr(m)     + rnet * vrloc
                  dvtheta(m)=dvtheta(m) + rnet * vthetaloc
                  dvphi(m)=  dvphi(m)   + rnet * vphiloc
               end do
            end do
         end do
         ider(m)=0              ! reset derivative switxh  to zero
      end do                    ! end m loop
      

      end



c     get_velocities
c     inter_vel
c     weights: for spatial interpolation
c     
c     $Id: velinterpol.F,v 1.9 2011/04/12 06:19:09 becker Exp $
c
c
c     obtain the time-interpolated velocities while index
c     specifies the 3-D position in the arrays that are
c     vr(nrntnp*nvtimes) long. the vtimes array is nvtimes*3 and has 
c     t_left t_mid t_right for each interval in a row
c
c
      subroutine get_velocities(vrloc,vthetaloc,vphiloc,
     &     vr,vtheta,vphi,vtimes,nvtimes,index,nrntnp,time)
#include "fortran_precision.h"

      VPRECF vr,vtheta,vphi

      dimension vr(*),vtheta(*),vphi(*),vtimes(*)
#ifdef FSTRACK_DEBUG
      if((index.lt.1).or.(index.gt.nrntnp*nvtimes))then
         print *,'get_velocities: error: index ',index
         print *,'get_velocities: bounds: nrntntp:',nrntnp, 
     &        'nvtimes:',nvtimes
         stop
      endif
#endif
      if(nvtimes.eq.1)then! only one time-step, steady-state calculation
         vrloc=      dble(vr(index))
         vthetaloc = dble(vtheta(index))
         vphiloc=    dble(vphi(index))
      else ! interpolate in time
         call inter_vel(time,vtimes,nvtimes,i1,i2,vf1,vf2)
         index1=(i1-1) * nrntnp + index
         index2=(i2-1) * nrntnp + index
         if(vf1.ne.0.0)then
            vrloc=      dble(vr(index1))     * vf1 
            vthetaloc = dble(vtheta(index1)) * vf1 
            vphiloc=    dble(vphi(index1))   * vf1 
         else
            vrloc=0.0d0
            vthetaloc=0.0d0
            vphiloc=0.0d0
         endif
         if(vf2.ne.0.0)then
            vrloc=      vrloc     + dble(vr(index2))     * vf2 
            vthetaloc = vthetaloc + dble(vtheta(index2)) * vf2 
            vphiloc=    vphiloc   + dble(vphi(index2))   * vf2 
         endif
      endif
      end
c
c     interpolate velocities given at vtimes(nvtimes * 3) to
c     time 'time'. returns factors and time levels, see above
c
c     WARNING: if the routine gets called with the same time twice, will
c     not recalculate the weights
c
      subroutine inter_vel(time,vtimes,nvtimes,i1,i2,vf1,vf2)
#include "fortran_precision.h"
c     dxlimit is the width of the transition between velocity stages
      parameter(dxlimit=0.01d0,
     &     xllimit= -dxlimit/2.0d0,
     &     xrlimit=  dxlimit/2.0d0)
      dimension vtimes(*)
c
c     these are the local, saved, quantities which are not determined anew
c     if the time doesn't change
c     
      save vf1_loc,vf2_loc,time_old
      save ileft_loc,iright_loc
      save icalled
      data icalled/0/

c
c
      tloc = time
c
c     special case: only one interval
c
      if(nvtimes.eq.1)then
         i1 = 1
         i2 = 1
         vf1 = 1.0d0
         vf2 = 0.0d0
      else ! more than one stage
         if((icalled.eq.0).or.(tloc.ne.time_old))then 
c     
c     heave to determine interpolation factors, time has changed
c     or not called at least once
c     
c     obtain the time intervals to the left and right of time
c     
            if(tloc.lt.vtimes(1))then
               if(abs(tloc-vtimes(1)) .lt. 0.1)then
                  tloc = vtimes(1)
               else
                  print *,'inter_vel: error, time: ',tloc, ' too small'
                  print *,'inter_vel: while vtimes are given from',
     &                 vtimes(1),' to ',vtimes(nvtimes*3)
                  stop  
               endif
            endif
            if(tloc.gt.vtimes(nvtimes*3))then
               if(abs(tloc-vtimes(nvtimes*3)) .lt. 0.1)then
                  tloc = vtimes(nvtimes*3)
               else
                  print *,'inter_vel: error, time: ',tloc, ' too large'
                  print *,'inter_vel: while vtimes are given from',
     &                 vtimes(1),' to ',vtimes(nvtimes*3)
                  stop
               endif
            endif
            iright_loc=1  ! right interval index
            i22=2         ! right interval midpoint
c     find the right interval such that its midpoint is larger or equal than time
            do while((tloc.gt.vtimes(i22)).and.(iright_loc.lt.nvtimes))
               iright_loc = iright_loc + 1
               i22 = i22 + 3
            enddo
            if(iright_loc.eq.1)then
               iright_loc = 2
               i22 = 5
            endif
            ileft_loc = iright_loc - 1 ! left interval index
c     distance from right boundary of left interval 
c     (=left boundary of right interval) normalized by the mean interval width 
            xll = 2.0d0 * (tloc-vtimes(i22-1))/
     &           (vtimes(i22) - vtimes(ileft_loc*3-2))
c     this will have xll go from -0.5 to 0.5 around the transition between stages
c     which is at xl1=0
c
c     vf1_loc and vf2_loc are the weights for velocities within the left and right
c     intervals respectively
            if(xll.lt.xllimit)then ! xllimit should be 1-dx, dx~0.1
               vf1_loc = 1.0d0
               vf2_loc = 0.0d0
            else
               if(xll.gt.xrlimit)then ! xrlimit should be 1+dx, dx~0.1
                  vf1_loc = 0.0d0
                  vf2_loc = 1.0d0
               else             ! in between 
                  xll =     (xll-xllimit)/dxlimit ! normalize by transition width
                  vf2_loc = ((1.0d0 - cos(xll * pi))/2.d0) ! this goes from 0 to 1
c     weight for left velocities
                  vf1_loc = 1.0d0 - vf2_loc
               endif
            endif
            icalled = 1
            time_old = tloc
         endif 
c
c     assign the local variables to the output variables
c     
         i1 = ileft_loc
         i2 = iright_loc
         vf1 = vf1_loc 
         vf2 = vf2_loc
      endif ! end ntimes>1 part
      end

c
c     determine the max horizontal divergence `vhdm' at time `time' from
c     the `vhdivmax' array. this array is based on vgm_f(theta,theta)
c     and vgm_f(phi,phi)
c
      subroutine get_vhdivmax(time,vtimes,nvtimes,vhdivmax,vhdm)
#include "fortran_precision.h"
      dimension vtimes(*),vhdivmax(*)
      save old_time,vhdm_loc
      save init
      data init/0/
      if((init.eq.0).or.(time.ne.old_time))then
         call inter_vel(time,vtimes,nvtimes,i1,i2,vf1,vf2)
         vhdm_loc = vhdivmax(i1) * vf1 + vhdivmax(i2) * vf2
         vhdm = vhdm_loc
         old_time = time
         init=1
      else
         vhdm = vhdm_loc
      endif
      end
c
c
c     determine if a tracer's strain is supposed to be frozen
c
c     returns istatus=1, if so, else 0
c
      subroutine determine_rs_status(time,vtimes,nvtimes,vhdivmax,vgm_f,
     &     xloc,istatus)
#include "fortran_precision.h"
      dimension vtimes(*),vhdivmax(*)
      dimension vgm_f(3,3),xloc(3)
      
      istatus = 0

c     horizontal component of the divergence
      vhdiv = vgm_f(2,2) + vgm_f(3,3) 
c     
c     get the max horizontal divergence at time `time'
c     (this is in spherical system)
c
      call get_vhdivmax(time,vtimes,nvtimes,vhdivmax,vhdm_loc)
c     
c     determine an empirical depth z of the T = 0.5 T_m isoline to check if 
c     the tracer is shallower than that depth. if so, freeze in the strain by
c     removing the symmetric part of the velocity gradient matrix so that 
c     there's only rotations affecting the strain matrix F
      if(xloc(1).gt.0.98901267310195d0)then 
c     tracer is shallower than ~70km
         vhdivnd = vhdiv/vhdm_loc ! normalized horizontal divergence
         if(vhdivnd.gt.0.91047619047619d0)then 
c     strongly positive divergence -> ridge
            rlim = 1.0d0        ! ie., all strains matter, z_T0.5 = 0
         else if(vhdivnd.lt.-0.00527472527472533d0) then 
c     strongly negative, away from ridge -> z_T0.5 = 70km
            rlim = 0.98901267310195d0
         else             
c     calculate depth of 0.5 T_m from an upper bound for z_T05 from 
c     `experimental' data. the approx is z(x) = -21+25/(x+0.28) [km], x is
c     relative horizontal divergence. this is converted into r=(1-z/6371)
            rlim = 1.0032961858421d0 - 0.00392403076440119d0 / 
     &           (vhdivnd + 0.28d0)
         endif
         if(xloc(1).gt.rlim)then 
c     if shallower than critical temperature, return true
            istatus = 1
         endif
      endif
      end
