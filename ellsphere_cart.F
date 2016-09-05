c
c $Id: ellsphere_cart.F,v 1.3 2011/04/12 06:19:09 becker Exp $
c
c     CARTESIAN version of ellsphere
c
c
c     runge kutta integration driver, calls odeintell
c     this routine also has ellderivs which determines velocities and their
c     derivatives. ellderivs uses routines in vel_interpol.F
c
c     we assume that the velocities are periodic and given on a grid
c     with x_i = i * dx and i=0 ... nx, i.e on 0, dx, 2dx, ... xmax-dx
c     and xmax = dx*nx
c
      subroutine ellsphere_cart(vx,vy,vz,nx,ny,nz,dx,dy,
     &     zlevels,vtimes,nvtimes,
     &     work,ti,tf,dt,ifellipse,eps,iflyap,renorm,rlyap,
     &     rlbailout,tbailout,ifremovestrain,vhdivmax)
c_______________________________________________________________________
c
c     ti  -> initial time 
c     tf  -> final time   both of these have to be in the right units
c                         for the time interpolation scheme to work
c     dt -> suggested first time step
c     work -> array with x,y,z coordinates and F components, if 
c     ifellipse is set, ie. work is work(3) or work(12)
c     vtimes(nvtimes*2) times at which velocity arrays are given
c     iflyap -> =0 - do not calculate lyapunov exponents
c               =1 calculate lyapunov exponents
c     renorm -> renormalization condition for deformation matrix
c            (if def(i,j)>renorm, renormalize)
c     rlyap(i) -> ith lyapunov exponent of tracer (should be 
c            dimensionalized to 3)
c     rlbailout -> bailout criterion for changes in lyapunov exponents
c     tbailout -> output, time of bailout on exit
c     ifremovestrain: remove the symmtetric part of the velocity gradient
c                    matrix depending on the rules laid out in ellderiv
c     vhdivmax(nvtimes): 
c               max horizontal divergence at surface, needed in ellderiv
c               for the strain removal to decide at what depth the critical 
c               lies, given at each time like vtimes
c_______________________________________________________________________

#include "fortran_precision.h"
      dimension vz(*), vy(*), vx(*), zlevels(*),rlyap(3),
     &     alpha(3),vtimes(*),vhdivmax(*),xymax(2)
      parameter (smaxstep = 1.0)
      dimension work(12)

      save xymax
      save init
      data init/0/

      if(init.eq.0)then
         xymax(1)=nx*dx
         xymax(2)=ny*dy
         init=1
      endif
      if (ifellipse.eq.0) then
c no strain tracing and vel gradient matrix needed
         nwork=3 
         if (iflyap.ne.0) then
            write (6,*) 'You cannot calculate lyapunov exponents'
            write (6,*) 'without calculating ellipses!!'
            write (6,*) 'ifellipse,iflyap',ifellipse,iflyap
            stop
         end if
      else
         nwork=12
      end if
      
      t1 = ti   ! start time
      t2 = tf   ! final time
      step = dt ! time step
c
c     max step is internally limited to smaxstep
c
      if(abs(step) .gt. smaxstep)then
         step=(abs(step)/step) * smaxstep
      endif
      stepmin=0.d0
c
c     integrate to obtain position and deformation matrix at t2
c
      call odeintell_cart(work,nwork,t1,t2,eps,step,stepmin,
     &     vx,vy,vz,nx,ny,nz,dx,dy,xymax,zlevels,vtimes,
     &     nvtimes,iflyap,renorm,rlyap,rlbailout,tbailout,
     &     ifremovestrain,vhdivmax)
      
      if (iflyap.ne.0) then
      !increment lyapunov exponents
         call gramschmidt(work,alpha)
         !and increment lyapunov exponents
         do i=1,3
            rlyap(i)=rlyap(i)+log(alpha(i))
         end do
      end if

c
c     check if z out of bounds
c
      if(XTRACER(3).gt.zlevels(nz))then
         print *,'ellsphere_cart: tracer x_z ',XTRACER(3), 'z_max:',
     &        zlevels(nz)
         stop
      endif
      if(XTRACER(3).lt.zlevels(1))then
         print *,'ellsphere_cart: tracer x_z ',XTRACER(3), 'z_min:',
     &        zlevels(1)
         stop
      endif

 
      return
      end

c
c     calculate time derivatives of work array (output: workder) at time t
c     given the velocity arrays vx,vy,vz (dimensions: nz*nx*ny*nvtimes)
c     specified at times vtimes(nvtimes) and rlevels(nz)
c     
c     the time derivative of work(1...3), the location x, is the velocity at x
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
      subroutine ellderivs_cart(t,work,workder,nwork,
     &     vx,vy,vz,nx,ny,nz,dx,dy,xymax,zlevels,
     &     vtimes,nvtimes,ifremovestrain,vhdivmax) 
c_______________________________________________________________________
c
c     Note:  passed muti-dimensional arrays are in c-order
c
c_______________________________________________________________________

#include "fortran_precision.h"
      dimension vx(*), vy(*), vz(*), zlevels(nz),
     &     vtimes(*),vhdivmax(*)   
      dimension work(nwork),workder(nwork)
      dimension grid(0:3,3),weight(0:3,0:1,3)
      dimension igrid(0:3,3),ider(0:3), dv(1:3,0:3),
     &     vgm(3,3),xymax(2),dxy(2),vloc(3),nxy(2)
      logical icalc_vgm
c     debugging
      save dxy,nxny,nxnynz,nxy
      save init
      data init/0/

      if(init.eq.0)then
c
c     initialize helper variables
c     
         dxy(1) = dx
         dxy(2) = dy

         nxy(1) = nx
         nxy(2) = ny
         nxny = nx * ny
         nxnynz = nxny * nz
         
c
c     do some checks : x_i = (i-1) * dx; i=1....nx; x_i = 0 ... xmax-dx
c
         if((xymax(1).le.0.0).or.(xymax(2).le.0.0).or.
     &        (xymax(1).ne.dx*nx).or.(xymax(2).ne.dy*ny))then
            print *,'ellderivs_cart: dimension error, xmax/ymax:',
     &           xymax(1),xymax(2)
            print *,'ellderivs_cart: nx/ny:',nx,ny
            print *,'ellderivs_cart: dx/dy:',dx,dy
            stop
         endif
         if(nz.lt.4)then
            print *,'ellderivs_cart: need at least four depth levels'
            stop
         endif
         do i=2,nz
            if(zlevels(i).le.zlevels(i-1))then
               print *,'ellderivs_cart: zlevels have to be ascending'
               stop
            endif
         enddo
         if(ifremovestrain.eq.1)then
            print *,'ellderivs_cart: removestrain not implemented'
            stop
            
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

c
c
c     find the right grid cells for tracers
c
      do j=1,2
      !
      ! do x and y, assuming periodic
      !
         do while (XTRACER(j).ge.xymax(j))
            XTRACER(j) = XTRACER(j) - xymax(j)
         end do
         do while (XTRACER(j).lt.0.d0)
            XTRACER(j) = XTRACER(j) + xymax(j)
         end do
         
         ixtracer=nint(XTRACER(j)/dxy(j))+1
c     pick grid points for weights
         do i=0,3
            igrid(i,j) = ixtracer-1+i
            grid(i,j)  = (igrid(i,j)-1) * dxy(j)
         end do
         call weights(XTRACER(j),grid(0,j),ngrid,iorder,weight(0,0,j))
c     now make sure grid points exist for velocities
         do i=0,3
            if (igrid(i,j).lt.1) igrid(i,j) = igrid(i,j) + nxy(j)
            if (igrid(i,j).gt.nxy(j)) igrid(i,j) = igrid(i,j) - nxy(j)
         enddo

      enddo
      !
      ! z component
      !
      ixtracer=0
      i=0
      do while ((ixtracer.eq.0).and.(i.lt.nz))
         i=i+1
         if (XTRACER(3).le.zlevels(i)) ixtracer=i
      end do
      if(ixtracer .eq. 0)then  ! no depth levels found, tracer is above surface
         ixtracer = nz ! assign last layer, x_z should be corrected by the RK routines
      endif
      do i=0,3
         igrid(i,3)=ixtracer-1+i
      end do
      do while (igrid(0,3).le.0)
         do i=0,3
            igrid(i,3)=igrid(i,3)+1
         end do
      end do
      do while (igrid(3,3).gt.nz)
         do i=0,3
            igrid(i,3)=igrid(i,3)-1
         end do
      end do
      do i=0,3
         grid(i,3)=zlevels(igrid(i,3))
      end do
      ngrid=3 
      iorder=1
      call weights(XTRACER(3),grid(0,3),ngrid,iorder,weight(0,0,3))

#ifdef FSTRACK_DEBUG
c     check if all indices are ok
      do i=0,3
         if((igrid(i,1).lt.1).or.(igrid(i,1).gt.nx))then
            print *,'ellsphere_cart: x ',i,'index',igrid(i,1),'error'
            stop
         endif
         if((igrid(i,2).lt.1).or.(igrid(i,2).gt.ny))then
            print *,'ellsphere_cart: y ',i,'index ',igrid(i,2),'error'
            stop
         endif
         if((igrid(i,3).lt.1).or.(igrid(i,3).gt.nz))then
            print *,'ellsphere_cart: z ',i,' index',igrid(i,3),'error'
            stop
         endif
      enddo
#endif
c
c     first calculate velocities only (idindex=0) or vel and  derivatives 
c     of velocity (e.g. v_(x,x)) if needed (idindex=3)
c
      do m=0, idindex           !m=0 -> no derivative
         do i=1,3
            dv(i,m)=0.d0         !m=1,2,3 -> derivative wrt x,y,z
         enddo
         do msub=0,3
            ider(msub)=0
         end do
         ider(m)=1 
         do i=0,3
            do j=0,3
               do k=0,3
                  rnet = weight(i,ider(1),1)*
     &                 weight(j,ider(2),2)*
     &                 weight(k,ider(3),3)
c     
c     the storage scheme is i_z * nxny + i_y *nx + i_x
c     
                  index=igrid(i,1) +        ! i_x
     &                 (igrid(j,2)-1)*nx +  ! i_y
     &                 (igrid(k,3)-1)*nxny  ! i_z
                  call get_velocities(vloc(1),vloc(2),vloc(3),
     &                 vx,vy,vz,vtimes,nvtimes,index,nxnynz,t)
                  do l=1,3
                     dv(l,m) = dv(l,m) + rnet * vloc(l)
                  enddo
               end do
            end do
         end do         
      end do
         

c     repack velocities
      do i=1,3
         workder(i)=dv(i,0)
      enddo

      if(icalc_vgm)then
c
c     build velocity gradient matrix (transpose of grad v)
c
         do i=1,3
            do j=1,3
               vgm(i,j) = dv(i,j)
            enddo
         enddo

c     calculate divergence and adjust vgm matrix so that the divergence
c     is zero (not ideal, but functional)
c
         div=(vgm(1,1)+vgm(2,2)+vgm(3,3))/3.0d0
         do i=1,3
            vgm(i,i)=vgm(i,i)-div
         end do
c
c     remove the strain part of the velocity gradient tensor?
c     
         if(ifremovestrain.eq.1)then 
            print *,'removestrain not implemented'
            stop
         endif
c     
c
c     the time derivative of F is L dot F, repack into workder array
c
         index=3
         do i=1,3
            do j=1,3
               index=index+1
               workder(index) = vgm(i,1) * DEF1(j) + 
     &                          vgm(i,2) * DEF2(j) +
     &                          vgm(i,3) * DEF3(j)
            end do
         end do
      endif

      return
      end

