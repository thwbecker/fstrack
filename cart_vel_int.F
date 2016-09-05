      subroutine cart_vel_int(t,work,workder,nwork,
     &     vx,vy,vz,nx,ny,nz,dx,dy,xymax,zlevels,
     &     vtimes,nvtimes,ifremovestrain,vhdivmax) 
c_______________________________________________________________________
c
c     Note:  passed muti-dimensional arrays are in c-order
c
c_______________________________________________________________________

#include "fortran_precision.h"
      VPRECF vx,vy,vz
      
      dimension vx(*), vy(*), vz(*), zlevels(nz),
     &     vtimes(*),vhdivmax(*)   
      dimension work(nwork),workder(nwork)
      dimension grid(0:3,3),weight(0:3,0:1,3)
      dimension igrid(0:3,3),ider(0:3), dv(1:3,0:3),
     &     xymax(2),dxy(2),vloc(3),nxy(2)
c     debugging
      save dxy,nxny,nxnynz,nxy
      save init
      data init/0/
c     



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
c     do some checks: x_i = (i-1) * dx; i=1....nx
c
         if((xymax(1).le.0.0).or.(xymax(2).le.0.0).or.
     &        (xymax(1).ne.nx*dx).or.(xymax(2).ne.dy*ny))then
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
c
c
c     find the right grid cells for tracers
c     
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
c     pick grid points
         do i=0,3
            igrid(i,j) = ixtracer - 1 + i
            grid(i,j) = (igrid(i,j)-1) * dxy(j)
         end do
         call weights(XTRACER(j),grid(0,j),ngrid,iorder,weight(0,0,j))
c     now make sure grid points exist for velocities
         do i=0,3
            if (igrid(i,j).lt.1) igrid(i,j) = igrid(i,j) + nxy(j)
            if (igrid(i,j).gt.nxy(j)) igrid(i,j) = igrid(i,j) - nxy(j)
         enddo
      enddo
c
c     z component
c      
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
      do m=0, 3           !m=0 -> no derivative
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

      end
