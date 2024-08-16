c.....................................................................
      subroutine eulera(ps,th,ph,c)
c.....................................................................

c **  make the rotation matrix from the three euler angles.

      implicit real*8(a-h,o-z)

      dimension c(3,3)

      sps=sin(ps)
      cps=cos(ps)
      sth=sin(th)
      cth=cos(th)
      sph=sin(ph)
      cph=cos(ph)
      c(1,1) =  -sps*sph-cps*cph*cth
      c(2,1) =   cps*sph-sps*cph*cth
      c(3,1) =   cph*sth
      c(1,2) =   cph*sps-sph*cps*cth
      c(2,2) =  -cps*cph-sps*sph*cth
      c(3,2) =  sph*sth
      c(1,3) =  cps*sth
      c(2,3) =  sps*sth
      c(3,3) =  cth

      return
      end
