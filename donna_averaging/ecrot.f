      implicit real*8(a-h,o-z)
      DIMENSION cc(3,3,3,3),ccp(3,3,3,3),a(3,3)
      dimension sp(6,6),cp(6,6),cct(3,3,3,3)
      character *40 fname, fout
	integer R,S
      pi=datan(1.d0)*4.d0
	rpd = pi/180.d0

      open(7,file = 'ecrot.parm')
	read(7,*) iav
	read(7,*) T
	read(7,*) P
	read(7,98) fout

c temperature in degrees Kelvin
c pressure in GPa
	if (T.eq.300.d0 .and. P.eq.0d0) then
	  ioltyp=1
	  iopxtyp=2
	else if (T.eq.1000.d0 .and. P.eq.0d0) then
	  ioltyp=4
	  iopxtyp=9
	else if (T.eq.300.d0 .and. P.eq.7.2d0) then
	  ioltyp=5
	  iopxtyp=8
	else if (T.eq.1500.d0 .and. P.eq.0.d0) then
	  ioltyp=6
	  iopxtyp=9
	  print*,'approximating opx with T=1000K'
	else
	  print*,' unlisted T,P condition input'
	  stop
	endif	

	open(unit=8,file=fout,status='unknown')

	iop = 0
80    if (iop .ne.0) close(1) 
	read(7,98,end=1110) fname
 98   format(a40)
	
      open(1,file=fname)
	iop = 1
100	read(1,*,end=80) nxtl,ngrn2,xpos,zpos
	write(8,*)
	write(8,98)fname
	write(8,109) xpos, zpos
109	format(2f10.3)
C
      do 3 i=1,3
      do 3 j=1,3
      do 3 k=1,3
      do 3 l=1,3
      CCT(i,j,k,l) = 0.0d0
 3    continue

	do 90 ii=1,nxtl
	 nct = nxtl - ii + 1
	 if (ii.eq.1 .and. ngrn2.ne.nxtl) then
         call mineralec(ioltyp,cc,rho)
	 else if (nct.eq.ngrn2) then
         call mineralec(iopxtyp,cc,rho)
	 end if 
       read(1,*) ps,th,ph,wt
       ps = ps*rpd
       th = th*rpd
       ph = ph*rpd
c
      call eulera(ps,th,ph,a)

c
      DO 65 M=1,3
         DO 60 N=1,3
            DO 55 R=1,3
               DO 50 S=1,3
                  CSUM = 0.0
                  DO 45 I=1,3
                     DO 40 J=1,3
                        DO 35 K=1,3
                           DO 30 L=1,3
                             AA = A(M,I)*A(N,J)*A(R,K)*A(S,L)
                             CSUM = CSUM + AA * cc(I,J,K,L)
 30                        CONTINUE
 35                     CONTINUE
 40                  CONTINUE
 45               CONTINUE
c                 IF (ABS(CSUM).LT.10.d0) CSUM=0.d0
                  CCP(M,N,R,S) = CSUM
                  CCT(M,N,R,S) = CCT(M,N,R,S) + CSUM
 50            CONTINUE
 55         CONTINUE
 60      CONTINUE
 65   CONTINUE
c
 90	continue
c
C
      do 70 i=1,3
      do 70 j=1,3
      do 70 k=1,3
      do 70 l=1,3
      ccp(i,j,k,l) = CCT(i,j,k,l)/dble(nxtl)
 70   continue
c
c
      cp(1,1) = ccp(1,1,1,1) 
      cp(1,2) = ccp(1,1,2,2) 
      cp(1,3) = ccp(1,1,3,3) 
      cp(1,4) = ccp(1,1,2,3) 
      cp(1,5) = ccp(1,1,1,3) 
      cp(1,6) = ccp(1,1,1,2) 
      cp(2,2) = ccp(2,2,2,2) 
      cp(2,3) = ccp(2,2,3,3) 
      cp(2,4) = ccp(2,2,2,3) 
      cp(2,5) = ccp(2,2,1,3) 
      cp(2,6) = ccp(2,2,1,2) 
      cp(3,3) = ccp(3,3,3,3) 
      cp(3,4) = ccp(3,3,2,3) 
      cp(3,5) = ccp(3,3,1,3) 
      cp(3,6) = ccp(3,3,1,2) 
      cp(4,4) = ccp(2,3,2,3) 
      cp(4,5) = ccp(2,3,1,3) 
      cp(4,6) = ccp(2,3,1,2) 
      cp(5,5) = ccp(1,3,1,3) 
      cp(5,6) = ccp(1,3,1,2) 
      cp(6,6) = ccp(1,2,1,2) 

	call stiffcomp(cp,sp)
 
      call VRH(iav,cp,sp,rho)

      do i=1,6
       do j=i,6
      write(8,*) i,j,cp(i,j)
      enddo
      enddo

	goto 100
c
1110  close(7)
	close(8)
	stop
      end
c
c.....................................................................
      subroutine eulera(ps,th,ph,c)
c.....................................................................

c **  make the rotation matrix from the three euler angles.

      implicit real*8(a-h,o-z)

      dimension c(3,3)

      sps=dsin(ps)
      cps=dcos(ps)
      sth=dsin(th)
      cth=dcos(th)
      sph=dsin(ph)
      cph=dcos(ph)
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

      SUBROUTINE mineralec(ityp,cc,rhoo)
      implicit real*8 (a-h,o-z)
      dimension c(6,6),ct(6,6)
      dimension co(6,6), cc(3,3,3,3)
C Note : density normalized elastic coeff are returned.
c modified 3/00 by dkb to choose which mineral
c Data from AGU Handbook of Physical Const., T Ahrens (ed)
c p67 & p71, for Temperature=300 K
C
c
	if (ityp.eq.1) then
C  Data for olivine Fo90 (Isaak & Anderson, AGU Hndbk, 1999 @ 300C)
	print*,'ol T0 P0'

	  c(1,1) = 3.206d11
	  c(2,2) = 1.971d11
	  c(3,3) = 2.342d11
	  c(4,4) = 0.6372d11
	  c(5,5) = 0.776d11
	  c(1,2) = 0.698d11
	  c(1,3) = 0.712d11
	  c(2,3) = 0.748d11
	  c(6,6) = 0.7829d11 
        rhoo = 3353.d0

	else if (ityp.eq.2) then
C  Data for bronzite (Kumazawa, JGR 74, p5973, 1969)
	print*,'opx T0 P0'
	c(1,1)=2.3d11
	c(2,2)=1.65d11
	c(3,3)= 2.06d11
	c(4,4)=  0.83d11
	c(5,5)=  0.76d11
	c(6,6)=  0.79d11
	c(1,2)=  0.7d11
 	c(1,3)=  0.57d11
	c(2,3)=  0.5d11
	rhoo = 3335.d0

	else if (ityp.eq.3) then
c  data for forsterite (Anderson & Isaak, AGU Hndbk 1995, @ 300K)

	  c(1,1) = 3.300d11
	  c(2,2) = 2.000d11
	  c(3,3) = 2.360d11
	  c(4,4) = 0.672d11
	  c(5,5) = 0.815d11
	  c(1,2) = 0.662d11
	  c(1,3) = 0.680d11
	  c(2,3) = 0.721d11
	  c(6,6) = 0.812d11 
        rhoo = 3222.d0

	else if (ityp.eq.4) then
C  Data for forsterite @ 1000K (Anderson & Isaak, 1995)
	print*,'ol T1000 P0'
	c(1,1)= 2.931d11
	c(2,2)= 1.753d11
	c(3,3)= 2.104d11
	c(4,4)= 0.5468d11
	c(5,5)= 0.685d11
	c(6,6)= 0.6707d11
	c(1,2)= 0.703d11
 	c(1,3)= 0.647d11
	c(2,3)= 0.618d11
	rhoo = 3275.d0

	else if (ityp.eq.5) then
C  Data for olivine @ 7.2 GPa (Abramson et al., 1997)
	print*,'ol T0 P7.2'
	c(1,1)= 3.671d11
	c(2,2)= 2.348d11
	c(3,3)= 2.743d11
	c(4,4)= 0.767d11
	c(5,5)= 0.890d11
	c(6,6)= 0.930d11
	c(1,2)= 0.929d11
 	c(1,3)= 0.983d11
	c(2,3)= 1.004d11
	rhoo = 3531.d0

	else if (ityp.eq.6) then
C  Data for olivine @ 1500K (Anderson & Isaak, 1995)
	print*,'ol T1500 P0'
	c(1,1)= 2.720d11
	c(2,2)= 1.598d11
	c(3,3)= 1.921d11
	c(4,4)= 0.4857d11
	c(5,5)= 0.622d11
	c(6,6)= 0.5952d11
	c(1,2)= 0.664d11
 	c(1,3)= 0.598d11
	c(2,3)= 0.562d11
	rhoo = 3212.d0

	else if (ityp.eq.7) then
C  Data for olivine @ 3.1 GPa (Abramson et al.,1997)
c   linear interpolation from Table if not listed
	c(1,1)= 3.425d11
	c(2,2)= 2.149d11
	c(3,3)= 2.493d11
	c(4,4)= 0.6857d11
	c(5,5)= 0.820d11
	c(6,6)= 0.8482d11
	c(1,2)= 0.819d11
 	c(1,3)= 0.804d11
	c(2,3)= 0.878d11
	rhoo = 3434.d0

	else if (ityp.eq.8) then
C  Data for opx @ 7.2 GPa (Chai et al., 1997)
c   linear interpolation from Table 1
	print*,'opx T0 P7.2'
	c(1,1)= 2.964d11
	c(2,2)= 2.309d11
	c(3,3)= 2.942d11
	c(4,4)= 0.9347d11
	c(5,5)= 0.847d11
	c(6,6)= 0.9902d11
	c(1,2)= 1.153d11
 	c(1,3)= 1.053d11
	c(2,3)= 0.988d11
	rhoo = 3484.d0

	else if (ityp.eq.9) then
C  Data for opx @ 1000K (extrapolate using Frisillo & Barsch, 1972)
	print*,'opx T1000 P0'
	c(1,1)= 1.94d11
	c(2,2)= 1.42d11
	c(3,3)= 1.81d11
	c(4,4)= 0.73d11
	c(5,5)= 0.66d11
	c(6,6)= 0.70d11
	c(1,2)= 0.62d11
 	c(1,3)= 0.35d11
	c(2,3)= 0.35d11
	rhoo = 2900.d0


	end if
c
c Rotating coordinate axis by 90deg combinations in order to not have to
c  call subroutine ROT over and over.
c
          c11=c(1,1)
          c22=c(2,2)
          c33=c(3,3)
          c44=c(4,4)
          c55=c(5,5)
          c66=c(6,6)
          c12=c(1,2)
          c13=c(1,3)
          c23=c(2,3)

          ct(1,1)=c11
          ct(2,2)=c22
          ct(3,3)=c33
          ct(4,4)=c44
          ct(5,5)=c55
          ct(6,6)=c66
          ct(1,2)=c12
          ct(1,3)=c13
          ct(2,3)=c23

         ct(1,4) = 0.0d0
         ct(1,5) = 0.0d0
         ct(1,6) = 0.0d0
         ct(2,4) = 0.0d0
         ct(2,5) = 0.0d0
         ct(2,6) = 0.0d0
         ct(3,4) = 0.0d0
         ct(3,5) = 0.0d0
         ct(3,6) = 0.0d0
         ct(4,5) = 0.0d0
         ct(4,6) = 0.0d0
         ct(5,6) = 0.0d0
c
c
         do 15 i=1,5
            do 10 j=1+i,6
               ct(j,i) = ct(i,j)
 10         continue
 15      continue
        do 25 i=1,6
            do 20 j=1,6
               co(i,j) = ct(i,j)/rhoo
 20         continue
 25      continue

	  call ec(co,cc)

      return
      end
C
      SUBROUTINE VRH(iav,C,S,RHO)
C Determines elastic properties of a crystal (modified from JMK &
c Dan Raymer subroutines)
c  dkb 5/99  pass iav to select 1: Voigt; 2: Reuss; 3: Hill

      IMPLICIT REAL*8 (A-H,O-Z)
      real*8 C(6,6),S(6,6)
C
c  VOIGT
	if (iav.ne.2) then
       BMOD=(C(1,1)+C(2,2)+C(3,3)+2.d0*(C(1,2)+C(2,3)+C(1,3)))/9.d0
       SMOD=(C(1,1)+C(2,2)+C(3,3) -C(1,2)-C(2,3)-C(1,3) +
     *       3.d0*(C(4,4)+C(5,5)+C(6,6)))/15.d0
	 BMODV = BMOD
	 SMODV = SMOD
	end if

c  REUSS
	if (iav.ne.1) then
       BMOD =1.d0/(S(1,1)+S(2,2)+S(3,3) +2.d0*(S(1,2)+S(2,3)+S(1,3)))
       SMOD =15.d0/(4.d0*(S(1,1)+S(2,2)+S(3,3)-S(1,2)-S(2,3)-S(1,3))
     *     + 3.d0*(S(4,4)+S(5,5)+S(6,6)))
	 BMODR = BMOD
	 SMODR = SMOD
	end if

c  HILL
	if (iav.eq.3) then
       BMOD = (BMODV + BMODR)/2.d0
       SMOD = (SMODV + SMODR)/2.d0
	end if

c  uncomment if want density-normalized values
c       VP = DSQRT((BMOD + 4.d0*SMOD/3.d0)/RHO)
c       VS = DSQRT(SMOD/RHO)

         VP = DSQRT(BMOD + 4.d0*SMOD/3.d0)
         VS = DSQRT(SMOD)

         write(8,19) VP,VS
19	   format(1p,2e12.3)

      RETURN
      END

C
c
      SUBROUTINE ec(c,cc)
      implicit real*8 (a-h,o-z)
      dimension c(6,6),cc(3,3,3,3)
c
c
         do 20 i=1,3
            do 15 j=1,3
               do 10 k=1,3
                  do 5 l=1,3
                     cc(i,j,k,l) = 0.0d0
 5                continue
 10            continue
 15         continue
 20      continue
c
      cc(1,1,1,1) = c(1,1)
      cc(2,2,2,2) = c(2,2)
      cc(3,3,3,3) = c(3,3)
      cc(2,3,2,3) = c(4,4)
      cc(3,2,3,2) =cc(2,3,2,3)
      cc(2,3,3,2) =cc(2,3,2,3)
      cc(3,2,2,3) =cc(2,3,2,3)
      cc(1,3,1,3) = c(5,5)
      cc(3,1,1,3) =cc(1,3,1,3)
      cc(1,3,3,1) =cc(1,3,1,3)
      cc(3,1,3,1) =cc(1,3,1,3)
      cc(1,1,2,2) = c(1,2)
      cc(2,2,1,1) =cc(1,1,2,2)
      cc(1,1,3,3) = c(1,3)
      cc(3,3,1,1) =cc(1,1,3,3)
      cc(1,1,2,3) = c(1,4)
      cc(1,1,3,2) =cc(1,1,2,3)
      cc(2,3,1,1) =cc(1,1,2,3)
      cc(3,2,1,1) =cc(1,1,2,3)
      cc(1,1,1,3) = c(1,5)
      cc(1,1,3,1) =cc(1,1,1,3)
      cc(1,3,1,1) =cc(1,1,1,3)
      cc(3,1,1,1) =cc(1,1,1,3)
      cc(2,2,3,3) = c(2,3)
      cc(3,3,2,2) =cc(2,2,3,3)
      cc(2,2,2,3) = c(2,4)
      cc(2,2,3,2) =cc(2,2,2,3)
      cc(2,3,2,2) =cc(2,2,2,3)
      cc(3,2,2,2) =cc(2,2,2,3)
      cc(2,2,1,3) = c(2,5)
      cc(2,2,3,1) =cc(2,2,1,3)
      cc(1,3,2,2) =cc(2,2,1,3)
      cc(3,1,2,2) =cc(2,2,1,3)
      cc(2,3,1,2) = c(4,6)
      cc(3,2,1,2) =cc(2,3,1,2)
      cc(1,2,2,3) =cc(2,3,1,2)
      cc(1,2,3,2) =cc(2,3,1,2)
      cc(2,3,2,1) =cc(2,3,1,2)
      cc(3,2,2,1) =cc(2,3,1,2)
      cc(2,1,2,3) =cc(2,3,1,2)
      cc(2,1,3,2) =cc(2,3,1,2)
      cc(1,3,1,2) = c(5,6)
      cc(3,1,1,2) =cc(1,3,1,2)
      cc(1,2,1,3) =cc(1,3,1,2)
      cc(1,2,3,1) =cc(1,3,1,2)
      cc(1,3,2,1) =cc(1,3,1,2)
      cc(3,1,2,1) =cc(1,3,1,2)
      cc(2,1,1,3) =cc(1,3,1,2)
      cc(2,1,3,1) =cc(1,3,1,2)
      cc(1,2,1,2) = c(6,6)
      cc(2,1,1,2) =cc(1,2,1,2)
      cc(1,2,2,1) =cc(1,2,1,2)
      cc(2,1,2,1) =cc(1,2,1,2)
c
c
      return
      end
