c --- polefigew.f  
c ---  program to plot single axis pole figures using input Euler angles
c ---  for olivine axis orientation (relative to global ref).
c ---  Plots pole figs for several locations based on parm
c ---  file listing of filename and position. 
c ---  Link with pssub_dbl.f for postscript routines

c  modify polefigew to deal with becker x,y,z positions

	implicit double precision (a-h)
	implicit double precision (o-z)
	dimension dc(3)
	character*30 fparm, fin, fout, titl

        pi=datan(1.d0)*4.d0
	d2r = pi/180.d0
	
	size=1.0d0
	
	iax=2			! pole axis

	do j=1,100000
! read in direction cosines
	  read(*,*,end=999) dc(1),dc(2),dc(3)
	  call test_pltpole(dc,xp,yp)
	  print *,xp,yp
	enddo

 999	stop
	end

	
	
