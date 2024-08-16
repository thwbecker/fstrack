c --- polefigew.f  
c ---  program to plot single axis pole figures using input Euler angles
c ---  for olivine axis orientation (relative to global ref).
c ---  Plots pole figs for several locations based on parm
c ---  file listing of filename and position. 
c ---  Link with pssub_dbl.f for postscript routines

c  modify polefigew to deal with becker x,y,z positions

	implicit double precision (a-h)
	implicit double precision (o-z)
	dimension c(3,3)
	character*30 fparm, fin, fout, titl

        pi=datan(1.d0)*4.d0
	d2r = pi/180.d0
	

	fparm='polefigtb.parm'

	open(unit=7,file=fparm, status='old')
	read(7,9) titl
	read(7,9) fout
	read(7,*) xmn,xmx,ymin,ymax
	read(7,*) scl
	read(7,*) iax
	read(7,*) ipagesize,figsize, symbsize, iccl
	read(7,*) ve
	read(7,*) dist
9	format(a30)
	ymn = ve*ymin*scl
	ymx = ve*ymax*scl
	xmn = xmn*scl
	xmx = xmx*scl
	xb = 0.2d0*(xmx-xmn)
	yb = 0.2d0*(ymx-ymn)
	dx = (xmx-xmn)/10.d0
	dy = (ymx-ymn)/10.d0

	print*,'    ----->    POLE AXIS #',iax

	call initps(xmn-xb,xmx+xb,ymn-5.d0*yb,ymx+yb,ipagesize,fout)
      call psaxes(xmn,xmx,dx,ymn,ymx,dy,dx/10.d0,ve)

	call plotps(xmn,ymx+yb/2.d0,3)
      call lineps(1)
	call ctextps(titl)

c  iplt is how often to plot pole figures (skip iplt)
c500	  read(7,509,end=999) iplt,fin
500	  read(7,9,end=999) fin
	  xold = -1.e5
	  yold = -1.e5
	  ict = 0
509	  format(i1,1x,a30)
	  open(unit=9,file=fin,status='old')
600	  read(9,*,end=899) nnxl, ngrn2, xpos, ypos
	  xpos = xpos*scl
	  ypos = ypos*scl
	  d = dsqrt((xpos-xold)**2.d0+(ypos-yold)**2.d0)

	  ypos = ve*ypos
	  ict = ict + 1
c	  ichk = ict + iplt - 1
c	  rplt = real(ichk/iplt)
c	  ipltchk = iplt*aint(rplt)
c	  if ( ipltchk .ne. ichk .or.

	  if (d.le.dist .or.
     +       xpos.lt.xmn .or. xpos.gt.xmx .or.
     +      ypos.lt.ymn .or. ypos.gt.ymx) then
	    do 700 j=1,nnxl
	      read(9,*) 
700	    continue
	    goto 600
	  end if
	
	  xold = xpos
	  yold = ypos/ve

	  call colorps(9)
	  call circleps(xpos,ypos,figsize,iccl)
	  call colorps(1)

	do 800 j=1,nnxl
	  read(9,*,end=999) psi, theta, phi
	  ps = psi*d2r
	  th = theta*d2r
	  ph = phi*d2r
	  call eulera(ps,th,ph,c)
	  if (nnxl-j+1.eq.ngrn2) then
	    call colorps(4)
	    call lineps(1)
	  end if
	  call pltpole(iax,c,xpos,ypos,figsize,symbsize)
800	continue
	goto 600

899	close (9)
   	goto 500

999 	call printps
	close (7)
	stop
	end

	
	
