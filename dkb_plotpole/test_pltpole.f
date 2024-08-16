c
c
c code from Donna Blackman in July 2004
c
c this function take not the rotation matrix, but the 
c direction cosine, ie. one column of the matrix
c

c       slight changes by TWB
c
c       $Id: pltpole.f,v 1.2 2004/07/23 23:18:19 becker Exp $
c
c...............................................................
	subroutine test_pltpole(dc,xp,yp)
c...............................................................
	
	implicit double precision (a-h,o-z)
	dimension dc(3)
	
c--     1 = dsqrt(dc(1)**2+dc(2)**2+dc(3)**2)

c---    do projection onto x1-x3 plane for equal area net 
c	x1 positive to right on page
c	x2 positive into the page
c	x3 positive up on the page
c---    take end of axis in "upper" half plane (x2 > 0) for proj
	
	if (dc(2).lt.0.d0) then
	   c1 = -dc(1)
	   c2 = -dc(2)
	   c3 = -dc(3)
	else
	   c1 = dc(1)
	   c2 = dc(2)
	   c3 = dc(3)
	endif


	alpha = dkb_arctan(c3,c1)
	gamma = dacos(c2)
	
	rad = dsqrt(2.d0)*dsin(gamma/2.d0)
	xp = rad*dcos(alpha)
	yp = rad*dsin(alpha)
	
	return
	end


	
	
