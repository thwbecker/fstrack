c
c
c code from Donna Blackman in July 2004
c
c
c       slight changes by TWB
c
c       $Id: pltpole.f,v 1.2 2004/07/23 23:18:19 becker Exp $
c
c...............................................................
	subroutine pltpole(iax,c,xo,yo,R,ss)
c...............................................................
	
	implicit double precision (a-h,o-z)
	dimension c(3,3)
	
c--     1 = dsqrt(c(1,iax)**2+c(2,iax)**2+c(3,iax)**2)

c---    do projection onto x1-x3 plane for equal area net 
c	x1 positive to right on page
c	x2 positive into the page
c	x3 positive up on the page
c---    take end of axis in "upper" half plane (x2 > 0) for proj
	
	if (c(2,iax).lt.0.d0) then
	   c(1,iax) = -c(1,iax)
	   c(2,iax) = -c(2,iax)
	   c(3,iax) = -c(3,iax)
	end if
	
	c1 = c(1,iax)
	c2 = c(2,iax)
	c3 = c(3,iax)
	alpha = dkb_arctan(c3,c1)
	gamma = dacos(c2)
	
	rad = dsqrt(2.d0)*R*dsin(gamma/2.d0)
	xp = xo + rad*dcos(alpha)
	yp = yo + rad*dsin(alpha)
c
c this plots a cross at xp,yp, I think
c	
	call plotps(xp-ss,yp-ss,3)
	call plotps(xp+ss,yp+ss,2)
	call plotps(xp-ss,yp+ss,3)
	call plotps(xp+ss,yp-ss,2)
	
	return
	end

C       ***********************                                              
	double precision FUNCTION dkb_arctan(y,x)                                              
C  **********************************************************           
C  VERSION OF ARCTAN(Y/X) WHICH WILL NOT CHOKE ON (0,0) ENTRY           
C       **********************************************************           
	implicit double precision (a-h)
	implicit double precision (o-z)
	PI  =  3.1415926535898                                            
	PI2 =  1.5707963267949                                            
	IF(X) 2,1,2                                                       
 1	dkb_arctan = dSIGN(PI2,Y)                                              
	RETURN                                                            
 2	dkb_arctan = dATAN(Y/X)                                                
	IF(X) 3,3,4                                                       
 3	dkb_arctan = dbk_arctan + dSIGN(PI,Y)                                      
 4	RETURN                                                            
	END 

	
	
