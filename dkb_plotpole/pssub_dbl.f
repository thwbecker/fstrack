	subroutine initps(x1,x2,y1,y2,isize,fout)
	implicit double precision (a-h)
	implicit double precision (o-z)
	character*30 fout
      common/psblk/ xl,yl

      open(unit=18,file=fout,status='unknown')
      rewind(18)
                        
      write(18,100)
      write(18,101)
      write(18,90)
      write(18,91)               
      write(18,92)
      write(18,102)
      write(18,103)
      write(18,104)
      write(18,105)
      write(18,106)
      write(18,107)
      write(18,108)
      write(18,109)
      write(18,110)
      write(18,111)
      write(18,112)
      write(18,113)               
      write(18,114)
      write(18,115)
      write(18,116)
      write(18,117)
      write(18,118)
      write(18,119)
      write(18,120)
      write(18,121)
      write(18,122)
      write(18,123)               
      write(18,124)
      write(18,130)
      write(18,125)
      write(18,126)               
      write(18,127)
      write(18,128)

  100 format('%!')
  101 format('%% Define macros etc')
  102  format('/L {lineto} def')
  103  format('/D {rlineto} def')
  104  format('/M {sk moveto} def')
  105  format('/G {rmoveto} def')
  106  format('/slw {setlinewidth} def')
  107  format('/np {newpath} def')
  108  format('/sk {stroke} def')
  109  format('/cp {closepath} def')
  110  format('/R {rotate} def')
  111  format('/sh {show} def')
  112  format('/sd {setdash} def',/,'/dpi2point 24 def')
  113  format('% Initialize coordinate system and scales')
114    format('/Helvetica findfont',/,' 30 scalefont setfont')
115    format(2i7,' translate')
116	format('76 dict begin %Colortable dictionary')
117	format('/c0 { 0 0 0 sr} bdef')
118	format('/c1 { 1 1 1 sr} bdef')
119	format('/c2 { 1 0 0 sr} bdef')
120	format('/c3 { 0 1 0 sr} bdef')
121	format('/c4 { 0 0 1 sr} bdef')
122	format('/c5 { 1 1 0 sr} bdef')
123	format('/c6 { 1 0 1 sr} bdef')
124	format('/c7 { 0 1 1 sr} bdef')
130	format('/c8 { 0.5 0.5 0.5 sr} bdef')
125	format('/SO { [] 0 setdash } bdef')
126	format('/DO { [.5 dpi2point mul 4
     + dpi2point mul] 0 setdash } bdef')
127	format('/DA { [6 dpi2point mul] 0 setdash } bdef')
90	format('/bdef {bind def} bind def')
91	format('/ldef {load def} bind def')
92	format('/sr /setrgbcolor ldef')
128	format('0.35 0.35 scale')
129	format('90 rotate')

       xdif = x2 - x1
       ydif = y2 - y1

       xlen = dble(isize)/xdif
       ylen = dble(isize)/ydif
	 xl = dmin1(xlen,ylen)
	 yl = xl

       if (x1.lt.0.) then
         if(x2.lt.0.) then
           ixo = idnint((-x2 + xdif)*xl) + 100
         else
           ixo = idnint((-x1)*xl) + 100
         end if
       else
         ixo = -idnint(x1*xl - 100)
       end if

       if (y1.lt.0.) then
         if(y2.lt.0.) then
           iyo = idnint((-y2 + ydif)*yl) + 175
         else
           iyo = idnint((-y1)*yl) + 175
         end if
       else
         iyo = -idnint(y1*yl + 175)
       end if

       write(18,115) ixo, iyo
                           
       return
       end
       

      subroutine plotps(x,y,k)
	implicit double precision (a-h)
	implicit double precision (o-z)
      common/psblk/ xl,yl

      nx1=idnint(xl*x)
      ny1=idnint(yl*y)
      if (k .eq. 3) write(18,100) nx1,ny1
      if (k .ne. 3) then
        if (k .eq. 2) write(18,102) 0,nx1,ny1
        if (k .eq. -1) write(18,101) 2,0,nx1,ny1
        if (k .eq. 4) write(18,101) 6,0,nx1,ny1
      endif  
             
 100  format(i6,i6,' M')  
 101  format('[',i1,'] ',i1,' sd ',i6,i6,' L')
 102  format('[] ',i2,' sd ',i6,i6,' L')
      return
      end

      subroutine colorps(icval)
	implicit double precision (a-h)
	implicit double precision (o-z)
	character*2 icol(12)
	data icol/'c0','c1','c2','c3','c4',
     + 'c5','c6','c7','c8','SO','DO','DA'/

	write(18,109)
	write(18,100) icol(icval)
109	format('sk')
100	format(a2)

      return
      end
                              
      subroutine psaxes
     +    (xmin,xmax,dx,ymin,ymax,dy,tlen,ve)
	implicit double precision (a-h)
	implicit double precision (o-z)
      common/psblk/ xl,yl

      nx = idnint((xmax-xmin)/dx)+ 1
      ny = idnint((ymax-ymin)/dy)+ 1
      xax = xmin - 3.d0*tlen
      yax = ymin - tlen
	xlab = xax - 0.1d0*(xmax-xmin)
	ylab = yax - 0.1d0*(ymax-ymin)
	call plotps(xlab,ymin,3)
	call ftextps(ymin/ve)
	call plotps(xlab,ymax,3)
	call ftextps(ymax/ve)
	call plotps(xmin-dx,ylab,3)
	call ftextps(xmin)
	call plotps(xmax-dx,ylab,3)
	call ftextps(xmax)

      kp=3
      do 100 ii=1,nx
      x1=xmin + dble(ii-1)*dx
      call plotps(x1,yax,kp)
      kp=2    
      call plotps(x1,yax-tlen,3)
      call plotps(x1,yax+tlen,2)
      call plotps(x1,yax,3)
 100  continue

      kp=3
      do 200 ii=1,ny
      y1=ymin + dble(ii-1)*dy
      call plotps(xax,y1,kp)
      kp=2
      call plotps(xax-tlen,y1,3)
      call plotps(xax+tlen,y1,2)
      call plotps(xax,y1,3)
 200  continue

      return
      end

      subroutine ctextps(string)
      character*20 string
      write(18,114)string
114   format('(',a20,') show')

      return
      end

      subroutine ftextps(xval)   
      implicit double precision (a-h)
      implicit double precision (o-z)
      character*1 num(8)
      integer k(7)

	xnum = xval
	k(5) = 0
                  
      num(1) = ' '
      if (xnum.lt.0.d0) then
        num(1) = '-'
	  xnum = -xnum
      end if
      k(1) = idint(xnum/1.d3)
      k(2) = idint(dmod(xnum,1.d3)/1.d2)
      k(3) = idint(dmod(xnum,1.d2)/1.d1)
      k(4) = dmod(xnum,1.d1)
      frac = (xnum - dble(idnint(xnum)))*1.d2
      k(6) = idint(dmod(frac,1.d2)/1.d1)
      k(7) = dmod(frac,1.d1) 
      do 100 m=1,7
       if (m.eq.5) then
	  num(m+1) = '.'
       elseif (k(m) .eq. 0) then
	  num(m+1)='0'
       elseif (k(m) .eq. 1) then
	  num(m+1)='1'
       elseif (k(m) .eq. 2) then
	  num(m+1)='2'
       elseif (k(m) .eq. 3) then
	  num(m+1)='3'
       elseif (k(m) .eq. 4) then
	  num(m+1)='4'
       elseif (k(m) .eq. 5) then
	  num(m+1)='5'
       elseif (k(m) .eq. 6) then
	  num(m+1)='6'
       elseif (k(m) .eq. 7) then
	  num(m+1)='7'
       elseif (k(m) .eq. 8) then
	  num(m+1)='8'
       elseif (k(m) .eq. 9) then
	  num(m+1)='9'
	 end if
100   continue
      if (xnum .lt. 1.d3) num(2) = ' '
      if (xnum .lt. 1.d2) num(3) = ' '
      if (xnum .lt. 1.d1) num(4) = ' '
      write(18,110) (num(j),j=1,7)
110   format( '(',7a1,') show' )

      return
      end
      
             
      subroutine circleps(xo,yo,r,nn)
	implicit double precision (a-h)
	implicit double precision (o-z)

	dang = 2.d0*3.145926535897932/dble(nn)
	call plotps(xo+r,yo,3)
	do 100 i=1,nn+1
	  ang = dble(i-1)*dang
	  x = xo + r*dcos(ang)
	  y = yo + r*dsin(ang)
	  call plotps(x,y,2)
100	continue

      return
      end

      subroutine lineps(k)
      write(18,100) k
  100 format(i4,' slw')
      return
      end


      subroutine printps
      write(18,100)
100   format('sk showpage')
      close(18)

      return
      end

