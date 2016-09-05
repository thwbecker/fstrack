	program tensor_phase

!	splitting analysis in SKS/SKKS incidence angle range for single tensor

!input: files with cijkl in 81 component form in ascii format
!	assumes cijkl in GPa (axes x=1=E, y=2=N, z=3=up),
!	divided by density in g/cm^3 (sic) - output by c6x6to81_dens
!       names hardwired to depth_050,...,depth_300.cijkl (# is bottom of layer)
!	incidence angle range and increment in degrees: SKS/SKKS 3:2:17
!	azimuth range and increment in degrees: 0:5:355

!output:files with incidence angle, azimuth, fast direction as x and y
!	unit components, fast az, delay time slow-fast qS over 50 km layer thickness (s)
!	given incidence angle 
!       one file per tensor (layer), names hardwired to depth_050,...,depth_300.out

! Feb. 2005	
	
	implicit none
	
	character (len=256)::		tensorfilename,outfilename
	character (len=48)::		scalefactstr,densitystr,minazstr,maxazstr
	character (len=48)::		azintstr,minincstr,maxincstr,incintstr
	
	integer ::			itensor,iout
	integer ::			iargc,iaz,iinc			
	integer ::			fastcomp,slowcomp
	integer ::			i,j,k,l,ilay
	
	
	real ::				az,minaz,maxaz,azint
	real ::				inc,mininc,maxinc,incint
	real ::				scalefact,qh,laythick,sfastlen
	real ::				sfastx,sfasty,dt,density,faz
	real :: 			vsfast_phase,vsslow_phase
	real, dimension(3)::		qinc,slowmag,groupmag
	real, dimension(3,3)::		slowvec,polvec,groupvec
	real, dimension(3,3,3,3)::      c
	real ::				pi
	parameter			(pi=3.1415926535898)
			
!	i/o file numbers
	itensor = 	10
	iout    =       11
		
    	  scalefact = 1
    	  
	  density = 3.53
         
	  minaz = 0
          azint = 5
          maxaz = 360-azint
      
	  mininc = 3
          maxinc = 15
          incint = 2

          laythick = 50
	
! loop over 6 layers
        do ilay = 1,6
        if (ilay == 1) then
           tensorfilename = 'depth_050.cijkl'
           outfilename = 'depth_050.out'
         end if
	 if (ilay == 2) then
	   tensorfilename = 'depth_100.cijkl'
	   outfilename = 'depth_100.out'
	 end if
	 if (ilay == 3) then
	   tensorfilename = 'depth_150.cijkl'
	   outfilename = 'depth_150.out'
	 end if
	 if (ilay == 4) then
	   tensorfilename = 'depth_200.cijkl'
	   outfilename = 'depth_200.out'
	 end if
	 if (ilay == 5) then
	    tensorfilename = 'depth_250.cijkl'
	    outfilename = 'depth_250.out'
	  end if
	  if (ilay == 6) then
	    tensorfilename = 'depth_300.cijkl'
	    outfilename = 'depth_300.out'
	  end if
        	
!	read tensor; scale
        open(itensor,file = tensorfilename)
        read(itensor,*) &
               ((((c(i,j,k,l),i=1,3),j=1,3),k=1,3),l=1,3)
        close(itensor)
                	
! 	open output file for this layer
        open(iout,file = outfilename)
	write(iout,*) '  inc   az    sfastx     sfasty    faz    dt'
	
!	loop over incidence angles
	inc = mininc
	iinc = 0
	incloop: do 
	  if(inc > maxinc) exit incloop
	  
!	  loop over azimuth
	  az = minaz
	  iaz = 0
	  iinc = iinc + 1
	  azloop: do 
	    if(az > maxaz) exit azloop
	    iaz = iaz + 1
	    	  
!	    build unit vector of incident slowness
!	    length q = 1; sin inc = horizontal projection of q
!	    minus signs since we are looking at incoming waves from center
!	    of unit sphere
	    qh = sin(inc*pi/180.)
!	    z component
	    qinc(3) = - cos(inc*pi/180.)
!	    x component
	    qinc(1) = - qh * sin(az*pi/180.)
!	    y component
	    qinc(2) = - qh * cos(az*pi/180.)
	    
!	    solve for phase, group velocities, polarizations
!	    results are sorted as (component 1-3, phase qP, qSV, qSH)
	    call slowness(c,qinc,slowmag,slowvec,polvec,groupmag,groupvec)
	    	    	    
!	    find fast S wave (qSV or qSH)
            if(slowmag(2) < slowmag(3)) then
              fastcomp = 2
              slowcomp = 3
            else
              fastcomp = 3
              slowcomp = 2
            end if
	    
!	    save x and y direction of fast S polarization
            sfastlen = sqrt(polvec(1,fastcomp)**2+polvec(2,fastcomp)**2)
            sfastx = polvec(1,fastcomp)/sfastlen
            sfasty = polvec(2,fastcomp)/sfastlen
	    	    
!	    fast and slow S phase velocities
            vsfast_phase = 1/slowmag(fastcomp)
            vsslow_phase = 1/slowmag(slowcomp)
	    
!	    calculate delay time in layer
            dt = laythick / cos(inc*pi/180) * (1/vsslow_phase - 1/vsfast_phase)
	    
!           fast azimuth
            faz = 180/pi * atan2(sfastx,sfasty)	    
	    	    
!	    write inc, az, sfastx, sfasty, dt(layer) to file 
            write(iout,'(f5.1,1x,f5.1,1x,f10.7,1x,f10.7,1x,f6.1,1x,f5.2)') inc,az,sfastx,sfasty,faz,dt

            az = az + azint
          end do azloop
	  
          inc = inc + incint
	  
        end do incloop
	
        close(iout)
	
         end do 
	
	  
	stop
	end
