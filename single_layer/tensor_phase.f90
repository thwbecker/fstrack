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

! Feb. 2005	by Vera Schulte-Pelkum
! 
! minor modification by TWB
!
! $Id: tensor_phase.f90,v 1.4 2005/07/01 21:54:24 becker Exp becker $
!
	
	implicit none
	
	character (len=256)::		tensorfilename,outfilename
	character (len=48)::		scalefactstr,densitystr,minazstr,maxazstr
	character (len=48)::		azintstr,minincstr,maxincstr,incintstr
	
	integer ::			itensor,iout
	integer ::			iargc,iaz,iinc			
	integer ::			i,j,k,l,ilay
	
	
	double precision ::		az,minaz,maxaz,azint
	double precision ::             inc,mininc,maxinc,incint
	double precision ::             laythick
	double precision ::             dt,faz,sfastx,sfasty
        double precision ::             vsphase(2)
	double precision, dimension(3,3,3,3)::      c
			
!	i/o file numbers
	itensor = 	10
	iout    =       11
		
! limits for azimuth loop 
        minaz = 0d0
        azint = 5d0
        maxaz = 360

! limits for incidence
        mininc = 3d0
        incint = 2d0
        maxinc = 15d0
! layer thickness        
!        laythick = 50d0
	laythick = 350d0

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
           
           !	read scaled tensor
           open(itensor,file = tensorfilename)
           call vera_read_cijkl(c,itensor)
           close(itensor)
           ! 	open output file for this layer
           open(iout,file = outfilename)
           !           write(iout,*) '# inc   az    sfastx     sfasty    faz    dt'
           
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
                 ! compute splitting for a single layer
                 call vera_layer_split_from_tensor_ftrn(c,&
                      inc,az,laythick,sfastx,sfasty,faz,dt,vsphase)
                 !	    write inc, az, sfastx, sfasty, dt(layer) to file 
                 call vera_print_splitting_ftrn(inc,az,sfastx,sfasty,faz,dt,iout)
                 az = az + azint
              end do azloop
	  
              inc = inc + incint
              
           end do incloop
           
        close(iout)
	
     end do
     
     
     stop
   end program tensor_phase
   
