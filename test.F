      program main
      
      do i=1,10
         call s(i,3,5)
      enddo

      end


      subroutine s(icnt,nx,ny)
      save init
      data init/0/

      print *,'icnt:',icnt,' init:',init,'nxny:',nxny
      if(init.eq.0)then
         nxny=nx*ny
         init=1
      endif
      
      end
