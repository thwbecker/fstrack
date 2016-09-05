program cijkl2sav
  ! convert vera's cijkl format to thorsten's sav format

  implicit none 
  double precision, dimension(3,3,3,3) :: c
  double precision, dimension(6,6) :: c6
  double precision :: rho
  character*200 :: rho_string
  integer :: i,j

  ! density scaling
  if(iargc()/=1) then
     print*, 'usage: cijkl2sav rho'
     print*, '       reads cijkl from stdin, writes sav to stdout'
     stop
  end if
  call getarg(1,rho_string)
  read(rho_string,*)rho
  ! read from stdin
  call vera_read_cijkl(c,5)
  ! convert to 6,6
  call vera_c4c2(c,c6)
  ! rescale by density
  do i=1,6
     do j=1,6
        c6(i,j)  = c6(i,j)*rho
     enddo
  enddo

  ! print to stdout
  call vera_print_c2_ftrn(c6)   


end program cijkl2sav
