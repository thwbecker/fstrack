program main

#define DREX_VNORM(X) (sqrt(dot_product(X,X)))

implicit none   
double precision, dimension(:,:,:),allocatable :: acs
integer :: n,i,n3,irand

irand = 0

write(0,*)'reading n**1/3 from stdin...'
read(*,*)n3

n=n3**3;
allocate(acs(n,3,3))
write(0,*)'writing ',n,' dir cosines'
call drex_init_acs_random_ftrn(n,n3,acs,irand)

do i=1,n
   print *,acs(i,1,1),acs(i,1,2),acs(i,1,3)
enddo


deallocate(acs)

end program main


