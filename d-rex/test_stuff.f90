program test_stuff

implicit none

double precision, dimension (3,3):: l
double precision, dimension (3)::isa
double precision :: gol

gol = 0.
l(1,1)=1;l(1,2)=2;l(1,3)=3;
l(2,1)=4;l(2,2)=5;l(2,3)=6;
l(3,1)=7;l(3,2)=8;l(3,3)=9;


call drex_isacalc_ftrn(isa,gol,l)

print *,isa
print *,gol

end program test_stuff
