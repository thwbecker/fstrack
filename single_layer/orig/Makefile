OPT = -O -C 

tensor_phase: tensor_phase.f90 Makefile slowness.f slowness.o
	$(F90) $(OPT) tensor_phase.f90 slowness.o -o tensor_phase

slowness.o: slowness.f Makefile
	$(F90) $(OPT) -c slowness.f -o slowness.o

c6x6to81_dens: c6x6to81_dens.f Makefile
	$(F90) -O c6x6to81_dens.f -o c6x6to81_dens
