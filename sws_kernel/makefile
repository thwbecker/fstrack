BDIR = $(HOME)/progs/bin/$(ARCH)/

all: $(BDIR)/splitting_total


$(BDIR)/splitting_total: splitting_total.f90
	$(F90) $(F90FLAGS)  splitting_total.f90 -o $(BDIR)/splitting_total
