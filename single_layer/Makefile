#
# single layer splitting estimation libraries
#
# Vera Schulte-Pelkum (02/2005)
#
# slight changes by TWB 
#
#
#
# set architecture flag
ARCH=$(shell uname -m | awk '{print(tolower($$1))}')
-include machine_dependent.$(ARCH)


ODIR = objects/$(ARCH)/
BDIR = ../bin/$(ARCH)/

VERA_LIB_OBJS = $(ODIR)/vera_util.o $(ODIR)/vera_jacobi.o 
VERA_LIB_OBJS_DEBUG = $(VERA_LIB_OBJS:.o=.dbg.o)
VERA_LIB_OBJS_DFAST = $(VERA_LIB_OBJS:.o=.dfast.o)

all: 	dir lib tensor_phase  c6x6to81_dens c6x6to81_dens_stream cijkl2sav 

clean:
	rm $(ODIR)/*.o $(ODIR)/*.a 

dir:
	if [ ! -s objects/ ];then \
		mkdir objects/;\
	fi;
	if [ ! -s $(ODIR)/ ];then \
		mkdir $(ODIR)/;\
	fi;

lib: $(ODIR)/libvera.a $(ODIR)/libvera.dbg.a $(ODIR)/libvera.dfast.a

$(ODIR)/libvera.a: $(VERA_LIB_OBJS)
	$(AR) rv $(ODIR)/libvera.a $(VERA_LIB_OBJS)

$(ODIR)/libvera.dbg.a: $(VERA_LIB_OBJS_DEBUG)
	$(AR) rv $(ODIR)/libvera.dbg.a $(VERA_LIB_OBJS_DEBUG)

$(ODIR)/libvera.dfast.a: $(VERA_LIB_OBJS_DFAST)
	$(AR) rv $(ODIR)/libvera.dfast.a $(VERA_LIB_OBJS_DFAST)

tensor_phase: $(ODIR)/tensor_phase.o $(ODIR)/libvera.a
	$(F77) $(ODIR)/tensor_phase.o $(FFLAGS)  -L$(ODIR)/ \
	-o $(BDIR)/tensor_phase -lvera $(LDFLAGS)

c6x6to81_dens: c6x6to81_dens.F $(ODIR)/libvera.a
	$(F77) $(FFLAGS) c6x6to81_dens.F  -L$(ODIR)/ \
		-o $(BDIR)/c6x6to81_dens -lvera $(LDFLAGS)

c6x6to81_dens_stream: c6x6to81_dens.F $(ODIR)/libvera.a
	$(F77) $(FFLAGS) c6x6to81_dens.F -DSTREAM -L$(ODIR)/ \
		-o $(BDIR)/c6x6to81_dens_stream -lvera $(LDFLAGS)

cijkl2sav: cijkl2sav.f90 $(ODIR)/libvera.a
	$(F90) $(F90FLAGS) cijkl2sav.f90  -L$(ODIR)/ \
		-o $(BDIR)/cijkl2sav -lvera $(LDFLAGS)



$(ODIR)/%.o: %.f $(HDR_FILES)
	$(F77) $(FFLAGS) $(INCLUDES) -c $< -o $(ODIR)/$*.o

$(ODIR)/%.o: %.F $(HDR_FILES)
	$(F77) $(FFLAGS) $(INCLUDES) -c $< -o $(ODIR)/$*.o

$(ODIR)/%.o: %.f90 $(HDR_FILES)
	$(F90) $(F90FLAGS) $(INCLUDES) -c $< -o $(ODIR)/$*.o

$(ODIR)/%.dbg.o: %.f  $(HDR_FILES)
	$(F77) $(FFLAGS_DEBUG) $(INCLUDES) -c $< -o $(ODIR)/$*.dbg.o

$(ODIR)/%.dbg.o: %.F  $(HDR_FILES)
	$(F77)  $(FFLAGS_DEBUG) $(INCLUDES) -c $< -o $(ODIR)/$*.dbg.o

$(ODIR)/%.dbg.o: %.f90  $(HDR_FILES)
	$(F90)  $(F90FLAGS_DEBUG) $(INCLUDES) -c $< -o $(ODIR)/$*.dbg.o

$(ODIR)/%.dfast.o: %.f  $(HDR_FILES)
	$(F77)   $(FFLAGS) $(INCLUDES) -c $< -o $(ODIR)/$*.dfast.o

$(ODIR)/%.dfast.o: %.F  $(HDR_FILES)
	$(F77)   $(FFLAGS) $(INCLUDES) -c $< -o $(ODIR)/$*.dfast.o

$(ODIR)/%.dfast.o: %.f90  $(HDR_FILES)
	$(F90)   $(F90FLAGS) $(INCLUDES) -c $< -o $(ODIR)/$*.dfast.o
