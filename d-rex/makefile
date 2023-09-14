#
# makefile to make the drex library
#
# set architecture flag
ARCH=$(shell uname -m | awk '{print(tolower($$1))}')

-include machine_dependent.$(ARCH)

# output directoriess
ODIR = objects/$(ARCH)/
BDIR = bin/$(ARCH)/
#
# portable drex library routines
#
DLIB_OBJS = $(ODIR)/drex_deriv.o $(ODIR)/drex_decsym.o \
	$(ODIR)/drex_util.o $(ODIR)/drex_driver.o $(ODIR)/drex_lambert.o \
	$(ODIR)/drex_elast_const.o $(ODIR)/drex_matrix.o 

# modules
MODS = $(ODIR)/drex_nrmod.o
# other stuff
MODS_DEBUG =  $(ODIR)/drex_nrmod.dbg.o
MODS_DFAST =  $(ODIR)/drex_nrmod.dfast.o

# for debugging
DLIB_DEBUG_OBJS = $(DLIB_OBJS:.o=.dbg.o)
DLIB_DFAST_OBJS = $(DLIB_OBJS:.o=.dfast.o)
#
DLIB = $(ODIR)/libdrex.a
DLIB_DEBUG = $(ODIR)/libdrex.dbg.a
DLIB_DFAST = $(ODIR)/libdrex.dfast.a
# header files and such
HDR_FLS = drex.h
#
# scca routine in decmod may use eispack eigensystem routines
# 
EISLIB = -L../eispack/$(ARCH)/ -lmyeis -lm

#
# DREX original
#
# objects for dr_v1 (the new, better version of the original)
# all of these depend on the module comvar
#
#
V1_OBJS = $(ODIR)/init.o $(ODIR)/DRexV1.o $(ODIR)/comvar.o 
#
# objects for dr_cf (outdated)
#CF_OBJS = $(ODIR)/D-RexCF.o $(ODIR)/deriv.o $(ODIR)/init.o

-include machine_dependent.$(ARCH)

all: 	dirs $(MODS) libs dr_v1 testing drex_random_cos

really_all: all all_libs

testing: $(BDIR)/test_stuff  #$(BDIR)/test_stuff_c

libs: $(DLIB) $(ODIR)/libmisc_c.a

all_libs: $(DLIB_DEBUG) $(DLIB_DFAST) 

$(ODIR)/libmisc_c.a: $(ODIR)/misc_c.o
	$(AR) rv $(ODIR)/libmisc_c.a $(ODIR)/misc_c.o

$(DLIB): $(DLIB_OBJS)  $(MODS)
	$(AR) rv $(DLIB) $(MODS) $(DLIB_OBJS) 

$(DLIB_DEBUG): $(DLIB_DEBUG_OBJS) $(MODS_DEBUG)
	$(AR) rv $(DLIB_DEBUG) $(DLIB_DEBUG_OBJS) $(MODS_DEBUG)

$(DLIB_DFAST): $(DLIB_DFAST_OBJS) $(MODS_DFAST)
	$(AR) rv $(DLIB_DFAST) $(DLIB_DFAST_OBJS) $(MODS_DFAST)

dr_v1: $(V1_OBJS) $(DLIB) 
	$(F90) $(V1_OBJS)  -L$(ODIR) -o $(BDIR)/dr_v1 -ldrex \
	$(EISLIB) -lmisc_c

drex_random_cos: $(DLIB) $(ODIR)/drex_random_cos.o
	$(F90)  $(ODIR)/drex_random_cos.o \
	-L$(ODIR) -o $(BDIR)/drex_random_cos  \
	-ldrex $(EISLIB) -lmisc_c


$(BDIR)/test_stuff_c: $(DLIB) $(ODIR)/test_stuff_c.o  \
	$(ODIR)/libmisc_c.a
	$(CC)  -L$(ODIR)/   -ldrex -lmisc_c \
	-o $(BDIR)/test_stuff_c $(LDFLAGS) 

$(BDIR)/test_stuff: $(DLIB) $(ODIR)/test_stuff.o
	$(F90) $(ODIR)/test_stuff.o \
	-L$(ODIR) -o $(BDIR)/test_stuff -ldrex -lmisc_c $(F90LDFLAGS) 


dr_cf: $(CF_OBJS) 
	$(F90) $(CF_OBJS) $(F90FLAGS) -o $(BDIR)/dr_cf

clean:
	\rm $(ODIR)/* $(BDIR)/* *.mod


$(ODIR)/drex_nrmod.o: drex_nrmod.f90  $(HDR_FLS) 
	$(F90) $(F90FLAGS)  -c drex_nrmod.f90 \
	-o $(ODIR)/drex_nrmod.o

$(ODIR)/drex_nrmod.dbg.o: drex_nrmod.f90  $(HDR_FLS) 
	$(F90) $(F90FLAGS_DEBUG) -DDREX_DEBUG  -c drex_nrmod.f90 \
	 -o $(ODIR)/drex_nrmod.dbg.o

$(ODIR)/drex_nrmod.dfast.o: drex_nrmod.f90  $(HDR_FLS) 
	$(F90) $(F90FLAGS) -DDREX_DEBUG -c drex_nrmod.f90 \
	 -o $(ODIR)/drex_nrmod.dfast.o

$(ODIR)/%.o: %.f90  $(HDR_FLS) $(MODS)
	$(F90) $(F90FLAGS) -c $< -o $(ODIR)/$*.o

$(ODIR)/%.o: %.c  $(HDR_FLS) 
	$(CC) $(CFLAGS) -c $< -o $(ODIR)/$*.o
#
# debugging routines
$(ODIR)/%.dbg.o: %.f90  $(HDR_FLS) $(MODS)
	$(F90) $(F90FLAGS_DEBUG) -DDREX_DEBUG  -c $< -o $(ODIR)/$*.dbg.o

$(ODIR)/%.dbg.o: %.c  $(HDR_FLS) 
	$(CC) $(CFLAGS_DEBUG) -DDREX_DEBUG -c $< -o $(ODIR)/$*.dbg.o

$(ODIR)/%.dfast.o: %.f90  $(HDR_FLS) $(MODS)
	$(F90) $(F90FLAGS) -DDREX_DEBUG -c $< -o $(ODIR)/$*.dfast.o

$(ODIR)/%.dfast.o: %.c  $(HDR_FLS) 
	$(CC) $(CFLAGS) -DDREX_DEBUG -c $< -o $(ODIR)/$*.dfast.o

dirs:
	if [ ! -s ./objects/ ]; then\
		mkdir objects;\
	fi;
	if [ ! -s $(ODIR) ];then \
		mkdir $(ODIR);\
	fi;\
	if [ ! -s ./bin/ ];then\
		mkdir bin;\
	fi;\
	if [ ! -s bin/$(ARCH)/ ];then \
		mkdir bin/$(ARCH);\
	fi;	

#
# those need the common block from the original DREX
#
$(ODIR)/init.o: init.f90  $(HDR_FLS) $(ODIR)/comvar.o
	$(F90) $(F90FLAGS) -c $< -o $(ODIR)/init.o

$(ODIR)/DRexV1.o: DRexV1.f90  $(HDR_FLS) $(ODIR)/comvar.o
	$(F90) $(F90FLAGS) -c $< -o $(ODIR)/DRexV1.o

