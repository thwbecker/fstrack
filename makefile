#
# makefile for fstrack, a finite strain tracker and splitting/anisotropy tools
# will compile only anisotropy tools by default
#
# $Id: makefile,v 1.17 2006/08/14 23:51:27 becker Exp becker $
#
# requirements: GMT (for flow part), EISPACK, DREX
#
#
# if any of the subroutines argument lists are modified, run 'make proto'
# else, 'make' should compile all necessary sources
#

ARCH=$(shell uname -m | awk '{print(tolower($$1))}')
-include machine_dependent.$(ARCH)



ODIR = objects/$(ARCH)/
BDIR = bin/$(ARCH)/

#
# GMT and netcdf libs
#GMTLIBS = -L$(GMTHOME)/lib/ -lgmt -lpsl -L$(NETCDFDIR)/lib/ -lnetcdf
GMTLIBS = -L$(GMT4HOME)/lib/ -lgmt -L$(NETCDFDIR)/lib/ -lnetcdf 
#
# EISPACK
EISLIBS = -L../eispack/$(ARCH)/ -lmyeis
#
# DREX stuff
#
DREX_DIR = ../d-rex
# for detecting dependencies
DREX_LIB = $(DREX_DIR)/objects/$(ARCH)/libdrex.a 
DREX_LIBS = -L$(DREX_DIR)/objects/$(ARCH)/ -ldrex
DREX_HDR = $(DREX_DIR)/drex.h $(DREX_DIR)/drex_fconst.h
# debugging version of library
DREX_LIBS_DEBUG = -L$(DREX_DIR)/objects/$(ARCH)/ -ldrex.dbg
DREX_LIBS_DFAST = -L$(DREX_DIR)/objects/$(ARCH)/ -ldrex.dfast
DREX_INCLUDE = -I$(DREX_DIR)/ -DUSE_DREX
#
# Schulte-Pelkum's tensor and splitting stuff
#
VERA_DIR = ./single_layer
VERA_INCLUDES = -I$(VERA_DIR)/
VERA_HDR = $(VERA_DIR)/vera_util.h
VERA_LIB = $(VERA_DIR)/objects/$(ARCH)/libvera.a
VERA_LIBS = -L$(VERA_DIR)/objects/$(ARCH)/ -lvera
VERA_LIBS_DEBUG = -L$(VERA_DIR)/objects/$(ARCH)/ -lvera.dbg
VERA_LIBS_DFAST = -L$(VERA_DIR)/objects/$(ARCH)/ -lvera
#
# PREM STUFF
#
PREM_SRCS = prem_util.c
PREM_OBJS = $(ODIR)/prem_util.o
PREM_DEFINES = -DPREM_MODEL_FILE=\"$(PWD)/prem.dat\"
PREM_INCS = prem.h
#
# ggrd stuff from hc package, here only used to handle scalar grid interpolation within fstrack
# comment this out, if not needed. if left in, will need the hc package as well
#
GGRD_DIR = $(HC_HOME)/
GGRD_INCLUDES = -I$(GGRD_DIR)/ -DUSE_GGRD -DFSTRACK_USE_GGRD
GGRD_HDR = $(GGRD_DIR)/hc.h
GGRD_LIB =  $(GGRD_DIR)/objects/$(ARCH)/libggrd.a 
GGRD_LIBS = -L$(GGRD_DIR)/objects/$(ARCH)/ -lggrd -lhc 
GGRD_LIBS_DBG = -L$(GGRD_DIR)/objects/$(ARCH)/ -lggrd.dbg -lhc.dbg

FSTRACK_OBJS = $(ODIR)/input_para.o $(ODIR)/state_handling.o  \
	$(ODIR)/output.o  $(ODIR)/init.o

FSTRACK_OBJS_DEBUG = $(FSTRACK_OBJS:.o=.dbg.o)
FSTRACK_OBJS_DFAST = $(FSTRACK_OBJS:.o=.dfast.o)

#
# numerical recipes tools, will be included in the linalg library
#
NRTOOLOBJS = $(ODIR)/svd_util.o $(ODIR)/nr_util.o $(ODIR)/datafit_util.o 
#
# linear algebra library
#
LAOBJS = $(NRTOOLOBJS) $(ODIR)/gramschmidt.o $(ODIR)/dgpadm.o $(ODIR)/misc_F.o \
	$(ODIR)/linalg.o $(ODIR)/linalg_lapack.o \
	$(ODIR)/invert3x3c.o $(ODIR)/numrec_svd.o $(ODIR)/calc_strain.o
LAOBJS_DEBUG = $(LAOBJS:.o=.dbg.o)
LAOBJS_DFAST = $(LAOBJS:.o=.dfast.o)
# 
# derivative and Runge Kutta routines for flow
#
DER_OBJS = $(ODIR)/deriv.o $(ODIR)/fse_derivs.o \
	$(ODIR)/velinterpol.o $(ODIR)/weights.o $(ODIR)/spline.o \
	 $(ODIR)/odeint.o 	$(ODIR)/readgrds.o 
DER_OBJS_DEBUG = $(DER_OBJS:.o=.dbg.o)
DER_OBJS_DFAST = $(DER_OBJS:.o=.dfast.o)
#
# cartesian code objects
#
CARTOBJS = $(ODIR)/ellsphere_cart.o $(ODIR)/odeintell_cart.o \
	$(ODIR)/cart_vel_int.o $(ODIR)/cornerflow.o
CARTOBJS_DEBUG = $(CART_OBJS:.o=.dbg.o)
CARTOBJS_DFAST = $(CART_OBJS:.o=.dfast.o)
#
# miscellaneous and input output routines
#
MISCIOOBJS = $(ODIR)/indexx.o $(ODIR)/trig.o \
	$(ODIR)/misc.o $(ODIR)/series_analyze.o \
	$(ODIR)/output_simple.o  $(ODIR)/sens_handling.o \
	$(ODIR)/splitting_util.o

MISCIOOBJS_DEBUG = $(MISCIOOBJS:.o=.dbg.o)
MISCIOOBJS_DFAST = $(MISCIOOBJS:.o=.dfast.o)
#
# advection and strain related routine
#
AOBJS = $(ODIR)/calc_isa.o $(ODIR)/advect.o \
	$(ODIR)/bailoutcrit.o $(ODIR)/calc_divergence.o \
	$(ODIR)/calc_lyapunov.o $(ODIR)/advect.o \
	$(ODIR)/advect_and_search_for_state.o
AOBJS_DEBUG = $(AOBJS:.o=.dbg.o)
AOBJS_DFAST = $(AOBJS:.o=.dfast.o)
#
INCLUDES = -I$(GMT4HOME)/include/ -I$(NETCDFDIR)/include/ \
	$(DREX_INCLUDE) $(VERA_INCLUDES) $(PREM_DEFINES) $(GGRD_INCLUDES)
#
# p,T models
#
PTOBJS = $(ODIR)/pmodel.o $(ODIR)/tmodel.o
PTOBJS_DEBUG = $(PTOBJS:.o=.dbg.o)
PTOBJS_DFAST = $(PTOBJS:.o=.dfast.o)

#
# header files with constants
#
HDR_FILES = fstrack.h filenames.h \
	structures.h trig.h constants.h fortran_precision.h \
	$(DREX_HDR) $(VERA_HDR) $(PREM_INCS)

FLOW_HDR_FILES = $(HDR_FILES) $(GGRD_HDR) fstrack_flow.h

#
#
# list of targets
#
all:  dirs libs tools seis_tools split_tools lpo_tools fstrack

really_all: all all_libs fstrack.dbg stools \
	fstrack.dfast sav2decompose.dbg testing 

# strain-rate extraction tools
#
stools: extract_strain_field 
#
# splitting tools
split_tools: 
	cd prem; make; cd ../;\
	cd single_layer; make; cd ../;\
	cd multi_layer; make ; cd ../
# Menke splitting now doesn't compile on Fedora 28 because of the xdr library or something like that
#	cd menke_splitting; make ; cd ..;

#
# LPO tools
lpo_tools: drex calc_lpo_from_streamline generate_vgm

drex:
	cd ../d-rex/; make ; cd ../fstrack;


#
# general finite strain particle tracker tools
#
tools: lpo_tools average_tracers  average_rphi_tracers \
	cvec2ellipsoid tracerl2cevec sav2afactor make_random_tensors \
	make_var_tensor polvgm2cartvgm seis_tools  split_tools
#
# seismic anisotropy interpretation tools
seis_tools: sav2cijkl sav2splitting fazi2splitstat \
	sav2decompose cijklrotate sav2rotate \
	plot_kernel split_fit
#
# debugging programs
testing: test_stuff cFfromdiscreteG   sav2splitting.dbg calc_cornerflow # cart_tester



LIBS = $(ODIR)/liblinalg.a \
	$(ODIR)/libcart.a $(ODIR)/libmiscio.a \
	$(ODIR)/libpt.a $(DREX_LIB) $(VERA_LIB) 

FLOW_LIBS = $(ODIR)/libder.a $(ODIR)/libadvect.a $(GGRD_LIB)

FLOW_LIBS_DEBUG = $(ODIR)/libder.dbg.a $(ODIR)/libadvect.dbg.a $(GGRD_LIB)

FLOW_LIBS_DFAST = $(ODIR)/libder.dfast.a $(ODIR)/libadvect.dfast.a $(GGRD_LIB)

LIBS_DEBUG = $(LIBS:.a=.dbg.a)

LIBS_DFAST = $(LIBS:.a=.dfast.a)

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

libs: $(LIBS) 

all_libs: $(LIBS_DEBUG) $(LIBS_DFAST) $(FLOW_LIBS) $(FLOW_LIBS_DEBUG) $(FLOW_LIBS_DFAST)

#
# libraries
#
$(ODIR)/liblinalg.a: $(LAOBJS)
	$(AR) rv $(ODIR)/liblinalg.a $(LAOBJS)
$(ODIR)/libder.a: $(DER_OBJS)
	$(AR) rv $(ODIR)/libder.a $(DER_OBJS)
$(ODIR)/libcart.a: $(CARTOBJS)
	$(AR) rv $(ODIR)/libcart.a $(CARTOBJS)
$(ODIR)/libmiscio.a: $(MISCIOOBJS)
	$(AR) rv $(ODIR)/libmiscio.a $(MISCIOOBJS)
$(ODIR)/libadvect.a: $(AOBJS)
	$(AR) rv $(ODIR)/libadvect.a $(AOBJS)
$(ODIR)/libpt.a: $(PTOBJS)
	$(AR) rv $(ODIR)/libpt.a $(PTOBJS)
# debugging versions
$(ODIR)/liblinalg.dbg.a: $(LAOBJS_DEBUG)
	$(AR) rv $(ODIR)/liblinalg.dbg.a $(LAOBJS_DEBUG)
$(ODIR)/libder.dbg.a: $(DER_OBJS_DEBUG)
	$(AR) rv $(ODIR)/libder.dbg.a $(DER_OBJS_DEBUG)
$(ODIR)/libcart.dbg.a: $(CARTOBJS_DEBUG)
	$(AR) rv $(ODIR)/libcart.dbg.a $(CARTOBJS_DEBUG)
$(ODIR)/libmiscio.dbg.a: $(MISCIOOBJS_DEBUG)
	$(AR) rv $(ODIR)/libmiscio.dbg.a $(MISCIOOBJS_DEBUG)
$(ODIR)/libadvect.dbg.a: $(AOBJS_DEBUG)
	$(AR) rv $(ODIR)/libadvect.dbg.a $(AOBJS_DEBUG)
$(ODIR)/libpt.dbg.a: $(PTOBJS_DEBUG)
	$(AR) rv $(ODIR)/libpt.dbg.a $(PTOBJS_DEBUG)
# fast debug versions
$(ODIR)/liblinalg.dfast.a: $(LAOBJS_DFAST)
	$(AR) rv $(ODIR)/liblinalg.dfast.a $(LAOBJS_DFAST)
$(ODIR)/libder.dfast.a: $(DER_OBJS_DFAST)
	$(AR) rv $(ODIR)/libder.dfast.a $(DER_OBJS_DFAST)
$(ODIR)/libcart.dfast.a: $(CARTOBJS_DFAST)
	$(AR) rv $(ODIR)/libcart.dfast.a $(CARTOBJS_DFAST)
$(ODIR)/libmiscio.dfast.a: $(MISCIOOBJS_DFAST)
	$(AR) rv $(ODIR)/libmiscio.dfast.a $(MISCIOOBJS_DFAST)
$(ODIR)/libadvect.dfast.a: $(AOBJS_DFAST)
	$(AR) rv $(ODIR)/libadvect.dfast.a $(AOBJS_DFAST)
$(ODIR)/libpt.dfast.a: $(PTOBJS_DFAST)
	$(AR) rv $(ODIR)/libpt.dfast.a $(PTOBJS_DFAST)


$(VERA_LIB):
	cd $(VERA_DIR); \
	make ; \
	cd ..;

$(GGRD_LIB):
	cd $(GGRD_DIR); \
	make ; \
	cd ..;


#
#
# binaries
#
#
fstrack: $(LIBS)  $(FSTRACK_OBJS) $(ODIR)/main.o $(FLOW_LIBS)
	$(CC) -o $(BDIR)/fstrack $(FSTRACK_OBJS) $(ODIR)/main.o \
	-L$(ODIR)/  -lpt -ladvect   -lder  -lpt  -lmiscio \
	-llinalg $(DREX_LIBS) $(GGRD_LIBS) \
	$(GMTLIBS) $(EISLIBS)  $(MATHLIBS) $(FTRN_LIB) \
	$(LDFLAGS) 

calc_lpo_from_streamline: $(LIBS) $(ODIR)/calc_lpo_from_streamline.o
	$(CC) 	$(ODIR)/calc_lpo_from_streamline.o \
	-o $(BDIR)/calc_lpo_from_streamline \
	-L$(ODIR)/   -lpt -lmiscio \
	  -llinalg -lmiscio $(DREX_LIBS)  $(EISLIBS)  $(MATHLIBS)  $(DREX_LIBS) $(FTRN_LIB) $(LDFLAGS)


generate_vgm: $(LIBS) $(ODIR)/generate_vgm.o
	$(CC) $(ODIR)/generate_vgm.o -o $(BDIR)/generate_vgm \
	-L$(ODIR)/  $(DREX_LIBS)  -lpt -lmiscio -llinalg -lmiscio  \
	$(DREX_LIBS)  $(EISLIBS)  $(MATHLIBS) \
	$(FTRN_LIB) $(LDFLAGS) 

fstrack.dbg: $(LIBS_DEBUG) $(FSTRACK_OBJS_DEBUG) $(ODIR)/main.dbg.o
	$(CC) $(FSTRACK_OBJS_DEBUG) $(ODIR)/main.dbg.o -o $(BDIR)/fstrack.dbg \
	-L$(ODIR)/    -ladvect.dbg    -lder.dbg -lpt.dbg  -lmiscio.dbg   \
	-llinalg.dbg $(DREX_LIBS_DEBUG) $(GGRD_LIBS) \
	$(GMTLIBS) $(EISLIBS)  $(MATHLIBS) $(FTRN_LIB) $(LDFLAGS) 


fstrack.dfast: $(LIBS_DFAST) $(FSTRACK_OBJS_DFAST) $(ODIR)/main.dfast.o
	$(CC) $(FSTRACK_OBJS_DFAST) $(ODIR)/main.dfast.o -o $(BDIR)/fstrack.dfast \
	-L$(ODIR)/   -ladvect.dfast   -lder.dfast   -lpt.dfast -lmiscio.dfast \
	-llinalg.dfast $(DREX_LIBS_DFAST) $(GGRD_LIBS) \
	$(GMTLIBS) $(EISLIBS)  $(MATHLIBS) \
	$(FTRN_LIB) $(LDFLAGS)

cart_tester: $(LIBS) $(ODIR)/cart_tester.o 
	$(CC) $(ODIR)/cart_tester.o	$(FFLAGS) $(INCLUDES) \
	-o $(BDIR)/cart_tester -L$(ODIR)/    -lcart -ladvect -lder  -lpt \
	-llinalg -lmiscio  $(EISLIBS) $(MATHLIBS) $(FTRN_LIB) $(LDFLAGS) 

calc_cornerflow: $(LIBS) $(ODIR)/calc_cornerflow.o 
	$(CC) $(ODIR)/calc_cornerflow.o	$(FFLAGS) $(INCLUDES) \
	-o $(BDIR)/calc_cornerflow -L$(ODIR)/   -lcart -ladvect -lder  -lpt \
	-llinalg -lmiscio  $(EISLIBS) $(MATHLIBS) $(FTRN_LIB) $(LDFLAGS) 

average_tracers: $(LIBS) $(ODIR)/average_tracers.o 
	$(CC) $(ODIR)/average_tracers.o \
	-o $(BDIR)/average_tracers -L$(ODIR)/  -lpt   \
	-lmiscio -llinalg   $(DREX_LIBS) $(EISLIBS)  $(MATHLIBS) \
	$(FTRN_LIB) $(LDFLAGS) 

average_rphi_tracers: $(LIBS) $(ODIR)/average_rphi_tracers.o 
	$(CC) $(ODIR)/average_rphi_tracers.o \
	-o $(BDIR)/average_rphi_tracers -L$(ODIR)/ -lpt\
	 -lmiscio -llinalg  $(DREX_LIBS) $(EISLIBS)  $(MATHLIBS) $(FTRN_LIB) $(LDFLAGS) 

sav2afactor: $(LIBS) $(ODIR)/sav2afactor.o $(PREM_OBJS)
	$(CC) $(ODIR)/sav2afactor.o \
	-o $(BDIR)/sav2afactor -L$(ODIR)/  -lpt $(PREM_OBJS) \
	 -lmiscio -llinalg  $(DREX_LIBS) $(EISLIBS)  $(MATHLIBS) $(FTRN_LIB) $(LDFLAGS) 

cvec2ellipsoid: $(ODIR)/cvec2ellipsoid.o 
	$(CC) $(ODIR)/cvec2ellipsoid.o  $(INCLUDES) \
			-o $(BDIR)/cvec2ellipsoid -L$(ODIR)/ -llinalg $(FTRN_LIB) $(LDFLAGS) 

tracerl2cevec: $(ODIR)/tracerl2cevec.o $(LIBS)
	$(CC) $(ODIR)/tracerl2cevec.o  $(INCLUDES) \
		-o $(BDIR)/tracerl2cevec -L$(ODIR)/  -lpt   \
	 -lmiscio  -llinalg  -lmiscio  $(EISLIBS)  $(MATHLIBS)\
	 $(FTRN_LIB) $(LDFLAGS)

polvgm2cartvgm: $(ODIR)/polvgm2cartvgm.o $(LIBS)
	$(CC) $(ODIR)/polvgm2cartvgm.o  \
		-o $(BDIR)/polvgm2cartvgm -L$(ODIR)/ \
		 -lpt  -lmiscio  -llinalg  \
		$(LFLAGS) $(EISLIBS)  $(MATHLIBS) \
	$(FTRN_LIB) $(LDFLAGS)

test_stuff: $(LIBS) $(ODIR)/test_stuff.o
	$(CC) $(ODIR)/test_stuff.o $(INCLUDES) $(OBJS) \
		-o $(ODIR)/test_stuff -L$(ODIR)/   \
		 -lmiscio -lder -ladvect -lpt -llinalg   \
		$(GMTLIBS) $(DREX_LIBS) $(EISLIBS)  $(MATHLIBS) \
		$(LDFLAGS)

cFfromdiscreteG: $(LIBS) $(ODIR)/cFfromdiscreteG.o
	$(CC) $(ODIR)/cFfromdiscreteG.o $(INCLUDES) \
		-o $(BDIR)/cFfromdiscreteG -L$(ODIR)/  -ladvect -lpt \
		-llinalg -lmiscio $(GMTLIBS) $(EISLIBS)  $(MATHLIBS) $(LDFLAGS)

extract_strain_field: $(LIBS) $(ODIR)/extract_strain_field.o $(FSTRACK_OBJS)
	$(CC) $(ODIR)/extract_strain_field.o $(INCLUDES) $(FSTRACK_OBJS) \
		-o $(BDIR)/extract_strain_field -L$(ODIR)/  -ladvect  -lder   -lpt \
		 -lmiscio -llinalg  -lmiscio $(GGRD_LIBS) $(GMTLIBS) $(DREX_LIBS) $(EISLIBS)  $(MATHLIBS) \
		$(FTRN_LIB) $(LDFLAGS)

sav2cijkl: $(LIBS) $(ODIR)/sav2cijkl.o 
	$(CC) $(INCLUDES) $(ODIR)/sav2cijkl.o  \
		-o $(BDIR)/sav2cijkl -L$(ODIR)/  -lmiscio  -lpt  \
		-llinalg $(VERA_LIBS) $(DREX_LIBS) $(EISLIBS)  $(MATHLIBS) \
		$(FTRN_LIB) $(LDFLAGS)

sav2rotate: $(LIBS) $(ODIR)/sav2rotate.o 
	$(CC) $(INCLUDES) $(ODIR)/sav2rotate.o  \
		-o $(BDIR)/sav2rotate -L$(ODIR)/   -lpt \
		-llinalg $(VERA_LIBS) $(DREX_LIBS)  -lmiscio -llinalg \
		$(EISLIBS)  $(MATHLIBS) $(FTRN_LIB) $(LDFLAGS)

plot_kernel: $(LIBS) $(ODIR)/plot_kernel.o 
	$(CC) $(INCLUDES) $(ODIR)/plot_kernel.o  \
		-o $(BDIR)/plot_kernel -L$(ODIR)/    -lpt \
		-llinalg $(VERA_LIBS)  -lmiscio 	-llinalg \
		 $(DREX_LIBS) $(EISLIBS)  $(MATHLIBS) $(FTRN_LIB) $(LDFLAGS)

sav2splitting: $(LIBS) $(ODIR)/sav2splitting.o $(PREM_OBJS)
	$(CC) $(INCLUDES) $(ODIR)/sav2splitting.o $(PREM_OBJS) \
		-o $(BDIR)/sav2splitting -L$(ODIR)/  -lmiscio   -lpt    \
		-llinalg $(VERA_LIBS) \
		$(DREX_LIBS) $(EISLIBS)  $(MATHLIBS) \
		$(FTRN_LIB) $(LDFLAGS) $(VERA_LIBS) 

split_fit: $(LIBS) $(ODIR)/split_fit.o $(PREM_OBJS)
	$(CC) $(INCLUDES) $(ODIR)/split_fit.o $(PREM_OBJS) \
		-o $(BDIR)/split_fit -L$(ODIR)/  -lmiscio   -lpt    \
		-llinalg $(VERA_LIBS) \
		$(DREX_LIBS) $(EISLIBS)  $(MATHLIBS) \
		$(FTRN_LIB) $(LDFLAGS) $(VERA_LIBS) 

cijklrotate: $(LIBS) $(ODIR)/cijklrotate.o 
	$(CC) $(INCLUDES) $(ODIR)/cijklrotate.o  \
		-o $(BDIR)/cijklrotate -L$(ODIR)/  -lpt \
		-llinalg -lmiscio $(VERA_LIBS) $(DREX_LIBS) \
		$(EISLIBS)  $(MATHLIBS) \
		$(FTRN_LIB) $(LDFLAGS)

sav2splitting.dbg: $(LIBS) $(ODIR)/sav2splitting.dbg.o  $(PREM_OBJS)
	$(CC) $(INCLUDES) $(ODIR)/sav2splitting.dbg.o  $(PREM_OBJS) \
		-o $(BDIR)/sav2splitting.dbg -L$(ODIR)/  -lmiscio.dbg   -ladvect.dbg  \
		   -lder.dbg  -lpt.dbg \
		-llinalg.dbg $(VERA_LIBS_DEBUG) $(DREX_LIBS_DEBUG) $(CFLAGS_DEBUG) \
		$(GMTLIBS) $(EISLIBS)  $(MATHLIBS) \
		$(FTRN_LIB) $(LDFLAGS)

fazi2splitstat: $(LIBS) $(ODIR)/fazi2splitstat.o 
	$(CC) $(INCLUDES) $(ODIR)/fazi2splitstat.o  \
		-o $(BDIR)/fazi2splitstat -L$(ODIR)/  -lmiscio   -lpt \
		-llinalg  $(VERA_LIBS) $(DREX_LIBS) $(EISLIBS)  $(MATHLIBS) \
		$(FTRN_LIB) $(LDFLAGS)

sav2decompose: $(LIBS) $(ODIR)/sav2decompose.o  $(PREM_OBJS)
	$(CC) $(INCLUDES) $(ODIR)/sav2decompose.o  $(PREM_OBJS) \
		-o $(BDIR)/sav2decompose -L$(ODIR)/    -lpt \
		  $(VERA_LIBS) $(DREX_LIBS) \
		-llinalg  -lmiscio    $(EISLIBS)  $(MATHLIBS) \
		$(FTRN_LIB) $(LDFLAGS)

make_random_tensors: $(LIBS) $(ODIR)/make_random_tensors.o  $(PREM_OBJS)
	$(CC) $(INCLUDES) $(ODIR)/make_random_tensors.o  $(PREM_OBJS) \
		-o $(BDIR)/make_random_tensors -L$(ODIR)/    -lpt \
		$(VERA_LIBS) $(DREX_LIBS) \
		-lpt -lmiscio -llinalg $(EISLIBS)  $(MATHLIBS) \
		$(FTRN_LIB) $(LDFLAGS)

make_var_tensor: $(LIBS) $(ODIR)/make_var_tensor.o  $(PREM_OBJS)
	$(CC) $(INCLUDES) $(ODIR)/make_var_tensor.o  $(PREM_OBJS) \
		-o $(BDIR)/make_var_tensor -L$(ODIR)/     \
		$(VERA_LIBS) $(DREX_LIBS) \
		 -lpt -lmiscio -llinalg  $(EISLIBS)  $(MATHLIBS) \
		$(FTRN_LIB) $(LDFLAGS)

sav2decompose.dbg: $(LIBS) $(ODIR)/sav2decompose.dbg.o $(PREM_OBJS)
	$(CC) $(INCLUDES) $(ODIR)/sav2decompose.dbg.o  $(PREM_OBJS) \
		-o $(BDIR)/sav2decompose.dbg \
		-L$(ODIR)/  -lmiscio   -ladvect     -lder  -lpt \
		-llinalg  $(GGRD_LIBS_DBG) $(VERA_LIBS_DEBUG) \
		$(DREX_LIBS_DEBUG) \
		-llinalg  -lmiscio \
		$(GMTLIBS) $(EISLIBS)  $(MATHLIBS) \
		$(FTRN_LIB) $(LDFLAGS)


driver_test: $(OBJS) $(ODIR)/driver.o
	$(F77) $(ODIR)/driver.o $(FFLAGS) $(INCLUDES) $(OBJS) \
		-o $(ODIR)/driver_test $(LFLAGS)

proto: 
	rm auto_proto.h;\
	cproto  $(INCLUDES) -f2 -q *.c  | grep -v "void main("  | grep -v "int main(" > auto_proto.h

clean:
	rm $(ODIR)/*.o $(ODIR)/*.a $(BDIR)/*
#
# general rules 
#
$(ODIR)/%.o: %.c  $(HDR_FILES)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $(ODIR)/$*.o

$(ODIR)/%.o: %.F $(HDR_FILES)
	$(F77) $(FFLAGS) $(INCLUDES) -c $< -o $(ODIR)/$*.o

$(ODIR)/vera_util.o: $(VERA_DIR)/vera_util.f90 $(HDR_FILES)
	$(F90) $(F90FLAGS) $(INCLUDES) -c $(VERA_DIR)/vera_util.f90 \
	-o $(ODIR)/vera_util.o

#


#
# debugging
#
$(ODIR)/%.dbg.o: %.c  $(HDR_FILES)
	$(CC) -DFSTRACK_DEBUG   \
	$(CFLAGS_DEBUG) $(INCLUDES) -c $< -o $(ODIR)/$*.dbg.o

$(ODIR)/%.dbg.o: %.F  $(HDR_FILES)
	$(F77)  -DFSTRACK_DEBUG $(FFLAGS_DEBUG) $(INCLUDES) -c $< -o $(ODIR)/$*.dbg.o

# debug flags, but fast
$(ODIR)/%.dfast.o: %.c  $(HDR_FILES)
	$(CC)   $(CFLAGS) -DFSTRACK_DEBUG $(INCLUDES) -c $< -o $(ODIR)/$*.dfast.o

$(ODIR)/%.dfast.o: %.F  $(HDR_FILES)
	$(F77)  $(FFLAGS) -DFSTRACK_DEBUG $(INCLUDES) -c $< -o $(ODIR)/$*.dfast.o
