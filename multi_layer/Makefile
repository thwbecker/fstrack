#
# makefile for anisotropic reflectivity 
#
# from Vera Schulte - Pelkum in March 2005
#
# modification by TWB
#
# $Id: Makefile,v 1.4 2005/10/19 23:22:28 becker Exp becker $
#
# set architecture flag
ARCH=$(shell uname -m | awk '{print(tolower($$1))}')

-include machine_dependent.$(ARCH)

ODIR = objects/$(ARCH)/
BDIR = ../bin/$(ARCH)/
#
# general objects
#
OBJS=    $(ODIR)/aniplanesubs.o \
         $(ODIR)/anieigen.o \
         $(ODIR)/cinvert.o \
         $(ODIR)/utilcmat.o \
         $(ODIR)/utilcijkl.o \
         $(ODIR)/gethex.o \
         $(ODIR)/sgeev.o \
         $(ODIR)/sortpol.o \
         $(ODIR)/utilmat.o \
	../single_layer/$(ODIR)/vera_util.o \
	../single_layer/$(ODIR)/vera_jacobi.o

DBG_OBJS = $(OBJS:.o=.dbg.o)
#
# FFT objects
#
FFT_OBJS =	 $(ODIR)/ffa.o \
	 $(ODIR)/ffaffssubs.o \
	 $(ODIR)/ffs.o \
	 $(ODIR)/getspec.o \
	 $(ODIR)/getts.o
DBG_FFT_OBJS = $(FFT_OBJS:.o=.dbg.o)
#
# anicake objects
#
ANICAKE_LIB_OBJS = $(ODIR)/filter.o $(ODIR)/spectoseis.o $(ODIR)/vera_split.o \
	$(ODIR)/vera_xcorr.o $(ODIR)/ah_cross_conv_spectoseis.o 

DBG_ANICAKE_LIB_OBJS = $(ANICAKE_LIB_OBJS:.o=.dbg.o)


CINCLUDES = -I../menke_splitting/ 

SAC_INCLUDES = -I$(SACAUX)/../include/

HDR_FILES = 

all: 	dirs libs progs

really_all: dirs librs sac_progs progs dbg_libs pbg_progs

libs:	$(ODIR)/libref_fft.a $(ODIR)/libref_general.a $(ODIR)/libanicake.a \


dbg_libs: $(ODIR)/libref_fft.dbg.a $(ODIR)/libref_general.dbg.a $(ODIR)/libanicake.dbg.a 

progs: anicake anicake.dbg spectoseis spectoseis_stream \
	skssplit_xcorr skssplit_xcorr_stream 

sac_progs: sac_vera_split sac_vera_split_dbg

clean:
	rm $(ODIR)/*.o $(ODIR)/*.a 

dirs:
	if [ ! -s objects/ ];then \
		mkdir objects/;\
	fi;
	if [ ! -s $(ODIR)/ ];then \
		mkdir $(ODIR)/;\
	fi;\
	if [ ! -s ../bin/ ];then \
		mkdir ../bin/;\
	fi;
	if [ ! -s $(BDIR)/ ];then \
		mkdir $(BDIR)/;\
	fi;

#
# libraries
#
$(ODIR)/libref_general.a: $(OBJS)
	$(AR) rv $(ODIR)/libref_general.a $(OBJS)

$(ODIR)/libanicake.a: $(ANICAKE_LIB_OBJS)
	$(AR) rv $(ODIR)/libanicake.a $(ANICAKE_LIB_OBJS)


$(ODIR)/libref_fft.a: $(FFT_OBJS)
	$(AR) rv $(ODIR)/libref_fft.a $(FFT_OBJS)


$(ODIR)/libref_general.dbg.a: $(DBG_OBJS)
	$(AR) rv $(ODIR)/libref_general.dbg.a $(DBG_OBJS)

$(ODIR)/libanicake.dbg.a: $(DBG_ANICAKE_LIB_OBJS)
	$(AR) rv $(ODIR)/libanicake.dbg.a $(DBG_ANICAKE_LIB_OBJS)

$(ODIR)/libref_fft.dbg.a: $(DBG_FFT_OBJS)
	$(AR) rv $(ODIR)/libref_fft.dbg.a $(DBG_FFT_OBJS)

#
# programs
#
#
# splitting routines
#
skssplit_xcorr: $(ODIR)/skssplit_xcorr.o $(ODIR)/vera_xcorr.o $(ODIR)/vera_split.o 
	$(F90) $(ODIR)/skssplit_xcorr.o $(ODIR)/vera_xcorr.o $(ODIR)/vera_split.o -o \
	$(BDIR)/skssplit_xcorr $(F90LDFLAGS) $(LDFLAGS)
skssplit_xcorr_stream: $(ODIR)/skssplit_xcorr.stream.o $(ODIR)/vera_xcorr.o  $(ODIR)/vera_split.o 
	$(F90) $(ODIR)/skssplit_xcorr.stream.o $(ODIR)/vera_xcorr.o $(ODIR)/vera_split.o -o \
	$(BDIR)/skssplit_xcorr_stream $(F90LDFLAGS) $(LDFLAGS)
#
# reflectivity code 
#
anicake: $(ODIR)/anicake.o $(ODIR)/libref_general.a $(ODIR)/libanicake.a $(ODIR)/libref_fft.a
	$(F90)  $(ODIR)/anicake.o -o $(BDIR)/anicake \
	-L$(ODIR)/ -lref_general -lanicake -lref_fft  $(F90LDFLAGS) $(LDFLAGS)


anicake.dbg: $(ODIR)/anicake.dbg.o $(ODIR)/libref_general.dbg.a \
	$(ODIR)/libanicake.dbg.a  $(ODIR)/libref_fft.dbg.a
	$(F90)  $(ODIR)/anicake.dbg.o -o $(BDIR)/anicake.dbg \
	-L$(ODIR)/ -lref_general.dbg -lanicake.dbg -lref_fft.dbg  $(F90LDFLAGS) $(LDFLAGS)

#
# seismogram generation
#
# 
# general version
spectoseis: $(ODIR)/spectoseis_driver.o $(ODIR)/filter.o $(ODIR)/libref_fft.a \
	$(ODIR)/spectoseis.o
	$(F90) $(ODIR)/spectoseis.o $(ODIR)/spectoseis_driver.o \
		$(ODIR)/filter.o -o $(BDIR)/spectoseis \
		-L$(ODIR)/ -lref_fft	$(F90LDFLAGS) $(LDFLAGS)
#
# stream I/O version
spectoseis_stream: $(ODIR)/spectoseis_driver.stream.o \
	$(ODIR)/filter.o $(ODIR)/libref_fft.a $(ODIR)/spectoseis.o
	$(F90) $(ODIR)/spectoseis_driver.stream.o $(ODIR)/spectoseis.o \
		$(ODIR)/filter.o -o $(BDIR)/spectoseis_stream \
		-L$(ODIR)/ -lref_fft	$(F90LDFLAGS) $(LDFLAGS)

sac_vera_split: $(ODIR)/vera_xcorr.o $(ODIR)/vera_split.o $(ODIR)/sac_vera_split.o \
	$(ODIR)/menke_distaz.o
	$(F90) $(ODIR)/vera_xcorr.o $(ODIR)/vera_split.o $(ODIR)/sac_vera_split.o \
	-L$(SACAUX)/../lib/ 	$(ODIR)/menke_distaz.o \
	-o  $(BDIR)/sac_vera_split -lsacio $(F90LDFLAGS)  $(LDFLAGS) -lm

sac_vera_split.dbg: $(ODIR)/vera_xcorr.dbg.o $(ODIR)/vera_split.dbg.o $(ODIR)/sac_vera_split.dbg.o \
	$(ODIR)/menke_distaz.dbg.o
	$(F90) $(ODIR)/vera_xcorr.dbg.o $(ODIR)/vera_split.dbg.o \
	$(ODIR)/sac_vera_split.dbg.o \
	-L$(SACAUX)/../lib/ 	$(ODIR)/menke_distaz.o \
	-o  $(BDIR)/sac_vera_split.dbg -lsacio $(F90LDFLAGS)  $(LDFLAGS) -lm

#
# extended source format
#
# general, command line version
$(ODIR)/spectoseis_driver.o: spectoseis_driver.f90 $(HDR_FILES)
	$(F90) $(F90FLAGS) $(INCLUDES) \
	$(F_EXT_SOURCE_FLAG) -c spectoseis_driver.f90 \
	-DUSE_COMMAND_LINE_AND_FILES -o $(ODIR)/spectoseis_driver.o 
#
# streaming version
$(ODIR)/spectoseis_driver.stream.o: spectoseis_driver.f90 $(HDR_FILES)
	$(F90) $(F90FLAGS) $(INCLUDES) $(F_EXT_SOURCE_FLAG) \
	-c spectoseis_driver.f90 -o $(ODIR)/spectoseis_driver.stream.o

$(ODIR)/spectoseis.o: spectoseis.f90 $(HDR_FILES)
	$(F90) $(F90FLAGS) $(INCLUDES) $(F_EXT_SOURCE_FLAG) \
	-c spectoseis.f90 -o $(ODIR)/spectoseis.o

$(ODIR)/spectoseis.dbg.o: spectoseis.f90 $(HDR_FILES)
	$(F90) $(F90FLAGS_DEBUG) -DDEBUG $(INCLUDES) $(F_EXT_SOURCE_FLAG) \
	-c spectoseis.f90 -o $(ODIR)/spectoseis.dbg.o


#
# splitting (stream and regular)

$(ODIR)/skssplit_xcorr.o: skssplit_xcorr_driver.F $(HDR_FILES)
	$(F77) $(FFLAGS) $(INCLUDES) $(F_EXT_SOURCE_FLAG) \
	-c skssplit_xcorr_driver.F \
	-DUSE_COMMAND_LINE_AND_FILES -o $(ODIR)/skssplit_xcorr.o

# streaming version
$(ODIR)/skssplit_xcorr.stream.o: skssplit_xcorr_driver.F $(HDR_FILES)
	$(F77) $(FFLAGS) $(INCLUDES) $(F_EXT_SOURCE_FLAG) \
	-c skssplit_xcorr_driver.F \
	 -o $(ODIR)/skssplit_xcorr.stream.o 

#
# regular
#

$(ODIR)/%.o: %.f $(HDR_FILES)
	$(F77) $(F77FLAGS) $(FTRN_ZERO_INIT_FLAG)\
		$(F_EXT_SOURCE_FLAG) $(INCLUDES)  -c $< -o $(ODIR)/$*.o

$(ODIR)/%.o: %.F $(HDR_FILES)
	$(F77) $(FFLAGS) $(FTRN_ZERO_INIT_FLAG)	$(F_EXT_SOURCE_FLAG) \
	$(INCLUDES)  -c $< -o $(ODIR)/$*.o



$(ODIR)/%.o: %.c $(HDR_FILES)
	$(CC) $(CFLAGS) $(CINCLUDES)  -c $< -o $(ODIR)/$*.o

$(ODIR)/%.dbg.o: %.c $(HDR_FILES)
	$(CC) -DDEBUG $(CFLAGS_DEBUG) $(CINCLUDES)  -c $< -o $(ODIR)/$*.dbg.o



$(ODIR)/%.o: %.f90 $(HDR_FILES)
	$(F90) $(F90FLAGS)  $(FTRN_ZERO_INIT_FLAG) 	$(F_EXT_SOURCE_FLAG) $(INCLUDES) -c $< -o $(ODIR)/$*.o

$(ODIR)/%.dbg.o: %.f  $(HDR_FILES)
	$(F77) -g -DDEBUG $(FFLAGS_DEBUG) 	$(F_EXT_SOURCE_FLAG) \
		$(INCLUDES)  $(FTRN_ZERO_INIT_FLAG) -c $< -o $(ODIR)/$*.dbg.o

$(ODIR)/%.dbg.o: %.F  $(HDR_FILES)
	$(F77) -g -DDEBUG $(FFLAGS_DEBUG) 	$(F_EXT_SOURCE_FLAG) \
		$(INCLUDES)  $(FTRN_ZERO_INIT_FLAG) -c $< -o $(ODIR)/$*.dbg.o

$(ODIR)/%.dbg.o: %.f90  $(HDR_FILES)
	$(F90) -g -DDEBUG $(F90FLAGS_DEBUG) 	$(F_EXT_SOURCE_FLAG) \
		$(INCLUDES)  $(FTRN_ZERO_INIT_FLAG) -c $< -o $(ODIR)/$*.dbg.o
