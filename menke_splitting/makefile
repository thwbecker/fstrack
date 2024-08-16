BDIR = ../bin/$(ARCH)/
#CFLAGS = -g

INCS = -I/usr/include/tirpc
LIBS = -ltirpc -lm -lc -lnsl

all: 	dirs single_layer 

dirs:
	if [ ! -s ../bin/ ];then\
		mkdir ../bin;\
	fi;\
	if [ ! -s ../bin/$(ARCH)/ ];then \
		mkdir ../bin/$(ARCH);\
	fi;	


single_layer:
	make ah_cross_conv_spectoseis 


really_all:
	make ah_cross_conv_1
	make sac_cross_conv_1
	make sac_test
	make ah_cross_conv_spectoseis
	make ah_cross_conv_1_unk_pol
	make ah_cross_conv_2
	make ah_cross_conv_1_group
	make ah_cross_conv_1t_group
	make ah_cross_conv_2_group
	make ah_est1pulse
	make ah_est2pulse
	make ah_est3pulse
	make ah_est4pulse
	make ah_est5pulse

clean:
	rm -f *.o $(BDIR)/*

ah_cross_conv_1: ah_cross_conv_1.c SHahio.c menke_distaz.c
	$(CC) $(CFLAGS) $(INCS)  -o $(BDIR)/ah_cross_conv_1 ah_cross_conv_1.c SHahio.c menke_distaz.c $(LIBS)

sac_cross_conv_1: sac_cross_conv_1.c SHahio.c menke_distaz.c read_sac.c
	$(CC) $(CFLAGS) $(INCS)  -o $(BDIR)/sac_cross_conv_1 sac_cross_conv_1.c SHahio.c menke_distaz.c read_sac.c \
	-I$(SACAUX)/../include/ \
	-L$(SACAUX)/../lib/ -lsacio $(LIBS)

sac_test: sac_test.c SHahio.c menke_distaz.c read_sac.c
	$(CC) $(CFLAGS) $(INCS)  -o $(BDIR)/sac_test sac_test.c SHahio.c menke_distaz.c read_sac.c -I$(SACAUX)/../include/ \
	-L$(SACAUX)/../lib/ -lsacio $(LIBS)

ah_cross_conv_spectoseis: ah_cross_conv_spectoseis.c  ah_cross_conv_spectoseis_driver.c SHahio.c menke_distaz.c
	$(CC) $(CFLAGS) $(INCS)  -o $(BDIR)/ah_cross_conv_spectoseis ah_cross_conv_spectoseis_driver.c ah_cross_conv_spectoseis.c SHahio.c menke_distaz.c $(LIBS)

ah_cross_conv_1_unk_pol: ah_cross_conv_1_unk_pol.c SHahio.c menke_distaz.c
	$(CC) $(CFLAGS) $(INCS)  -o $(BDIR)/ah_cross_conv_1_unk_pol ah_cross_conv_1_unk_pol.c SHahio.c menke_distaz.c $(LIBS)

ah_cross_conv_2: ah_cross_conv_2.c SHahio.c menke_distaz.c
	$(CC) $(CFLAGS) $(INCS)  -o $(BDIR)/ah_cross_conv_2 ah_cross_conv_2.c SHahio.c menke_distaz.c $(LIBS)

ah_cross_conv_1_group: ah_cross_conv_1_group.c SHahio.c menke_distaz.c
	$(CC) $(CFLAGS) $(INCS)  -o $(BDIR)/ah_cross_conv_1_group ah_cross_conv_1_group.c SHahio.c menke_distaz.c $(LIBS)

ah_cross_conv_1t_group: ah_cross_conv_1t_group.c SHahio.c
	$(CC) $(CFLAGS) $(INCS)  -o $(BDIR)/ah_cross_conv_1t_group ah_cross_conv_1t_group.c SHahio.c $(LIBS)

ah_cross_conv_2_group: ah_cross_conv_2_group.c SHahio.c menke_distaz.c
	$(CC) $(CFLAGS) $(INCS)  -o $(BDIR)/ah_cross_conv_2_group ah_cross_conv_2_group.c SHahio.c menke_distaz.c $(LIBS)

ah_est1pulse: ah_est1pulse.c SHahio.c gauss.c
	$(CC) $(CFLAGS) $(INCS)  -o $(BDIR)/ah_est1pulse ah_est1pulse.c  SHahio.c gauss.c $(LIBS)

ah_est2pulse: ah_est2pulse.c SHahio.c gauss.c
	$(CC) $(CFLAGS) $(INCS)  -o $(BDIR)/ah_est2pulse ah_est2pulse.c  SHahio.c gauss.c $(LIBS)

ah_est3pulse: ah_est3pulse.c SHahio.c gauss.c
	$(CC) $(CFLAGS) $(INCS)  -o $(BDIR)/ah_est3pulse ah_est3pulse.c  SHahio.c gauss.c $(LIBS)

ah_est4pulse: ah_est4pulse.c SHahio.c gauss.c
	$(CC) $(CFLAGS) $(INCS)  -o $(BDIR)/ah_est4pulse ah_est4pulse.c  SHahio.c gauss.c $(LIBS)

ah_est5pulse: ah_est5pulse.c SHahio.c gauss.c
	$(CC) $(CFLAGS) $(INCS)  -o $(BDIR)/ah_est5pulse ah_est5pulse.c  SHahio.c gauss.c $(LIBS)
