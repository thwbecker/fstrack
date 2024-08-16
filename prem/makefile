BDIR = ../bin/$(ARCH)/

INC = -I$(HC_HOME)/

DFLAGS = -DPREM_MODEL_FILE=\"$(HC_HOME)/prem/prem.dat\"

all: readprem premtable readprem_z premtable_z prem_compute_gp

clean:
	rm $(BDIR)/readprem $(BDIR)/premtable

SRCS = readprem.c $(HC_HOME)/prem_util.c 

readprem: $(SRCS)
	$(CC) $(INC) $(CFLAGS) $(DFLAGS) $(SRCS) -o $(BDIR)/readprem -lm

prem_compute_gp: prem_compute_gp.c  $(HC_HOME)/prem_util.c
	$(CC) $(INC) $(CFLAGS) $(DFLAGS) prem_compute_gp.c \
	$(HC_HOME)/prem_util.c -o $(BDIR)/prem_compute_gp -lm

premtable: $(SRCS)
	$(CC) $(INC) $(CFLAGS) $(DFLAGS) $(SRCS) -DPRINT_TABLE -o $(BDIR)/premtable -lm

readprem_z: $(SRCS)
	$(CC) $(INC) $(CFLAGS) $(DFLAGS) $(SRCS) -DUSE_DEPTH -o $(BDIR)/readprem_z -lm

premtable_z: $(SRCS)
	$(CC) $(INC) $(CFLAGS) $(DFLAGS) $(SRCS) -DUSE_DEPTH -DPRINT_TABLE -o $(BDIR)/premtable_z -lm
