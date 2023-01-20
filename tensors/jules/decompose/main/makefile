# Compiler
# Intel Fortran
#FOR = ifort
#FOR = ifc
#PRODFLAGS = -automatic -O3
#PRODFLAGS = -B100 -s
#PRODFLAGS = -static -tpp7 -xW -ipo -O3
#PRODFLAGS = -tpp7 -xW -ipo -O3
FOR = $(F90)
PRODFLAGS = $(FFLAGS)

# Directory name
DLIB = ../lib/
# File names
FILE1 = nrmod
FILE2 = tensmod
FILE3 = decmod
F1 = $(DLIB)$(FILE1).f90
F2 = $(DLIB)$(FILE2).f90
F3 = $(DLIB)$(FILE3).f90
PROG = elastic_MM

# Executable file
EXECUT = xexe

# Objects
OBJECTS = $(FILE1).o $(FILE2).o $(FILE3).o $(PROG).o 

# Executable file
$(EXECUT)	: $(OBJECTS)
	$(FOR) $(PRODFLAGS) $(OBJECTS) -o $(EXECUT)
$(FILE1).o	: $(F1)
	$(FOR) $(PRODFLAGS) -c $(F1)
$(FILE2).o	: $(F2)
	$(FOR) $(PRODFLAGS) -c $(F2)
$(FILE3).o	: $(F3) $(F2) $(F1)
	$(FOR) $(PRODFLAGS) -c $(F3)
$(PROG).o	: $(PROG).f90 $(F1) $(F2) $(F3)
	$(FOR) $(PRODFLAGS) -c $(PROG).f90

clean:
	rm -f $(EXE) *.o core.*
