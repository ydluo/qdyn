#== BEGIN USER SETTINGS =============================================

# 1. Modify next line if you need to place the executable elsewhere:
# yingdi
#EXEC = $(HOME)/qdyn_svn/trunk/src/qdyn 
# pablo
EXEC = $(HOME)/bin/qdyn 

# 2. Uncomment/add/modify the lines corresponding to your compiler:

# GNU 
#F90 = g95
#OPT = 
#F90 = gfortran
#OPT = -O3 -Wall
#OPT = -O3 

#-- Intel Fortran --
F90 = ifort
OPT = -O3 -ip -ipo -unroll -openmp
#-no-inline-factor
#OPT = -O3
# OPT = -O0 -warn all -C
#OPT = -O0 -warn all -parallel -par-report3  -C
#OPT = -O0 -g

#OPT = -O3 -ip -ipo -arch SSE2 -tpp7 -xN # pablo.ethz.ch
#OPT = -O3 -ip -ipo -xB   # lapfabio
#OPT = -O3 -ip -ipo -unroll

#-- Lahey Fortran95 --
#F90 = lf95
#OPT = -O
#OPT = --chk aesux

#-- Digital --
#F90 = f95
#OPT = -fast

#== END USER SETTINGS =============================================

OBJS := $(patsubst %.f90,%.o,$(wildcard *.f90)) 
#OBJS := $($(wildcard *.f90):.f90=.o) 

all: $(EXEC)

$(EXEC): $(OBJS)
	$(F90) $(OPT) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(EXEC) $(OBJS) *.M *.mod 

.SUFFIXES: 
.SUFFIXES: .f90 .o

%.o: %.f90
	$(F90) $(OPT) -c $<

include Makefile.depend
