#== BEGIN USER SETTINGS =============================================

# 1. Modify the path to the executable file only if you really need to place it elsewhere.
#    If you change it, when you call qdyn.m you must set the input EXEC_PATH as set here
EXEC_PATH = .
#EXEC_PATH = $(HOME)/qdyn_svn/trunk/src/
#EXEC_PATH = $(HOME)/bin/

# 2. Comment/uncomment/add/modify the lines corresponding to your compiler:

# GNU, serial
#Only MPI no openMP
F90 = mpif90
OPT = -O3 -Wall
# For parallel runs with MPI+OpenMP:
#F90 = mpif90
#OPT = -O3 -Wall -fopenmp
#OPT = -O2 -w
#-- Intel Fortran --
# F90 = ifort
# OPT = -O3 -ip -ipo
# For parallel runs with OpenMP:
#OPT = -O3 -ip -ipo -parallel -openmp -par-threshold:50
#-no-inline-factor
#OPT = -O3 -pg -g
#OPT = -O3
#OPT = -Og -g -Wall -Wextra -pedantic -fimplicit-none -fcheck=all -fbacktrace
#OPT = -O0 -warn all -C
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

# 3. Set the compiler option that enables preprocessing in your compiler
# ifort
# PREPROC = -fpp
# gfortran
PREPROC = -cpp
#PG, KAUST server
# CFLAGS = -B/usr/lib/x86_64-linux-gnu
#CFLAGS = '' 

travis: F90 = mpif90
	OPT = -O3 -Wall
	PREPROC = -cpp

#== END USER SETTINGS =============================================

EXEC = $(EXEC_PATH)/qdyn

OBJS := $(patsubst %.f90,%.o,$(wildcard *.f90))
#OBJS := $($(wildcard *.f90):.f90=.o)

all: $(EXEC)
travis: $(EXEC)

$(EXEC): $(OBJS)
	$(F90) $(OPT) $(CFLAGS) -o $@ $(OBJS) $(LIBS)
#	$(F90) $(OPT) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(EXEC) $(OBJS) *.M *.mod


.SUFFIXES:
.SUFFIXES: .f90 .o

%.o: %.f90
	$(F90) $(OPT) -c $<

# A special compilation rule for constants.f90:
# Use a preprocessor directive (fpp) to define the path to the source files
# This is needed only to find the kernel files
constants.o: constants.f90
	$(F90) $(OPT) $(PREPROC) -D_FPP_SRC_PATH_="'$(CURDIR)'" -c $<

include Makefile.depend
# DO NOT DELETE
