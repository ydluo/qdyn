EXEC = $(HOME)/bin/qdyn

#F77  = g77 
#OPT = -funroll-loops -O2

F77  = ifort
#OPT  = -g -pg -O3 -ip -ipo  # for profiling
OPT  = -O3 -ip -ipo 
#OPT  = -O3 -ip -ipo -xN -no-prec-div # pablo
#OPT  = -O3 -ip -ipo -xB # lapfabio

OBJS = qdyn.o fftsg.o ode_bs.o

.f.o:
	$(F77) $(OPT) -c $?

qdyn:	$(OBJS) qdyn.h
	$(F77) $(OPT) -o $(EXEC) $(OBJS)

qdyn.o:	qdyn.f qdyn.h
	$(F77) $(OPT) -c qdyn.f

clean:
	rm -f $(EXEC) $(OBJS) *.M *.mod 
