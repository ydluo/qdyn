# Make qdyn
EXEC = $(HOME)/qdyn_nm/qdyn

F95  = ifort
OPT  = -O3 -ip -ipo


OBJS = qdyn.o
SS_Main = fftsg.f90 ode_bs.f90 problem_class.f90 input.f90 initialize.f90 solver_1d.f90 solve_master.f90 main.f90 

        

qdyn:	 
	$(F95) $(OPT) $(SS_Main) -o $(EXEC) 

qdyn.o: $(SS_Main)
	$(F95) $(OPT) -c $(SS_Main)

clean:
	rm -f $(EXEC) $(OBJS) *.M *.mod 
