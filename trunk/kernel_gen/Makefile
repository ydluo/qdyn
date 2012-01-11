EXEC = $(HOME)/qdyn_svn/trunk/kernel_gen/kernel_gen_3D

F90 = ifort

OPT = -O0 -C 

OBJS = calc_kernel.f90 dc3d_notag.f kernel_gen_3D.f90 

kernel_gen_3D: $(OBJS)
	$(F90) $(OPT) -o $(EXEC)  $(OBJS)
clean:
	rm $(EXEC) 
