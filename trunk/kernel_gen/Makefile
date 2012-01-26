EXEC1 = $(HOME)/qdyn_svn/trunk/kernel_gen/kernel_gen_3D
EXEC2 = $(HOME)/qdyn_svn/trunk/kernel_gen/kernel_gen_3D_c

F90 = ifort

OPT = -O0 -C 

OBJS1 = calc_kernel.f90 dc3d_notag.f kernel_gen_3D.f90
OBJS2 = calc_kernel.f90 dc3d_notag.f kernel_gen_3D_c.f90


kernel_gen_3D: $(OBJS1)
	$(F90) $(OPT) -o $(EXEC1)  $(OBJS1)

kernel_gen_3D_c: $(OBJS2)
	$(F90) $(OPT) -o $(EXEC2)  $(OBJS2)


clean:
	rm $(EXEC1)  $(EXEC2) $(EXEC3) 
