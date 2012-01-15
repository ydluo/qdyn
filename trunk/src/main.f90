!Quasi-Dynamic
! for 3D problem, kernel instability occurs while L/dw > 100

 program main

 use problem_class
 use input
 use initialize
 use solve_master

 type(problem_type) :: pb

 call read_main(pb)
 call init_field(pb)
 call init_kernel(pb)
 call solve(pb)

 end program main 
