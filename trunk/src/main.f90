!Quasi-Dynamic

 program main

 use problem_class
 use input
 use initialize
 use solve_master

 type(problem_type), intent(inout)  :: pb

 call read_main(pb)
 call init_field(pb)
 call init_kernel(pb)
 call solve(pb)

 end program main 
