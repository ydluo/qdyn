!Quasi-Dynamic
! for 3D problem, kernel instability occurs while L/dw > 100

program main

  use problem_class
  use input
  use mesh
  use initialize
  use calc
  use solve_master

  type(problem_type) :: pb

  call read_main(pb)
  write(6,*) 'Initializing parameters: ...'
  call init_mesh(pb%mesh)
  call init_field(pb)
  call init_kernel(pb%lam,pb%smu,pb%mesh,pb%kernel)
  call solve(pb)

end program main 
