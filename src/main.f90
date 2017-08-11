!Quasi-Dynamic
! for 3D problem, kernel instability occurs while L/dw > 100

program main

  use problem_class
  use input
  use initialize
  use solver
  use my_mpi, only: init_mpi
  use derivs_all

  type(problem_type), pointer :: pb
  allocate(pb)

  odepb => pb

  call init_mpi()
  call read_main(pb)
  call init_all(pb)
  call solve(pb)

end program main
