!Quasi-Dynamic
! for 3D problem, kernel instability occurs while L/dw > 100

program main

  use problem_class
  use input
  use initialize
  use solver
  use my_mpi, only: init_mpi
  use derivs_all
  use unittests, only: init_tests, kernel_export
  use logger, only : log_screen

  type(problem_type), pointer :: pb
  character(len=32) :: arg = "none"
  character(len=255) :: msg
  allocate(pb)

  odepb => pb

  ! Look for command line arguments
  if (iargc() == 1) then
    call getarg(1, arg)

    if (arg == "test") then
      ! Initiate unit tests
      call init_mpi()
      call init_tests(pb)
      stop
    else if (arg == "kernel") then
      call kernel_export(pb)
      stop
    else
      write(msg, *) "Argument not recognised: ", arg
      call log_screen(msg)
      stop
    endif
  endif

  call init_mpi()
  call read_main(pb)
  !call read_lastox(pb)
  call init_all(pb)
  call solve(pb)

end program main
