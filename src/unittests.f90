module unittests

  use constants, only : SOLVER_TYPE
  use problem_class, only : problem_type
  use mesh, only : mesh_get_size
  use solver, only : init_rk45
  use friction
  ! use friction_cns
  use fault_stress, only : init_kernel, export_kernel

  use unittests_rsf
  ! use unittests_cns
  use unittests_aux

  implicit none
  private

  public :: init_tests, kernel_export

contains

!===============================================================================
! Initiate the unit testing framework
!===============================================================================
subroutine init_tests(pb)
  type(problem_type) :: pb

  pb%test%test_mode = .true.
  pb%test%test_passed = .true.

  write(6,*) ""
  write(6,*) "-----------------------------------------------------------------"
  write(6,*) "Initiating QDYN unit test suite"
  write(6,*) ""

  ! Allocate mesh variables/parameters
  call allocate_mesh(pb)
  call initiate_parameters(pb)

  ! Test kernels
  call test_kernel(pb)

  ! Spring-block configuration
  call initiate_springblock(pb)
  call allocate_mesh(pb)
  call initiate_parameters(pb)
  call init_kernel( pb%lam, pb%smu, pb%mesh, pb%kernel, pb%D, pb%H, &
                    pb%i_sigma_cpl, pb%finite, pb%test%test_mode)

  ! write(6,*) ""
  ! write(6,*) "Test suite initialised successfully"

  ! Test RSF
  call test_rsf_friction(pb)

  ! Test CNS (analytical solutions Chen?)
  ! call test_cns_friction(pb)

  ! Test TP (analytical soln. Noda)

  ! Test stress calculations (ask analytical soln. Pablo)

  ! Test MPI

  write(6,*) ""
  write(6,*) "Finalised QDYN unit test suite"
  write(6,*) "-----------------------------------------------------------------"

  if (pb%test%test_passed) then
    call exit(0)
  else
    call exit(1)
  endif
end subroutine init_tests

!===============================================================================
! Export stress transfer kernels
!===============================================================================
subroutine test_kernel(pb)

  use constants, only : SRC_PATH

  type(problem_type) :: pb
  integer :: i, n, num_passed
  character(len=200) :: filename, dir
  double precision :: kernel_n, e
  logical :: pass, all_pass

  dir = SRC_PATH//"/../test/kernels/"
  all_pass = .true.
  num_passed = 4

  ! Discretisation parameters
  pb%mesh%nn = 1024
  pb%mesh%Lfault = 10e3
  pb%mesh%W = 100 * pb%mesh%Lfault
  pb%mesh%dim = 1
  pb%kernel%kind = pb%mesh%dim + 1

  write(6,*) ""
  write(6,*) "Testing kernels ..."
  write(6,*) ""

  do i = 0, 3
    pb%finite = i
    call init_kernel( pb%lam, pb%smu, pb%mesh, pb%kernel, pb%D, pb%H, &
                      pb%i_sigma_cpl, pb%finite, pb%test%test_mode)

    select case (pb%finite)
    case(0)
      filename = "kernel_infinite.dat"
    case(1)
      filename = "kernel_finite.dat"
    case(2)
      filename = "kernel_infinite_symmetric.dat"
    case(3)
      filename = "kernel_finite_symmetric.dat"
    end select

    pass = .true.

    open(unit=99, file=trim(dir)//filename, action="read")
    do n = 1, pb%kernel%k2f%nnfft
      read(99, *) kernel_n
      e = (pb%kernel%k2f%kernel(n) - kernel_n) / (kernel_n + 1e-15)
      if (abs(e) > 1e-12) then
        pass = .false.
        all_pass = .false.
        num_passed = num_passed - 1
      endif
    end do

    call print_subresult(trim(filename), pass)
  end do

  pb%test%test_passed = pb%test%test_passed .and. all_pass
  write(6, '(A, I0, A, I0, A)') " Kernel tests: ", num_passed, " / 4 passed"
  write(6,*) ""

end subroutine test_kernel

!===============================================================================
! Set the parameters corresponding to a 0-D spring-block configuration
!===============================================================================
subroutine initiate_springblock(pb)
  type(problem_type) :: pb
  ! Spring-block configuration
  pb%mesh%dim = 0
  pb%mesh%nn = 1
  pb%mesh%Lfault = 1
  pb%mesh%W = 1

  write(6,*) " * Spring-block configuration initialised"
end subroutine initiate_springblock

!===============================================================================
! Allocate mesh parameters for both RSF and CNS
!===============================================================================
subroutine allocate_mesh(pb)
  type(problem_type) :: pb
  integer :: n

  n = mesh_get_size(pb%mesh) ! number of nodes in this processor

  ! Allocate general parameters
  if (allocated(pb%tau)) then
    deallocate (  pb%tau, pb%sigma, pb%v, pb%theta, pb%theta2, pb%v_star, &
                  pb%ot%iot, pb%ot%iasp, pb%dc, pb%coh, pb%v_pl , pb%a, pb%b, &
                  pb%v1, pb%v2, pb%mu_star, pb%dtau_dt, pb%slip, pb%theta_star )
  endif
  allocate ( pb%tau(n), pb%sigma(n), pb%v(n), pb%theta(n), pb%theta2(n),  &
             pb%v_star(n), pb%ot%iot(n), pb%ot%iasp(n), &
             pb%dc(n), pb%coh(n), pb%v_pl(n) )

  ! Allocate RSF parameters
  allocate ( pb%a(n), pb%b(n), pb%v1(n), pb%v2(n), pb%mu_star(n))

  ! Allocate CNS parameters
  ! allocate( pb%cns_params%a_tilde(n), pb%cns_params%mu_tilde_star(n), &
  !           pb%cns_params%y_gr_star(n), pb%cns_params%H(n), &
  !           pb%cns_params%phi_c(n), pb%cns_params%phi0(n), &
  !           pb%cns_params%IPS_const_diff(n), pb%cns_params%IPS_const_diss1(n), &
  !           pb%cns_params%IPS_const_diss2(n), pb%cns_params%L(n), &
  !           pb%cns_params%lambda(n) )

  allocate ( pb%dtau_dt(n), pb%slip(n), pb%theta_star(n) )

  write(6,*) " * Mesh variables/parameters allocated"
end subroutine allocate_mesh

!===============================================================================
! Allocate simulation settings and mesh parameters
!===============================================================================
subroutine initiate_parameters(pb)
  type(problem_type) :: pb

  ! Solver settings
  pb%dt_try = 1e-1
  pb%dt_max = 0.d0
  pb%acc = 1e-7

  ! Simulation features
  pb%features%localisation = 0
  pb%features%stress_coupling = 0
  pb%i_sigma_cpl = 0
  pb%features%tp = 0
  pb%neqs = 2
  pb%finite = 1

  ! Periodic loading
  pb%Omper = 0.d0
  pb%Aper = 0.d0
  pb%Tper = 31536000.d0

  ! Damage model
  pb%D = 0.d0
  pb%H = 0.d0

  ! Elastic moduli and radiation damping
  pb%lam = 30e9
  pb%smu = 30e9
  pb%beta = 3e3
  pb%zimpedance = 0.5d0*pb%smu/pb%beta

  !
  pb%coh = 0.d0
  pb%time = 0.d0

  write(6,*) " * Global parameters initialised"

end subroutine initiate_parameters

!===============================================================================
! Export stress transfer kernels
!===============================================================================
subroutine kernel_export(pb)
  type(problem_type) :: pb
  integer :: i

  write(6,*) 'Intializing kernel: ...'

  ! Allocate mesh variables/parameters
  call initiate_parameters(pb)

  pb%test%test_mode = .true.
  pb%mesh%dim = 1
  pb%kernel%kind = pb%mesh%dim + 1

  do i = 0, 3
    pb%finite = i
    call export_kernel( pb%lam, pb%smu, pb%mesh, pb%kernel, pb%D, pb%H, &
                        pb%i_sigma_cpl, pb%finite)
  end do

  call test_kernel(pb)

end subroutine kernel_export

end module unittests
