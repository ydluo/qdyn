module unittests_rsf

  use constants, only : SOLVER_TYPE
  use problem_class, only : problem_type
  use mesh, only : mesh_get_size
  use solver, only : init_rk45
  use friction
  use fault_stress, only : init_kernel
  use unittests_aux

  implicit none
  private

  public :: test_rsf_friction

contains

!===============================================================================
! Instantiate variables for RSF testing
!===============================================================================
subroutine initiate_RSF(pb)
  type(problem_type) :: pb

  pb%a = 0.001
  pb%b = 0.0015
  pb%dc = 1e-5
  pb%v_star = 1e-6
  pb%mu_star = 0.6
  pb%v1 = 1e-3
  pb%v2 = 1e-1

  pb%i_rns_law = 0
  pb%itheta_law = 1
  pb%slip = 0.d0
  pb%v = pb%v_star
  pb%theta = pb%dc / pb%v

  call set_theta_star(pb)
  pb%tau = pb%sigma * friction_mu(pb%v, pb%theta, pb) + pb%coh

  ! Initiate solvers
  call initiate_solver(pb)

  write(6,*) " * RSF model set-up"

end subroutine initiate_RSF

!===============================================================================
! Subroutine to test code units using the rate-and-state friction framework
!===============================================================================
subroutine test_rsf_friction(pb)

  type(problem_type) :: pb
  double precision, dimension(pb%mesh%nn) :: a, b, dc, v_star, v1, v2, mu_star
  double precision, dimension(pb%mesh%nn) :: dtheta, dv, dtau, mu, mu2
  double precision, dimension(pb%mesh%nn) :: dmu_dv, dmu_dv2
  double precision, dimension(pb%mesh%nn) :: dmu_dtheta, dmu_dtheta2
  double precision, dimension(pb%mesh%nn) :: x, zero, mu_truth, dmu_truth
  double precision :: atol, rtol, randno
  integer :: num_tests, num_passed, i
  logical :: pass, subpass1, subpass2, subpass3

  num_tests = 0
  num_passed = 0
  atol = pb%acc
  rtol = pb%acc
  zero = 0.d0

  write(6,*) ""
  write(6,*) "Testing rate-and-state friction..."
  write(6,*) ""

  call initiate_RSF(pb)

  ! - Initiate solver, etc. (initial values)
  ! - Solve ODE with a = 0 or b = 0
  ! - Compare with analytical solutions
  ! - Solve ODE with a, b > 0
  ! - Compare with interpolated data
  ! - Request output at fixed t, compare with data
  ! - Test both solvers
  ! - V-steps: test steady-state mu/theta

  ! Steady-state tests: ensure that dtheta/dt = 0 for theta = Dc/V
  ! Subtest ageing law
  pb%itheta_law = 1
  call dtheta_dt(pb%v, x, x, pb%theta, x, dtheta, x, pb)
  subpass1 = abs_assert_close(dtheta, zero, atol)

  ! Subtest slip law
  pb%itheta_law = 2
  call dtheta_dt(pb%v, x, x, pb%theta, x, dtheta, x, pb)
  subpass2 = abs_assert_close(dtheta, zero, atol)

  ! Collect results of subtests and print to screen
  pass = subpass1 .and. subpass2
  pb%test%test_passed = pb%test%test_passed .and. pass
  call print_result("Steady-state evolution laws", pass, num_passed, num_tests)
  call print_subresult("Steady-state ageing law", subpass1)
  call print_subresult("Steady-state slip law", subpass2)

  ! Steady-state friction tests: ensure that dmu = (a-b)*log(V1/V0)
  mu = friction_mu(pb%v, pb%theta, pb)
  pb%v = 10*pb%v_star
  pb%theta = pb%dc / pb%v
  mu2 = friction_mu(pb%v, pb%theta, pb)
  dmu_truth = (pb%a - pb%b)*log(10d0)
  pass = abs_assert_close(mu2-mu, dmu_truth, atol)
  pb%test%test_passed = pb%test%test_passed .and. pass
  call print_result("Steady-state friction (classical RSF)", pass, num_passed, num_tests)

  ! Steady-state friction tests for cut-off velocity model
  ! For V1 = V2: mu = mu* - (a-b)*log(V1/V + 1)
  pb%v = pb%v_star
  pb%theta = pb%dc / pb%v
  pb%v1 = pb%v2
  pb%i_rns_law = 1
  call set_theta_star(pb)
  mu = friction_mu(pb%v, pb%theta, pb)
  mu_truth = pb%mu_star - (pb%a - pb%b)*log(pb%v1/pb%v + 1d0)
  subpass1 = abs_assert_close(mu, mu_truth, atol)

  ! For V1 = V2 << V: mu -> mu*
  pb%v = 1e5*pb%v1
  pb%theta = pb%dc / pb%v
  mu = friction_mu(pb%v, pb%theta, pb)
  mu_truth = pb%mu_star
  subpass2 = abs_assert_close(mu, mu_truth, atol)

  ! For V = V2 >> V1: mu -> mu* + b*log(2)
  pb%v1 = 1e-9
  pb%v2 = 1e0
  pb%v = pb%v2
  pb%theta = pb%dc / pb%v
  call set_theta_star(pb)
  mu = friction_mu(pb%v, pb%theta, pb)
  mu_truth = pb%mu_star + pb%b*log(2d0)
  subpass3 = abs_assert_close(mu, mu_truth, atol)

  ! Collect results of subtests and print to screen
  pass = subpass1 .and. subpass2 .and. subpass3
  pb%test%test_passed = pb%test%test_passed .and. pass
  call print_result("Steady-state friction (cut-off velocities)", pass, num_passed, num_tests)
  call print_subresult("V1 = V2", subpass1)
  call print_subresult("V1 = V2 << V", subpass2)
  call print_subresult("V1 << V2 = V", subpass3)

  ! Ensure that the classical rate-and-state friction law (i_rns_law = 0)
  ! is identical to the regularised RSF law (i_rns_law = 2), for various
  ! V and theta (random) -- NOTE that the random seed is set to default,
  ! so that randno is deterministic between consecutive test runs
  pass = .true.
  subpass1 = .true.
  subpass2 = .true.
  subpass3 = .true.
  do i = 1, 1000
    call random_number(randno)
    randno = randno + 0.5

    pb%i_rns_law = 0
    call set_theta_star(pb)
    mu = friction_mu(pb%v, pb%theta*randno, pb)
    call dmu_dv_dtheta(dmu_dv, dmu_dtheta, pb%v, pb%theta*randno, pb)
    pb%i_rns_law = 2
    call set_theta_star(pb)
    mu2 = friction_mu(pb%v, pb%theta*randno, pb)
    call dmu_dv_dtheta(dmu_dv2, dmu_dtheta2, pb%v, pb%theta*randno, pb)
    ! Test if mu, dmu/dv, and dmu/dtheta are ok for all random values of theta
    ! If one of the assertions fails, subpass will become .false.
    subpass1 = subpass1 .and. abs_assert_close(mu, mu2, atol)
    subpass2 = subpass2 .and. rel_assert_close(dmu_dv, dmu_dv2, rtol)
    subpass3 = subpass3 .and. rel_assert_close(dmu_dtheta, dmu_dtheta2, rtol)
  enddo

  pass = subpass1 .and. subpass2 .and. subpass3
  pb%test%test_passed = pb%test%test_passed .and. pass
  call print_result("Compare classical / regularised RSF", pass, num_passed, num_tests)
  call print_subresult("Friction", subpass1)
  call print_subresult("dmu/dV", subpass2)
  call print_subresult("dmu/dtheta", subpass3)

  write(6, '(A, I0, A, I0, A)') " Rate-and-state friction: ", num_passed, " / ", num_tests, " passed"

end subroutine test_rsf_friction

end module unittests_rsf
