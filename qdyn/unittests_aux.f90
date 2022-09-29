module unittests_aux

  use problem_class, only : problem_type
  use constants, only : SOLVER_TYPE
  use solver, only : init_rk45

  implicit none
  public

contains

!===============================================================================
! Assert whether or not a equals b within a specified absolute tolerance
!===============================================================================
function abs_assert_close(a, b, atol) result(pass)
  double precision, dimension(:), intent(in) :: a, b
  double precision, intent(in) :: atol
  logical :: pass

  if (any(abs(a-b) > atol)) then
    ! At least one element of a and b differ by more than
    ! the specified tolerance. Return false (test failed)
    pass = .false.
  else
    ! All elements of a and b are equal within the specified
    ! tolerance. Return true (test passed)
    pass = .true.
  endif

end function abs_assert_close

!===============================================================================
! Assert whether or not a equals b within a specified relative tolerance
!===============================================================================
function rel_assert_close(a, b, rtol) result(pass)
  double precision, dimension(:), intent(in) :: a, b
  double precision, intent(in) :: rtol
  logical :: pass

  if (any(abs((a-b)/(a+b)) > 2*rtol)) then
    ! At least one element of a and b differ by more than
    ! the specified tolerance. Return false (test failed)
    pass = .false.
  else
    ! All elements of a and b are equal within the specified
    ! tolerance. Return true (test passed)
    pass = .true.
  endif

end function rel_assert_close

!===============================================================================
! Print test results to screen
!===============================================================================
subroutine print_result(test, pass, num_passed, num_tests)
  character(*) :: test
  logical :: pass
  integer :: num_passed, num_tests

  num_tests = num_tests + 1

  if (pass .eqv. .true.) then
    write(6,*) " * ", test, " =>"//achar(27)//"[32m passed"//achar(27)//"[0m"
    num_passed = num_passed + 1
  else
    write(6,*) " * ", test, " =>"//achar(27)//"[31m FAILED"//achar(27)//"[0m"
  endif

end subroutine print_result

!===============================================================================
! Print subtest results to screen
!===============================================================================
subroutine print_subresult(test, pass)
  character(*) :: test
  logical :: pass

  if (pass .eqv. .true.) then
    write(6,*) "    - ", test, " =>"//achar(27)//"[32m passed"//achar(27)//"[0m"
  else
    write(6,*) "    - ", test, " =>"//achar(27)//"[31m FAILED"//achar(27)//"[0m"
  endif

end subroutine print_subresult

!===============================================================================
! Set-up the ODE solver(s)
!===============================================================================
subroutine initiate_solver(pb)
  type(problem_type) :: pb

  ! Init RK45 requires:
  ! - pb%neqs, pb%mesh%nn
  ! - all pb%features set
  ! - pb%i_rns_law
  ! - pb%tau, pb%v, pb%theta, pb%sigma, pb%theta2
  ! - pb%time, pb%acc

  ! RK45 must be re-initialised each test

  if (SOLVER_TYPE == 2) then
    call init_rk45(pb)
  endif

  write(6,*) " * ODE solver (re)initialised"

end subroutine initiate_solver

end module unittests_aux
