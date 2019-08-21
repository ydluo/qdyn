! Collection of auxiliary functions

module ode_rk45_2

  use problem_class

  implicit none
  public

contains

!===============================================================================
! Helper routine to pack variables for solver
!===============================================================================
subroutine init_rk45_2(pb)

  type(problem_type), intent(inout) :: pb
  integer :: neqs

  neqs = pb%mesh%nn * pb%neqs

  pb%dt_did = pb%dt_try
  if (pb%dt_max == 0) then
    pb%dt_max = pb%tmax
  endif

  pb%rk45_2%t_coeffs = (/ 1.0/4.0, 3.0/8.0, 12.0/13.0, 1.0, 1.0/2.0 /)
  pb%rk45_2%coeffs = (/ &
            0.25, &
            3.0/32.0, 9.0/32.0, &
            1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0, &
            439.0/216.0, -8.0, 3680.0/513.0, -845.0/4104.0, &
            -8.0/27.0, 2.0, -3544.0/2565.0, 1859.0/4104.0, -11.0/40.0, &
            25.0/216.0, 1408.0/2565.0, 2197.0/4104.0, -1.0/5.0 &
            /)
  pb%rk45_2%error_coeffs = (/ &
            16.0/135.0 - 25/216.0, &
            6656.0/12825.0 - 1408/2565.0, &
            28561.0/56430.0 - 2197/4104.0, &
            -9.0/50.0 + 1.0/5.0, &
            2.0/55.0 &
            /)

  ! allocate k1, k2, dy1, dy2, etc.
  allocate( pb%rk45_2%k1(neqs), &
            pb%rk45_2%k2(neqs), &
            pb%rk45_2%k3(neqs), &
            pb%rk45_2%k4(neqs), &
            pb%rk45_2%k5(neqs), &
            pb%rk45_2%k6(neqs), &
            pb%rk45_2%dy1(neqs), &
            pb%rk45_2%dy2(neqs), &
            pb%rk45_2%dy3(neqs), &
            pb%rk45_2%dy4(neqs), &
            pb%rk45_2%dy5(neqs), &
            pb%rk45_2%dy6(neqs), &
            pb%rk45_2%e(neqs), &
            pb%rk45_2%y_new(neqs) &
  )

end subroutine init_rk45_2

!===============================================================================
! Helper routine to pack variables for solver
!===============================================================================
subroutine rkf45_d2(f, y, t, dtmax, rtol, atol, pb)

  external f
  type(problem_type), intent(inout) :: pb
  double precision, dimension(pb%neqs*pb%mesh%nn) :: y, y_new, e
  double precision :: t, dt, dtmax, rtol, atol, e_ratio, s, power
  integer :: i

  i = 0
  ! UPDATE dt after successful step
  dt = pb%dt_did
  y_new = pb%rk45_2%y_new
  e = pb%rk45_2%e

  do
    ! Increment loop counter
    i = i + 1

    ! Check if current dt > dtmax
    if (dt > dtmax) then
      ! Make step with size dtmax
      call rk45_step(f, y, t, dtmax, y_new, e, pb)
    else
      ! Make step with size dt
      call rk45_step(f, y, t, dt, y_new, e, pb)
    endif
    ! Calculate error
    e_ratio = maxval(e / (rtol * abs(y_new) + atol)) + 1e-12
    ! Check for NaNs
    if (isnan(e_ratio)) then
      ! Decimate step size and try again
      dt = dt * 0.1
      cycle
    endif

    ! Check for oscillations around e_ratio = 1
    ! Accept suboptimal answer when this is the case
    if ((i >= 100) .and. (e_ratio < 1)) then
      e_ratio = 1d0
    endif
    ! Check if current dt > dtmax and e_ratio < 1 (suboptimal)
    if ((dt > dtmax) .and. (e_ratio < 1)) then
      dt = dtmax
      e_ratio = 1d0
    endif

    ! If error is within a factor 10**0.1 = 25% from rtol, accept result
    if (abs(log10(e_ratio)) < 0.1) then
      y = y_new
      t = t + dt
      return
    else
      ! If error is too small or too big, select optimal exponent
      if (e_ratio > 1) then
        power = -0.2
      else
        power = -0.25
      endif
      ! Modifief for step size adjustment
      s = 0.999 * e_ratio**power
      ! Update step size
      dt = dt * s
    endif
  enddo

end subroutine rkf45_d2

!===============================================================================
! Perform a single integration step t -> t + dt
! Return y(t+dt) and error estimate
!===============================================================================
subroutine rk45_step(f, y, t, dt, y_new, e, pb)

  external f
  type(problem_type), intent(inout) :: pb
  double precision, dimension(pb%neqs*pb%mesh%nn) :: y, y_new, e, h
  double precision :: t, dt

  call f(t, y, pb%rk45_2%k1, pb)
  pb%rk45_2%dy1 =   pb%rk45_2%coeffs(1) * pb%rk45_2%k1

  call f(t + pb%rk45_2%t_coeffs(1) * dt, y + pb%rk45_2%dy1 * dt, pb%rk45_2%k2, pb)
  pb%rk45_2%dy2 =   pb%rk45_2%coeffs(2) * pb%rk45_2%k1 &
                  + pb%rk45_2%coeffs(3) * pb%rk45_2%k2

  call f(t + pb%rk45_2%t_coeffs(2) * dt, y + pb%rk45_2%dy2 * dt, pb%rk45_2%k3, pb)
  pb%rk45_2%dy3 =   pb%rk45_2%coeffs(4) * pb%rk45_2%k1 &
                  + pb%rk45_2%coeffs(5) * pb%rk45_2%k2 &
                  + pb%rk45_2%coeffs(6) * pb%rk45_2%k3

  call f(t + pb%rk45_2%t_coeffs(3) * dt, y + pb%rk45_2%dy3 * dt, pb%rk45_2%k4, pb)
  pb%rk45_2%dy4 =   pb%rk45_2%coeffs(7) * pb%rk45_2%k1 &
                  + pb%rk45_2%coeffs(8) * pb%rk45_2%k2 &
                  + pb%rk45_2%coeffs(9) * pb%rk45_2%k3 &
                  + pb%rk45_2%coeffs(10) * pb%rk45_2%k4

  call f(t + pb%rk45_2%t_coeffs(4) * dt, y + pb%rk45_2%dy4 * dt, pb%rk45_2%k5, pb)
  pb%rk45_2%dy5 =   pb%rk45_2%coeffs(11) * pb%rk45_2%k1 &
                  + pb%rk45_2%coeffs(12) * pb%rk45_2%k2 &
                  + pb%rk45_2%coeffs(13) * pb%rk45_2%k3 &
                  + pb%rk45_2%coeffs(14) * pb%rk45_2%k4 &
                  + pb%rk45_2%coeffs(15) * pb%rk45_2%k5

  call f(t + pb%rk45_2%t_coeffs(5) * dt, y + pb%rk45_2%dy5 * dt, pb%rk45_2%k6, pb)
  pb%rk45_2%dy6 =   pb%rk45_2%coeffs(16) * pb%rk45_2%k1 &
                  + pb%rk45_2%coeffs(17) * pb%rk45_2%k3 &
                  + pb%rk45_2%coeffs(18) * pb%rk45_2%k4 &
                  + pb%rk45_2%coeffs(19) * pb%rk45_2%k5

  y_new = y + pb%rk45_2%dy6 * dt
  h =   pb%rk45_2%error_coeffs(1) * pb%rk45_2%k1 &
      + pb%rk45_2%error_coeffs(2) * pb%rk45_2%k3 &
      + pb%rk45_2%error_coeffs(3) * pb%rk45_2%k4 &
      + pb%rk45_2%error_coeffs(4) * pb%rk45_2%k5 &
      + pb%rk45_2%error_coeffs(5) * pb%rk45_2%k6
  e = abs(h) * dt

end subroutine rk45_step


end module ode_rk45_2
