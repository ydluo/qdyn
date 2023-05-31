module friction

  ! This is the only module that needs modifications to implement a new friction law
  ! Follow the instructions in the comment blocks below that start with "! new friction law:"
  !
  ! Assumptions:
  !   The friction coefficient (mu) and the state variable rate (dtheta/dt)
  !   depend on slip velocity (v) and state variable (theta)
  !     mu = f(v,theta)
  !     dtheta/dt = g(v,theta)
  !   All friction properties can be spatially non-uniform

  use problem_class, only : problem_type
  use logger, only : log_msg

  implicit none
  private

  public  :: set_theta_star, compute_velocity_RSF, RSF_derivs, dtheta_dt

contains

!--------------------------------------------------------------------------------------
subroutine set_theta_star(pb)

  type(problem_type), intent(inout) :: pb

  select case (pb%i_rns_law)

  case (0)
    pb%theta_star = pb%dc/pb%v_star

  case (1)
    pb%theta_star = pb%dc/pb%v2

  case (2) ! 2018 SCEC Benchmark
    pb%theta_star = pb%dc/pb%v_star

  case (3) ! SEISMIC: the CNS friction law does not use theta_star
    pb%theta_star = 1

! new friction law:
!  case(xxx)
!    implement here your definition of theta_star (could be none)
!    pb%theta_star = ...

  case default
    stop 'set_theta_star: unknown friction law type'
  end select

end subroutine set_theta_star

!--------------------------------------------------------------------------------------
function compute_velocity_RSF(tau, sigma, theta, pb) result(v)

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: sigma, tau, theta
  double precision, dimension(pb%mesh%nn) :: mu, v

  mu = tau / sigma

  select case (pb%i_rns_law)

  case (0)
    v = pb%v_star * exp(pb%inv_a * (mu - pb%mu_star - pb%b * log(theta / pb%theta_star)))

  case (1)
    v = pb%v1 / (exp(- pb%inv_a * (mu - pb%mu_star - pb%b * log(theta / pb%theta_star + 1d0))) - 1)

  case (2) ! SCEC 2018 benchmark
    v = 2 * pb%v_star * sinh(mu * pb%inv_a) * exp(-pb%inv_a * (pb%mu_star + pb%b * log(theta / pb%theta_star)))

  ! Add viscous creep
  v = v + pb%inv_visc * tau

! new friction law:
!  case(xxx)
!    implement here your friction coefficient: mu = f(v,theta)
!    mu = ...

  case default
    stop 'compute_velocity_RSF: unknown friction law type'
  end select

end function compute_velocity_RSF

!--------------------------------------------------------------------------------------
subroutine dtheta_dt(v, theta, dth_dt, pb)

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: v, theta
  double precision, dimension(pb%mesh%nn) :: dth_dt, omega

  omega = v * theta / pb%dc
  select case (pb%itheta_law)

  case(0) ! "aging" in the no-healing approximation
    dth_dt = -omega

  case(1) ! "aging" law
    dth_dt = 1.d0-omega

  case(2) ! "slip" law
    dth_dt = -omega*log(omega)

! new friction law:
!  case(xxx)
!    implement here your state evolution law: dtheta/dt = g(v,theta)
!    dth_dt = ...

  case default
    stop 'dtheta_dt: unknown state evolution law type'
  end select

end subroutine dtheta_dt

!--------------------------------------------------------------------------------------
subroutine RSF_derivs(dV_dtau, dV_dtheta, dV_dsigma, v, theta, tau, sigma, pb)

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: v, theta, tau, sigma
  double precision, dimension(pb%mesh%nn) :: dV_dtau, dV_dtheta, dV_dsigma
  double precision, dimension(pb%mesh%nn) :: dummy

  select case (pb%i_rns_law)

  case(0)
    dummy = v * pb%inv_a
    dV_dtau = dummy / sigma
    dV_dtheta = - dummy * pb%b / theta

  case(1)
    dummy = v * (pb%v1 + v) / (pb%a * pb%v1)
    dV_dtau = dummy / sigma
    dV_dtheta = -dummy * pb%b / (theta + pb%theta_star)

  case(2) ! 2018 SCEC Benchmark
    dummy = 2 * pb%v_star * exp(-pb%inv_a * (pb%mu_star + pb%b * log(theta/pb%theta_star)))
    dV_dtau = sqrt(dummy**2 + v**2) / (pb%a * sigma)
    dV_dtheta = - v * pb%b / (pb%a * theta)

  case(3) ! SEISMIC: CNS model
    call log_msg("friction.f90::dmu_dv_dtheta is deprecated for the CNS model")
    stop

  case default
    call log_msg("dmu_dv_dtheta: unkown friction law type")
    stop
  end select

  ! Acceleration due to (effective) stress changes
  dV_dsigma = -dV_dtau * tau / sigma

  ! Viscous creep is state-independent, so it only contributes to dV_dtau
  dV_dtau = dV_dtau + pb%inv_visc

end subroutine RSF_derivs

end module friction
