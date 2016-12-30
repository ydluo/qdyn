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

! <SEISMIC>
! Rate-and-state friction law i_rns_law = 3 refers to Chen's microphysical model
! In this interpretation, rate-and-state friction results from the interplay
! between compaction by pressure solution creep and dilation by granular flow
! At low sliding velocities, pressure solution creep is relatively fast compared
! to granular dilatation, so that the gouge densifies (porosity reduces). By
! contrast, at high sliding velocities, the gouge dilates and 'softens', i.e.
! supports lower shear stress. At a given velocity, a steady-state porosity can
! be attained that corresponds to the porosity at which compaction balances
! dilatation, and hence results in a steady-state shear strength. Velocity-weakening
! behaviour emerges naturally from this model. To make the transition to Velocity-
! strengthening, a dislocation-glide type of model is adopted.
!
! For details and derivation of the model, see:
!   - Chen & Spiers (2016), JGR: Solid Earth, doi:10.1002/2016JB013470
!   - Chen & Spiers (subm.)
!
! In the view of rate-and-state friction, we take the porosity as the
! (microstructural) state, so that the constitutive relations can be implemented
! in the current QDyn framework by defining the partial- and time-derivatives of
! the friction coefficient and the porosity (state). The (spatially distributed)
! values of the material/gouge properties and initial porosity are read from the
! input script, and stored in pb%chen_params. Parameters defined specifically
! in this file are:
!   - tan_psi: dilatation angle
!   - mu_tilde: friction coefficient of the grain-boundary (rate-strengthening)
!   - e_ps_dot: compaction strain rate due to pressure solution creep
! </SEISMIC>

  use problem_class, only : problem_type

  implicit none
  private

  public  :: set_theta_star, friction_mu, dmu_dv_dtheta, dtheta_dt
  public  :: calc_tan_psi, calc_mu_tilde

contains

!--------------------------------------------------------------------------------------
subroutine set_theta_star(pb)

  type(problem_type), intent(inout) :: pb

  select case (pb%i_rns_law)

  case (0)
    pb%theta_star = pb%dc/pb%v_star

  case (1)
    pb%theta_star = pb%dc/pb%v2

! new friction law:
!  case(xxx)
!    implement here your definition of theta_star (could be none)
!    pb%theta_star = ...

  case default
    stop 'set_theta_star: unknown friction law type'
  end select

end subroutine set_theta_star

!--------------------------------------------------------------------------------------
function friction_mu(v,theta,pb) result(mu)

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: v, theta
  double precision, dimension(pb%mesh%nn) :: mu
  ! SEISMIC: additional parameters for Chen model
  double precision, dimension(pb%mesh%nn) :: mu_tilde, tan_psi


  select case (pb%i_rns_law)

  case (0)
    mu = pb%mu_star - pb%a*log(pb%v_star/v) + pb%b*log(theta/pb%theta_star)

  case (1)
    mu = pb%mu_star - pb%a*log(pb%v1/v+1d0) + pb%b*log(theta/pb%theta_star+1d0)

  case (3) ! SEISMIC: Chen's model
    mu_tilde = calc_mu_tilde(v, pb)
    tan_psi = calc_tan_psi(theta, pb)
    mu = (mu_tilde + tan_psi)/(1.0 - mu_tilde*tan_psi)

! new friction law:
!  case(xxx)
!    implement here your friction coefficient: mu = f(v,theta)
!    mu = ...

  case default
    stop 'friction_mu: unknown friction law type'
  end select

end function friction_mu

!--------------------------------------------------------------------------------------
function dtheta_dt(v,theta,pb) result(dth_dt)

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: v, theta
  double precision, dimension(pb%mesh%nn) :: dth_dt, omega
  double precision, dimension(pb%mesh%nn) :: tan_psi, e_ps_dot

  omega = v * theta / pb%dc
  select case (pb%itheta_law)

  case(0) ! "aging" in the no-healing approximation
    dth_dt = -omega

  case(1) ! "aging" law
    dth_dt = 1.d0-omega

  case(2) ! "slip" law
    dth_dt = -omega*log(omega)

  case (3) ! SEISMIC: Chen's model
    tan_psi = calc_tan_psi(theta, pb)
    e_ps_dot = pb%chen_params%IPS_const*pb%sigma/(2*pb%chen_params%phi0 - 2*theta)**2
    dth_dt = (1.0 - theta)*(tan_psi*v/pb%chen_params%w - e_ps_dot)

! new friction law:
!  case(xxx)
!    implement here your state evolution law: dtheta/dt = g(v,theta)
!    dth_dt = ...

  case default
    stop 'dtheta_dt: unknown state evolution law type'
  end select

end function dtheta_dt

!--------------------------------------------------------------------------------------
subroutine dmu_dv_dtheta(dmu_dv,dmu_dtheta,v,dv_dt,theta,pb)

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: v, dv_dt, theta
  double precision, dimension(pb%mesh%nn), intent(out) :: dmu_dv, dmu_dtheta
  double precision, dimension(pb%mesh%nn) :: tan_psi, mu_tilde, dth_dt, denom

  select case (pb%i_rns_law)

  case(0)
    dmu_dtheta = pb%b / theta
    dmu_dv = pb%a / v

  case(1)
    dmu_dtheta = pb%b * pb%v2 / ( pb%v2*theta + pb%dc )
    dmu_dv = pb%a * pb%v1 / v / ( pb%v1 + v )

  case(3) ! SEISMIC: Chen's model
    tan_psi = calc_tan_psi(theta, pb)
    mu_tilde = calc_mu_tilde(v, pb)
    dth_dt = dtheta_dt(v, theta, pb)
    ! SEISMIC: pre-define the denominator
    denom = 1.0/(1.0 - mu_tilde*tan_psi)**2
    dmu_dtheta = 2.0*pb%chen_params%H*( mu_tilde*(2.0*tan_psi + mu_tilde) - 1.0 )*denom
    dmu_dv = &
      (1.0 - mu_tilde*tan_psi)*(pb%chen_params%a/v - 2.0*pb%chen_params%H*dth_dt/dv_dt) &
      - (mu_tilde + tan_psi)*(2.0*pb%chen_params%H*mu_tilde*dth_dt/dv_dt - tan_psi*pb%chen_params%a/v) &
      * denom

! new friction law:
!  case(xxx)
!    implement here the partial derivatives of the friction coefficient
!    dmu_dtheta = ...
!    dmu_dv = ...

  case default
    stop 'dmu_dv_dtheta: unknown friction law type'
  end select

end subroutine dmu_dv_dtheta

!--------------------------------------------------------------------------------------
! SEISMIC: calculate the dilatation angle, which is a geometric consequence
! of the instantaneous porosity (theta)
!--------------------------------------------------------------------------------------
function calc_tan_psi(theta,pb) result(tan_psi)

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: theta
  double precision, dimension(pb%mesh%nn) :: tan_psi

  tan_psi = 2*pb%chen_params%H*(pb%chen_params%phi0 - theta)

end function calc_tan_psi

!--------------------------------------------------------------------------------------
! SEISMIC: calculate the grain-boundary coefficient of friction, which is
! logarithmically rate-strengthening, assuming deformation of the contact
! asperities by dislocation-glide
!--------------------------------------------------------------------------------------
function calc_mu_tilde(v,pb) result(mu_tilde)

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: v
  double precision, dimension(pb%mesh%nn) :: mu_tilde

  ! SEISMIC: slowness_star is defined as 1/V_star in the input script
  mu_tilde = pb%chen_params%mu_tilde_star + pb%chen_params%a * log(v*pb%chen_params%slowness_star)

end function calc_mu_tilde

end module friction
