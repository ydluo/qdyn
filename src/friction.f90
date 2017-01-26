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
! Rate-and-state friction law i_rns_law = 3 refers to Chen's microphysical model,
! with itheta_law = 3 to diffusion controlled pressure solution, and itheta_law = 4
! to dissolution controlled pressure solution creep.
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
!   - Niemeijer & Spiers (2007), J. Geophys. Res., doi:10.1029/2007JB00500
!   - Chen & Spiers (2016), JGR: Solid Earth, doi:10.1002/2016JB013470
!
! In the view of rate-and-state friction, we take the porosity as the
! (microstructural) state, and the constitutive relations are solved for the
! shear stress (rather than velocity; see derivs_all.f90). The (spatially non-uniform)
! values of the material/gouge properties and initial porosity are read from the
! input script, and stored in pb%chen_params.
!
! </SEISMIC>

  use problem_class, only : problem_type

  implicit none
  private

  public  :: set_theta_star, friction_mu, dmu_dv_dtheta, dtheta_dt
  public  :: compute_velocity, dalpha_dt
  private  :: calc_tan_psi, calc_mu_tilde, calc_e_ps

contains

!--------------------------------------------------------------------------------------
subroutine set_theta_star(pb)

  type(problem_type), intent(inout) :: pb

  select case (pb%i_rns_law)

  case (0)
    pb%theta_star = pb%dc/pb%v_star

  case (1)
    pb%theta_star = pb%dc/pb%v2

  case (3) ! SEISMIC: Chen's friction law does not use theta_star
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
function friction_mu(v,theta,pb) result(mu)

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: v, theta
  double precision, dimension(pb%mesh%nn) :: mu
  ! SEISMIC: additional parameters for Chen model
  double precision, dimension(pb%mesh%nn) :: mu_tilde, tan_psi, mu_gr, mu_ps


  select case (pb%i_rns_law)

  case (0)
    mu = pb%mu_star - pb%a*log(pb%v_star/v) + pb%b*log(theta/pb%theta_star)

  case (1)
    mu = pb%mu_star - pb%a*log(pb%v1/v+1d0) + pb%b*log(theta/pb%theta_star+1d0)

  case (3) ! SEISMIC: Chen's model

    write (6,*) "friction.f90::friciton_mu is deprecated for Chen's model"

    ! SEISMIC: this function should no longer be called to determine the
    ! state of shear stress, but is kept here for reference. The state of stress
    ! is solved for from a kinematic balance equation. This function calculates
    ! the steady-state friction corresponding to the slip velocity and the porosity
    ! as read from the input file.
    ! TODO: replace the initiation through this function entirely by the mesh
    ! input (i.e. set state of stress from input file)

    tan_psi = calc_tan_psi(theta, pb)   ! Dilatation angle
    ! SEISMIC: calculate the grain-boundary friction, assuming all deformation
    ! is accommodated by granular flow (i.e. v_total = v_gr)
    mu_tilde = calc_mu_tilde(v/pb%chen_params%w, pb)

    ! SEISMIC: Calculate the frictional strength controlled by granular flow
    mu_gr = (mu_tilde + tan_psi)/(1 - mu_tilde*tan_psi)

    ! SEISMIC: Next, calculate the shear strength (divided by normal stress)
    ! when the strength is controlled by viscous flow (no granular flow)
    select case (pb%itheta_law)

    case (3) ! SEISMIC: Chen's model, diffusion controlled
      mu_ps = v*((2*pb%chen_params%phi0 - 2*theta)**2)/(pb%chen_params%w*pb%chen_params%IPS_const_diff*pb%sigma)
    case (4) ! SEISMIC: Chen's model, dissolution controlled
      mu_ps = (pb%chen_params%phi0 - theta)/(pb%sigma*pb%chen_params%phi0*pb%chen_params%IPS_const_diss2)* &
              log(abs(v)/(pb%chen_params%w*pb%chen_params%IPS_const_diss1) + 1)
    case default
      write(6,*) "friction_mu: Chen's friction model is selected (i_rns_law == 3),"
      write(6,*) "but itheta_law is unsupported (must be either 3 or 4)"
      write(6,*) "3 = diffusion-, 4 = dissolution controlled pressure solution creep"
      stop

    end select

    ! SEISMIC: this weird looking function approximates the min(mu_gr, mu_ps) function,
    ! but is continuously differentiable, so that the solver doesn't complain.
    ! Higher order powers yield a closer approximation
    ! By approximation, the mechanism that offers the lowest shear strength controls
    ! the overall shear strength of the gouge
    mu = 1/(1/mu_gr**5 + 1/mu_ps**5)**(1.0/5.0)

! new friction law:
!  case(xxx)
!    implement here your friction coefficient: mu = f(v,theta)
!    mu = ...

  case default
    stop 'friction_mu: unknown friction law type'
  end select

end function friction_mu

!--------------------------------------------------------------------------------------
function dtheta_dt(v,tau,sigma,theta,pb) result(dth_dt)

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: v, tau, sigma, theta
  double precision, dimension(pb%mesh%nn) :: dth_dt, omega
  double precision, dimension(pb%mesh%nn) :: tan_psi, e_ps_dot
  double precision, dimension(pb%mesh%nn) :: y_ps, y_gr

  omega = v * theta / pb%dc
  select case (pb%itheta_law)

  case(0) ! "aging" in the no-healing approximation
    dth_dt = -omega

  case(1) ! "aging" law
    dth_dt = 1.d0-omega

  case(2) ! "slip" law
    dth_dt = -omega*log(omega)

  case (3, 4) ! SEISMIC: Chen's model
    tan_psi = calc_tan_psi(theta, pb)   ! Dilatation angle
    e_ps_dot = calc_e_ps(sigma, theta, .true., pb)   ! Compaction rate
    y_ps = calc_e_ps(tau, theta, .false., pb)   ! Press. soln. shear rate
    y_gr = v/pb%chen_params%w - y_ps
    ! SEISMIC: the nett rate of change of porosity is given by the relation
    ! phi_dot = -(1 - phi)*strain_rate, where the strain rate equals
    ! the compaction rate by pressure solution (e_ps_dot) minus the dilatation
    ! rate by granular flow (tan_psi*v/w). Note that compaction is measured
    ! positive, dilatation negative
    dth_dt = (1 - theta)*(tan_psi*y_gr - e_ps_dot)


! new friction law:
!  case(xxx)
!    implement here your state evolution law: dtheta/dt = g(v,theta)
!    dth_dt = ...

  case default
    stop 'dtheta_dt: unknown state evolution law type'
  end select

end function dtheta_dt


function dalpha_dt(v,tau,sigma,theta,alpha,pb) result(da_dt)

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: tau, sigma, theta, alpha, v
  double precision, dimension(pb%mesh%nn) :: da_dt, psi, tan_psi, sin_psi, cos_psi
  double precision, dimension(pb%mesh%nn) :: sigma_i, q, power_term
  double precision, dimension(pb%mesh%nn) :: y_ps, y_gr

  q = 2*(pb%chen_params%phi0 - theta)
  tan_psi = pb%chen_params%H*q
  psi = atan(tan_psi)
  cos_psi = cos(psi)
  sin_psi = tan_psi * cos_psi

  y_ps = calc_e_ps(tau, theta, .false., pb)   ! Press. soln. shear rate
  y_gr = v/pb%chen_params%w - y_ps

  ! 1.9 is approximately 6/pi, where 6 is the average grain coordination number
  sigma_i = 1.9*(sigma*cos_psi + tau*sin_psi)/(alpha*q)
  power_term = ((pb%coh_params%alpha_c - alpha)/(pb%coh_params%alpha_c - pb%coh_params%alpha0))**1.3
  da_dt = (pb%coh_params%NG_const*power_term/q)*(pb%coh_params%E_surf - 0.5*pb%coh_params%compl*sigma_i**2) &
          - y_gr*(alpha - pb%coh_params%alpha0)

end function dalpha_dt

!--------------------------------------------------------------------------------------
subroutine dmu_dv_dtheta(dmu_dv,dmu_dtheta,v,theta,pb)

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: v, theta
  double precision, dimension(pb%mesh%nn), intent(out) :: dmu_dv, dmu_dtheta

  select case (pb%i_rns_law)

  case(0)
    dmu_dtheta = pb%b / theta
    dmu_dv = pb%a / v

  case(1)
    dmu_dtheta = pb%b * pb%v2 / ( pb%v2*theta + pb%dc )
    dmu_dv = pb%a * pb%v1 / v / ( pb%v1 + v )

  case default
    write (6,*) "dmu_dv_dtheta: unkown friction law type"
    write (6,*) "This function should only be called for i_rns_law = 0 and 1"
    write (6,*) "Chen's model (i_rns_law = 3) does not compute these partial derivatives"
    stop
  end select

end subroutine dmu_dv_dtheta

!--------------------------------------------------------------------------------------
! SEISMIC: calculate the dilatation angle, which is a geometric consequence
! of the instantaneous porosity (theta)
!--------------------------------------------------------------------------------------
function compute_velocity(tau,sigma,theta,alpha,pb) result(v)

  use constants, only : PI

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: tau, sigma, theta, alpha
  double precision, dimension(pb%mesh%nn) :: v

  ! SEISMIC: define some extra parameters for Chen's friction law
  double precision, dimension(pb%mesh%nn) :: tan_psi, mu_tilde, denom, y_ps
  double precision, dimension(pb%mesh%nn) :: mu_star, y_gr, tau_gr
  double precision, dimension(pb%mesh%nn) :: cos_psi, cohesion

  tan_psi = calc_tan_psi(theta, pb)         ! Dilatation angle

  ! Pre-compute the denominator for efficiency
  denom = 1.0/(pb%chen_params%a*(sigma+tau*tan_psi))
  mu_star = pb%chen_params%mu_tilde_star

  ! The granular flow strain rate for the case when the shear stress approaches
  ! the shear strength of the grain contacts
  y_gr = pb%chen_params%y_gr_star*exp((tau*(1 - mu_star*tan_psi) - sigma*(mu_star + tan_psi))*denom)

  ! Shear strain rate [1/s] due to viscous flow (pressure solution)
  y_ps = calc_e_ps(tau, theta, .false., pb)
  mu_tilde = calc_mu_tilde(y_gr, pb)           ! Grain-boundary friction

  cohesion = 0
  if (pb%features%cohesion == 1) then
    cos_psi = cos(atan(tan_psi))
    cohesion = (PI/pb%chen_params%H)*(tan_psi*alpha*pb%coh_params%C_star)/(cos_psi*(1-mu_tilde*tan_psi))
  endif

  ! Calculate the shear strength of the grain contacts, if the gouge were
  ! to deform by granular flow. To allow for granular flow (y_gr != 0),
  ! tau_gr should be near tau (frictional 'yield')
  tau_gr = (mu_tilde + tan_psi)/(1 - mu_tilde*tan_psi)*sigma + cohesion
  ! Use the error function to approximate the frictional yielding of the
  ! contacts. If tau_gr < 0.9*tau, this function quickly approaches zero
  ! NOTE: the use of an error function is more expensive than an if-statement
  ! but is continuously differentiable, which is required for the solver
  y_gr = 0.5*(1 + erf(100*(tau-0.9*tau_gr)/tau_gr))*y_gr

  ! The total slip velocity is the combined contribution of granular flow
  ! and pressure solution (parallel processes)
  v = pb%chen_params%w*(y_gr + y_ps)


end function compute_velocity

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
! logarithmically rate-strengthening, assuming some atomic scale jump process
!--------------------------------------------------------------------------------------
function calc_mu_tilde(y_gr,pb) result(mu_tilde)

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: y_gr
  double precision, dimension(pb%mesh%nn) :: mu_tilde

  ! SEISMIC NOTE: the granular flow strain rate should be used here,
  ! not the total strain rate
  mu_tilde = pb%chen_params%mu_tilde_star + pb%chen_params%a * log(abs(y_gr)/pb%chen_params%y_gr_star)

end function calc_mu_tilde

!--------------------------------------------------------------------------------------
! SEISMIC: calculate the rate of pressure solution as a function of effective
! normal or shear stress, and porosity. To prevent porosities from going negative
! compaction rates can be truncated by an error function. The logical 'truncate'
! indicates if truncation is needed. This only applies to compaction in the
! normal direction, and not to shear deformation (viscous flow), hence truncate
! should be false when calculating shear strain rate
!--------------------------------------------------------------------------------------
function calc_e_ps(sigma,theta,truncate,pb) result(e_ps_dot)

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: sigma, theta
  double precision, dimension(pb%mesh%nn) :: e_ps_dot
  logical, intent(in) :: truncate

  select case (pb%itheta_law)

  case (3) ! Chen's model, diffusion controlled
    e_ps_dot = pb%chen_params%IPS_const_diff*sigma/(2*pb%chen_params%phi0 - 2*theta)**2
  case (4) ! Chen's model, dissolution controlled
    e_ps_dot =  pb%chen_params%IPS_const_diss1*(exp(pb%chen_params%IPS_const_diss2*sigma* &
                pb%chen_params%phi0/(pb%chen_params%phi0 - theta)) - 1)
  case default
    stop 'calc_e_ps: unsupported pressure solution law (3 = diffusion, 4 = dissolution controlled)'
  end select

  ! If truncation is required, use an error function to set the strain rate
  ! to zero. A cut-off porosity of about 3% is chosen here, which corresponds
  ! to the percolation limit of aggregates
  if (truncate .eqv. .true.) e_ps_dot = 0.5*(1 + erf(100*(theta-0.03)))*e_ps_dot
  !!if (truncate .eqv. .true.) e_ps_dot = erf(2*theta/0.05)*e_ps_dot

end function calc_e_ps

end module friction
