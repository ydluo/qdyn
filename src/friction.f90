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
!
! Rate-and-state friction law i_rns_law = 3 refers to the CNS (Chen-Niemeijer-Spiers)
! microphysical model, with itheta_law = 3 to diffusion controlled pressure solution,
! and itheta_law = 4 to dissolution controlled pressure solution creep.
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
! input script, and stored in pb%cns_params.
!
! </SEISMIC>

  use problem_class, only : problem_type

  implicit none
  private

  public  :: set_theta_star, friction_mu, dmu_dv_dtheta, dtheta_dt
  public  :: compute_velocity, dalpha_dt, CNS_derivs
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
function friction_mu(v,theta,pb) result(mu)

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: v, theta
  double precision, dimension(pb%mesh%nn) :: mu

  select case (pb%i_rns_law)

  case (0)
    mu = pb%mu_star - pb%a*log(pb%v_star/v) + pb%b*log(theta/pb%theta_star)

  case (1)
    mu = pb%mu_star - pb%a*log(pb%v1/v+1d0) + pb%b*log(theta/pb%theta_star+1d0)

  case (3) ! SEISMIC: CNS model
    write (6,*) "friction.f90::friciton_mu is deprecated for the CNS model"
    stop

! new friction law:
!  case(xxx)
!    implement here your friction coefficient: mu = f(v,theta)
!    mu = ...

  case default
    stop 'friction_mu: unknown friction law type'
  end select

end function friction_mu

!--------------------------------------------------------------------------------------
subroutine dtheta_dt(v,tau,sigma,theta,theta2,dth_dt,dth2_dt,pb)

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: v, tau, sigma
  double precision, dimension(pb%mesh%nn), intent(in) :: theta, theta2
  double precision, dimension(pb%mesh%nn) :: dth_dt, dth2_dt, omega
  double precision, dimension(pb%mesh%nn) :: tan_psi, e_ps_dot, e_ps_dot_bulk
  double precision, dimension(pb%mesh%nn) :: y_ps, y_gr, y_ps_bulk

  omega = v * theta / pb%dc
  select case (pb%itheta_law)

  case(0) ! "aging" in the no-healing approximation
    dth_dt = -omega

  case(1) ! "aging" law
    dth_dt = 1.d0-omega

  case(2) ! "slip" law
    dth_dt = -omega*log(omega)

  case (3) ! SEISMIC: CNS model
    tan_psi = calc_tan_psi(theta, pb)   ! Dilatation angle
    e_ps_dot = calc_e_ps(sigma, theta, .true., .false., pb)   ! Compaction rate
    y_ps = calc_e_ps(tau, theta, .false., .false., pb)   ! Press. soln. shear rate

    e_ps_dot_bulk = 0
    y_ps_bulk = 0

    if (pb%features%localisation == 1) then
      e_ps_dot_bulk = calc_e_ps(sigma, theta2, .true., .true., pb)
      y_ps_bulk = calc_e_ps(tau, theta2, .false., .true., pb)
    endif

    y_gr = (1.0/pb%cns_params%lambda)*(v/pb%cns_params%w - (1-pb%cns_params%lambda)*y_ps_bulk) - y_ps
    ! SEISMIC: the nett rate of change of porosity is given by the relation
    ! phi_dot = -(1 - phi)*strain_rate, where the strain rate equals
    ! the compaction rate by pressure solution (e_ps_dot) minus the dilatation
    ! rate by granular flow (tan_psi*y_gr). Note that compaction is measured
    ! positive, dilatation negative
    dth_dt = (1 - theta)*(tan_psi*y_gr - e_ps_dot)
    dth2_dt = -(1 - theta2)*e_ps_dot_bulk


! new friction law:
!  case(xxx)
!    implement here your state evolution law: dtheta/dt = g(v,theta)
!    dth_dt = ...

  case default
    stop 'dtheta_dt: unknown state evolution law type'
  end select

end subroutine dtheta_dt

!--------------------------------------------------------------------------------------
function dalpha_dt(v,tau,sigma,theta,alpha,pb) result(da_dt)

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: tau, sigma, theta, alpha, v
  double precision, dimension(pb%mesh%nn) :: da_dt, psi, tan_psi, sin_psi, cos_psi
  double precision, dimension(pb%mesh%nn) :: sigma_i, q, power_term, delta_alpha
  double precision, dimension(pb%mesh%nn) :: y_ps, y_gr, x, y
  double precision, parameter :: small = 1e-15

  write(6,*) "ERROR: The grain-boundary healing feature is work in progress..."
  stop

  q = 2*(pb%cns_params%phi0 - theta)
  tan_psi = pb%cns_params%H*q
  psi = atan(tan_psi)
  cos_psi = cos(psi)
  sin_psi = tan_psi * cos_psi

  y_ps = calc_e_ps(tau, theta, .false., .false., pb)   ! Press. soln. shear rate
  y_gr = v/pb%cns_params%w - y_ps

  delta_alpha = alpha - pb%coh_params%alpha0

  ! 1.9 is approximately 6/pi, where 6 is the average grain coordination number (doesn't have to be precise)
  sigma_i = 1.9*(sigma*cos_psi + tau*sin_psi)/(alpha*q)
  power_term = ((pb%coh_params%alpha_c - alpha)/(pb%coh_params%alpha_c - pb%coh_params%alpha0))**1.3
  da_dt = 0*(pb%coh_params%NG_const*power_term/q)*(pb%coh_params%E_surf - 0.5*pb%coh_params%compl*sigma_i**2) &
          - y_gr

  ! Define logical expressions, which should yield either 0 or 1 depending on
  ! the sign of da_dt and delta_alpha. The addition of small should prevent
  ! zero division
  x = 0.5*(1.0 + (da_dt + small)/(abs(da_dt) + small))
  y = 0.5*(1.0 + (delta_alpha + small)/(abs(delta_alpha) + small))

  ! If da_dt < 0 and delta_alpha < 0, set da_dt to zero. Leave as is otherwise
  da_dt = (x + (1-x)*y)*da_dt

end function dalpha_dt

!--------------------------------------------------------------------------------------
subroutine dmu_dv_dtheta(dmu_dv,dmu_dtheta,v,tau,sigma,theta,theta2,pb)

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: v, tau, sigma, theta, theta2
  double precision, dimension(pb%mesh%nn), intent(out) :: dmu_dv, dmu_dtheta

  ! SEISMIC: define some extra parameters for Chen's friction law
  double precision, dimension(pb%mesh%nn) :: tan_psi, mu_tilde, dth_dt, denom, y_ps, delta
  double precision, dimension(pb%mesh%nn) :: mu_star, dummy_var, f_phi_sqrt, V_gr, y_gr
  double precision, dimension(pb%mesh%nn) :: dv_dtau_ps_d, dv_dtheta_ps_d, dv_dtau_ps_s, dv_dtheta_ps_s
  double precision, dimension(pb%mesh%nn) :: dv_dtau, dv_dtheta, L_sb, y_ps_bulk

  select case (pb%i_rns_law)

  case(0)
    dmu_dtheta = pb%b / theta
    dmu_dv = pb%a / v

  case(1)
    dmu_dtheta = pb%b * pb%v2 / ( pb%v2*theta + pb%dc )
    dmu_dv = pb%a * pb%v1 / v / ( pb%v1 + v )

  case(3) ! SEISMIC: CNS model

    ! Thickness of the localised zone
    L_sb = pb%cns_params%lambda*pb%cns_params%w

    ! Shear strain rate [1/s] due to viscous flow (pressure solution)
    y_ps = calc_e_ps(tau, theta, .false., .false., pb)
    y_ps_bulk = 0
    if (pb%features%localisation == 1) then
      y_ps_bulk = calc_e_ps(tau, theta2, .false., .true., pb)
    endif

    ! Granular flow strain rate [1/s] is the total strain rate minus the pressure
    ! solution rate (taking localisation into account)
    y_gr = (1.0/pb%cns_params%lambda)*(v/pb%cns_params%w - (1-pb%cns_params%lambda)*y_ps_bulk) - y_ps

    tan_psi = calc_tan_psi(theta, pb)         ! Dilatation angle
    mu_tilde = calc_mu_tilde(y_gr, pb)        ! Grain-boundary friction
    mu_star = pb%cns_params%mu_tilde_star     ! Reference GB friction

    ! Pre-compute the denominator for efficiency
    denom = 1.0/(pb%cns_params%a*(sigma+tau*tan_psi))
    ! Pre-compute a dummy variable that we will re-use a few times
    dummy_var = (tau*(1 - mu_star*tan_psi) - sigma*(mu_star + tan_psi))*denom

    ! Pre-compute the square root of the porosity function. This will serve as
    ! a basis for calculating the full porosity functions for either itheta_law
    f_phi_sqrt = 1.0/(2*(pb%cns_params%phi0 - theta))

    ! Granular flow velocity (not strain rate) [units m/s]
    V_gr = L_sb*pb%cns_params%y_gr_star*exp(dummy_var)

    ! Next, calculate the partial derivatives dv_dtau and dv_dtheta for the
    ! pressure solution creep. The functional form depends on the rate-limiting
    ! mechanism (diffusion or dissolution). Both mechanisms are calculated, but
    ! if diffusion is rate controlling, then dissolution kinetics are set to zero
    ! and vice-versa, so addition of the two mechanisms results in only one rate-
    ! controlling mechanism for each fault element.
    ! Note that the pressure solution strain rate is multiplied by the
    ! the gouge layer thickness (pb%cns_params%w) to obtain a velocity

    ! Diffusion controlled pressure solution
    dv_dtau_ps_d = L_sb*pb%cns_params%IPS_const_diff*f_phi_sqrt**2
    dv_dtheta_ps_d = 4*dv_dtau_ps_d*f_phi_sqrt*tau

    ! Dissolution controlled pressure solution
    dv_dtau_ps_s =  L_sb*(y_ps + pb%cns_params%IPS_const_diss1)* &
                    pb%cns_params%IPS_const_diss2*2*pb%cns_params%phi0*f_phi_sqrt
    dv_dtheta_ps_s = 2*dv_dtau_ps_s*f_phi_sqrt*tau

    ! The partial derivatives dv_dtau and dv_dtheta of the overall slip velocity
    ! are the sum of the partial derivatives of the pressure solution and
    ! granular flow velocities.

    dv_dtau = dv_dtau_ps_d + dv_dtau_ps_s + &
              V_gr*((1-mu_star*tan_psi)*denom - tan_psi*dummy_var*pb%cns_params%a*denom)
    dv_dtheta = dv_dtau_ps_d + dv_dtau_ps_s + &
                V_gr*(2*pb%cns_params%H*(sigma + mu_star*tau)*denom + &
                2*pb%cns_params%H*tau*dummy_var*pb%cns_params%a*denom)

    ! For the CNS model, this function should return dV/dtau and dV/dtheta instead
    ! of dmu/dV and dmu/dtheta. For compatibility with classical rate-and-state,
    ! the same variable names (dmu_dv and dmu_dtheta) are used to store these quantities
    dmu_dv = dv_dtau
    dmu_dtheta = dv_dtheta

  case default
    write (6,*) "dmu_dv_dtheta: unkown friction law type"
    stop
  end select

end subroutine dmu_dv_dtheta

!--------------------------------------------------------------------------------------
! SEISMIC: the subroutine below was written to optimise/reduce function calls, and to
! improve code handling (less code duplication, hopefully fewer bugs). This subroutine
! should replace calls to dtheta_dt, compute_velocity, and dmu_dv_dtheta
! The subroutine only handles the CNS model
!--------------------------------------------------------------------------------------
subroutine CNS_derivs(v, dth_dt, dth2_dt, dv_dtau, dv_dtheta, tau, sigma, theta, theta2, alpha, pb)

  use constants, only : PI

  ! Input variables
  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: tau, sigma, theta, theta2, alpha

  ! Output variables
  double precision, dimension(pb%mesh%nn) :: v, dth_dt, dth2_dt, dv_dtau, dv_dtheta

  ! Internal variables
  double precision, dimension(pb%mesh%nn) :: tan_psi, mu_tilde, denom, dummy_var
  double precision, dimension(pb%mesh%nn) :: mu_star, y_gr, tau_gr
  double precision, dimension(pb%mesh%nn) :: cos_psi, cohesion, f_phi_sqrt, L_sb
  double precision, dimension(pb%mesh%nn) :: e_ps, y_ps, e_ps_bulk, y_ps_bulk
  double precision, dimension(pb%mesh%nn) :: dv_dtau_ps_d, dv_dtau_ps_s
  double precision, dimension(pb%mesh%nn) :: dv_dtheta_ps_d, dv_dtheta_ps_s

  if (pb%i_rns_law /= 3) then
    write(6,*) "The friction.f90::CNS_derivs subroutine only applies to the CNS model"
    write(6,*) "This routine should only be called when i_rns_law == 3"
    stop
  endif

  ! Pre-define some internal variables
  tan_psi = calc_tan_psi(theta, pb)             ! Dilatation angle
  mu_star = pb%cns_params%mu_tilde_star         ! Reference GB friction
  L_sb = pb%cns_params%lambda*pb%cns_params%w   ! Thickness of localised zone

  ! Pre-compute the denominator for efficiency
  denom = 1.0/(pb%cns_params%a*(sigma+tau*tan_psi))
  ! Pre-compute a dummy variable that we will re-use a few times
  dummy_var = (tau*(1 - mu_star*tan_psi) - sigma*(mu_star + tan_psi))*denom

  ! The granular flow strain rate for the case when the shear stress approaches
  ! the shear strength of the grain contacts
  y_gr = pb%cns_params%y_gr_star*exp(dummy_var)

  mu_tilde = calc_mu_tilde(y_gr, pb)           ! Grain-boundary friction

  ! Pre-compute the square root of the porosity function. This will serve as
  ! a basis for calculating the porosity functions
  f_phi_sqrt = 1.0/(2*(pb%cns_params%phi0 - theta))

  ! Shear strain rate [1/s] due to viscous flow (pressure solution)
  e_ps = calc_e_ps(sigma, theta, .true., .false., pb)   ! Compaction rate
  y_ps = calc_e_ps(tau, theta, .false., .false., pb)    ! Shear creep rate

  ! If localisation is requested, calculate pressure solution rates
  ! in the bulk zone, otherwise set to zero
  e_ps_bulk = 0
  y_ps_bulk = 0
  if (pb%features%localisation == 1) then
    e_ps_bulk = calc_e_ps(sigma, theta2, .true., .true., pb)
    y_ps_bulk = calc_e_ps(tau, theta2, .false., .true., pb)
  endif

  ! Cohesion is still in progress...
  cohesion = 0
  if (pb%features%cohesion == 1) then
    cos_psi = cos(atan(tan_psi))
    cohesion = (PI/pb%cns_params%H)*(tan_psi*alpha*pb%coh_params%C_star)/(cos_psi*(1-mu_tilde*tan_psi))
  endif

  ! The total slip velocity is the combined contribution of granular flow
  ! and pressure solution (parallel processes)
  v = pb%cns_params%w*(pb%cns_params%lambda*(y_gr + y_ps) + (1-pb%cns_params%lambda)*y_ps_bulk)

  ! SEISMIC: the nett rate of change of porosity is given by the relation
  ! phi_dot = -(1 - phi)*strain_rate, where the strain rate equals
  ! the compaction rate by pressure solution (e_ps) minus the dilatation
  ! rate by granular flow (tan_psi*v/w). Note that compaction is measured
  ! positive, dilatation negative
  dth_dt = (1 - theta)*(tan_psi*y_gr - e_ps)
  dth2_dt = -(1 - theta2)*e_ps_bulk

  ! Next, calculate the partial derivatives dv_dtau and dv_dtheta for the
  ! pressure solution creep. The functional form depends on the rate-limiting
  ! mechanism (diffusion, dissolution, or precipitation)
  ! Note that the pressure solution strain rate is multiplied by the
  ! the gouge layer thickness (pb%cns_params%w) to obtain a velocity

  dv_dtau_ps_d = L_sb*pb%cns_params%IPS_const_diff*f_phi_sqrt**2
  dv_dtheta_ps_d = 4*dv_dtau_ps_d*f_phi_sqrt*tau

  dv_dtau_ps_s =  L_sb*(y_ps + pb%cns_params%IPS_const_diss1)* &
                  pb%cns_params%IPS_const_diss2*2*pb%cns_params%phi0*f_phi_sqrt
  dv_dtheta_ps_s = 2*dv_dtau_ps_s*f_phi_sqrt*tau

  ! The partial derivatives dv_dtau and dv_dtheta of the overall slip velocity
  ! are the sum of the partial derivatives of the pressure solution and
  ! granular flow velocities.

  dv_dtau = dv_dtau_ps_d + dv_dtau_ps_s + &
            L_sb*y_gr*((1-mu_star*tan_psi)*denom - tan_psi*dummy_var*pb%cns_params%a*denom)
  dv_dtheta = dv_dtau_ps_d + dv_dtau_ps_s + &
              L_sb*y_gr*(2*pb%cns_params%H*(sigma + mu_star*tau)*denom + &
              2*pb%cns_params%H*tau*dummy_var*pb%cns_params%a*denom)

end subroutine CNS_derivs

!--------------------------------------------------------------------------------------
! SEISMIC: calculate the dilatation angle, which is a geometric consequence
! of the instantaneous porosity (theta)
!--------------------------------------------------------------------------------------
function compute_velocity(tau,sigma,theta,theta2,alpha,pb) result(v)

  use constants, only : PI

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: tau, sigma, theta
  double precision, dimension(pb%mesh%nn), intent(in) :: theta2, alpha
  double precision, dimension(pb%mesh%nn) :: v

  ! SEISMIC: define some extra parameters for the CNS friction law
  double precision, dimension(pb%mesh%nn) :: tan_psi, mu_tilde, denom, y_ps
  double precision, dimension(pb%mesh%nn) :: mu_star, y_gr, tau_gr
  double precision, dimension(pb%mesh%nn) :: cos_psi, cohesion
  double precision, dimension(pb%mesh%nn) :: y_ps_bulk

  tan_psi = calc_tan_psi(theta, pb)         ! Dilatation angle

  ! Pre-compute the denominator for efficiency
  denom = 1.0/(pb%cns_params%a*(sigma+tau*tan_psi))
  mu_star = pb%cns_params%mu_tilde_star

  ! The granular flow strain rate for the case when the shear stress approaches
  ! the shear strength of the grain contacts
  y_gr = pb%cns_params%y_gr_star*exp((tau*(1 - mu_star*tan_psi) - sigma*(mu_star + tan_psi))*denom)

  ! Shear strain rate [1/s] due to viscous flow (pressure solution)
  y_ps = calc_e_ps(tau, theta, .false., .false., pb)
  y_ps_bulk = 0
  if (pb%features%localisation == 1) then
    y_ps_bulk = calc_e_ps(tau, theta2, .false., .true., pb)
  endif
  mu_tilde = calc_mu_tilde(y_gr, pb)           ! Grain-boundary friction

  cohesion = 0
  if (pb%features%cohesion == 1) then
    cos_psi = cos(atan(tan_psi))
    cohesion = (PI/pb%cns_params%H)*(tan_psi*alpha*pb%coh_params%C_star)/(cos_psi*(1-mu_tilde*tan_psi))
  endif

  ! The total slip velocity is the combined contribution of granular flow
  ! and pressure solution (parallel processes)
  v = pb%cns_params%w*(pb%cns_params%lambda*(y_gr + y_ps) + (1-pb%cns_params%lambda)*y_ps_bulk)

end function compute_velocity

!--------------------------------------------------------------------------------------
! SEISMIC: calculate the dilatation angle, which is a geometric consequence
! of the instantaneous porosity (theta)
!--------------------------------------------------------------------------------------
function calc_tan_psi(theta,pb) result(tan_psi)

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: theta
  double precision, dimension(pb%mesh%nn) :: tan_psi

  tan_psi = 2*pb%cns_params%H*(pb%cns_params%phi0 - theta)

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
  mu_tilde = pb%cns_params%mu_tilde_star + pb%cns_params%a * log(abs(y_gr)/pb%cns_params%y_gr_star)

end function calc_mu_tilde

!--------------------------------------------------------------------------------------
! SEISMIC: calculate the rate of pressure solution as a function of effective
! normal or shear stress, and porosity. To prevent porosities from going negative
! compaction rates can be truncated by an error function. The logical 'truncate'
! indicates if truncation is needed. This only applies to compaction in the
! normal direction, and not to shear deformation (viscous flow), hence truncate
! should be false when calculating shear strain rate
! Important NOTE: if the cut-off by the error function is too sharp, then
! the solver may not converge or require extremely high accuracy (ODE becomes stiff)
!--------------------------------------------------------------------------------------
function calc_e_ps(sigma,theta,truncate,bulk,pb) result(e_ps_dot)

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: sigma, theta
  double precision, dimension(pb%mesh%nn) :: e_ps_dot_d, e_ps_dot_s, e_ps_dot
  double precision, dimension(pb%mesh%nn) :: Z_diff, Z_diss1, Z_diss2
  logical, intent(in) :: truncate, bulk

  if (bulk .eqv. .true.) then
    Z_diff = pb%cns_params%IPS_const_diff_bulk
    Z_diss1 = pb%cns_params%IPS_const_diss1_bulk
    Z_diss2 = pb%cns_params%IPS_const_diss2_bulk
  else
    Z_diff = pb%cns_params%IPS_const_diff
    Z_diss1 = pb%cns_params%IPS_const_diss1
    Z_diss2 = pb%cns_params%IPS_const_diss2
  endif

  e_ps_dot_d = Z_diff*sigma/(2*pb%cns_params%phi0 - 2*theta)**2
  e_ps_dot_s =  Z_diss1*Z_diss2*sigma*pb%cns_params%phi0/(pb%cns_params%phi0 - theta)

  e_ps_dot = e_ps_dot_d + e_ps_dot_s

  ! If truncation is required, use an error function to set the strain rate
  ! to zero. A cut-off porosity of about 3% is chosen here, which corresponds
  ! to the percolation limit of aggregates
  if (truncate .eqv. .true.) e_ps_dot = 0.5*(1 + erf(50*(theta-0.06)))*e_ps_dot

end function calc_e_ps

end module friction
