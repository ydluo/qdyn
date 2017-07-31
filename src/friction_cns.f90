!===============================================================================
!
! SEISMIC: CNS microphysical model
!
! Rate-and-state friction law i_rns_law = 3 refers to the CNS (Chen-Niemeijer-
! Spiers) microphysical model, with itheta_law = 3 to diffusion controlled
! pressure solution, and itheta_law = 4 to dissolution controlled pressure
! solution creep. In this interpretation, rate-and-state friction results from
! the interplay between compaction by pressure solution creep and dilation by
! granular flow. At low sliding velocities, pressure solution creep is
! relatively fast compared to granular dilatation, so that the gouge densifies
! (porosity reduces). By contrast, at high sliding velocities, the gouge dilates
! and 'softens', i.e. supports lower shear stress. At a given velocity, a
! steady-state porosity can be attained that corresponds to the porosity at
! which compaction balances dilatation, and hence results in a steady-state
! shear strength. Velocity-weakening behaviour emerges naturally from this
! model. This model is compatible with the thermal pressurisation model
!
! References:
!
!   [NS]  Niemeijer & Spiers (2007), J. Geophys. Res., doi:10.1029/2007JB00500
!   [CS]  Chen & Spiers (2016), JGR: Solid Earth, doi:10.1002/2016JB013470
!   [VdE] Van den Ende et al. (subm.), Tectonophysics
!
! In the view of rate-and-state friction, we take the porosity as the
! (microstructural) state, and the constitutive relations are solved for the
! shear stress (rather than velocity; see derivs_all.f90). The (spatially non-
! uniform) values of the material/gouge properties and initial porosity are
! read from the input script, and stored in pb%cns_params.
!
!===============================================================================

module friction_cns

  use problem_class, only : problem_type

  implicit none
  private

  public :: dphi_dt, compute_velocity, dalpha_dt, CNS_derivs
  private :: calc_tan_psi, calc_mu_tilde, calc_e_ps

contains

!===============================================================================
! SEISMIC: calculate compaction/dilatation rate (dphi/dt)
!===============================================================================
subroutine dphi_dt(v,tau,sigma,phi,phi2,dth_dt,dth2_dt,pb)

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: v, tau, sigma
  double precision, dimension(pb%mesh%nn), intent(in) :: phi, phi2
  double precision, dimension(pb%mesh%nn) :: dth_dt, dth2_dt
  double precision, dimension(pb%mesh%nn) :: tan_psi, e_ps_dot, e_ps_dot_bulk
  double precision, dimension(pb%mesh%nn) :: y_ps, y_gr, y_ps_bulk

  tan_psi = calc_tan_psi(phi, pb)   ! Dilatation angle
  e_ps_dot = calc_e_ps(sigma, phi, .true., .false., pb)   ! Compaction rate
  y_ps = calc_e_ps(tau, phi, .false., .false., pb)   ! Press. soln. shear rate

  e_ps_dot_bulk = 0
  y_ps_bulk = 0

  ! If localisation is requested, calculate the values in the bulk
  if (pb%features%localisation == 1) then
    e_ps_dot_bulk = calc_e_ps(sigma, phi2, .true., .true., pb)
    y_ps_bulk = calc_e_ps(tau, phi2, .false., .true., pb)
  endif

  ! Calculate granular flow rate (assumed to operate only within the localised
  ! zone) as y_gr = y_t - y_ps
  y_gr =  (1.0/pb%cns_params%lambda)*(v/pb%cns_params%L - &
          (1-pb%cns_params%lambda)*y_ps_bulk) - y_ps

  ! The nett rate of change of porosity is given by the relation
  ! phi_dot = -(1 - phi)*strain_rate [CS], where the strain rate equals
  ! the compaction rate by pressure solution (e_ps_dot) minus the dilatation
  ! rate by granular flow (tan_psi*y_gr). Note that compaction is measured
  ! positive, dilatation negative
  dth_dt = (1 - phi)*(tan_psi*y_gr - e_ps_dot)
  dth2_dt = -(1 - phi2)*e_ps_dot_bulk

end subroutine dphi_dt

!===============================================================================
! SEISMIC: time-dependent cohesion, work in progess...
!===============================================================================
function dalpha_dt(v,tau,sigma,phi,alpha,pb) result(da_dt)

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: tau, sigma, phi, alpha, v
  double precision, dimension(pb%mesh%nn) :: da_dt, psi, tan_psi, sin_psi, cos_psi
  double precision, dimension(pb%mesh%nn) :: sigma_i, q, power_term, delta_alpha
  double precision, dimension(pb%mesh%nn) :: y_ps, y_gr, x, y
  double precision, parameter :: small = 1e-15

  write(6,*) "ERROR: The grain-boundary healing feature is in development..."
  stop

  q = 2*(pb%cns_params%phi_c - phi)
  tan_psi = pb%cns_params%H*q
  psi = atan(tan_psi)
  cos_psi = cos(psi)
  sin_psi = tan_psi * cos_psi

  y_ps = calc_e_ps(tau, phi, .false., .false., pb)   ! Press. soln. shear rate
  y_gr = v/pb%cns_params%L - y_ps

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

!===============================================================================
! SEISMIC: the subroutine below was written to optimise/reduce function calls,
! and to improve code handling (less code duplication, hopefully fewer bugs).
!===============================================================================
subroutine CNS_derivs(  v, dth_dt, dth2_dt, dv_dtau, dv_dphi, dv_dP, &
                        tau, sigma, phi, phi2, alpha, pb )

  use constants, only : PI

  ! Input variables
  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: tau, sigma, phi, phi2, alpha

  ! Output variables
  double precision, dimension(pb%mesh%nn) :: v, dth_dt, dth2_dt, dv_dtau, dv_dphi, dv_dP

  ! Internal variables
  double precision, dimension(pb%mesh%nn) :: tan_psi, mu_tilde, denom, dummy_var
  double precision, dimension(pb%mesh%nn) :: mu_star, y_gr, tau_gr
  double precision, dimension(pb%mesh%nn) :: cos_psi, cohesion, f_phi_sqrt, L_sb
  double precision, dimension(pb%mesh%nn) :: e_ps, y_ps, e_ps_bulk, y_ps_bulk
  double precision, dimension(pb%mesh%nn) :: dv_dtau_ps, dv_dphi_ps

  ! Pre-define some internal variables
  tan_psi = calc_tan_psi(phi, pb)               ! Dilatation angle
  mu_star = pb%cns_params%mu_tilde_star         ! Reference GB friction
  L_sb = pb%cns_params%lambda*pb%cns_params%L   ! Thickness of localised zone

  ! Pre-compute the denominator for efficiency
  denom = 1.0/(pb%cns_params%a_tilde*(sigma+tau*tan_psi))
  ! Pre-compute a dummy variable that we will re-use a few times
  dummy_var = (tau*(1 - mu_star*tan_psi) - sigma*(mu_star + tan_psi))*denom

  ! The granular flow strain rate [CS, Eqn. 33b]
  y_gr = pb%cns_params%y_gr_star*exp(dummy_var)

  mu_tilde = calc_mu_tilde(y_gr, pb)           ! Grain-boundary friction

  ! Pre-compute the square root of the porosity function. This will serve as
  ! a basis for calculating the porosity functions
  ! f_phi_sqrt = 1.0/(2*(pb%cns_params%phi_c - phi))

  ! Shear strain rate [1/s] due to viscous flow (pressure solution)
  e_ps = calc_e_ps(sigma, phi, .true., .false., pb)   ! Compaction rate
  y_ps = calc_e_ps(tau, phi, .false., .false., pb)    ! Shear creep rate

  ! If localisation is requested, calculate pressure solution rates
  ! in the bulk zone, otherwise set to zero
  e_ps_bulk = 0
  y_ps_bulk = 0
  if (pb%features%localisation == 1) then
    e_ps_bulk = calc_e_ps(sigma, phi2, .true., .true., pb)
    y_ps_bulk = calc_e_ps(tau, phi2, .false., .true., pb)
  endif

  ! If cohesion is requested, calculate the contribution of cohesion to the
  ! gouge strength. See [NS], Eqn. 20. This feature still in development...
  cohesion = 0
  if (pb%features%cohesion == 1) then
    cos_psi = cos(atan(tan_psi))
    cohesion = (PI/pb%cns_params%H)*(tan_psi*alpha*pb%coh_params%C_star)/(cos_psi*(1-mu_tilde*tan_psi))
  endif

  ! The total slip velocity is the combined contribution of granular flow
  ! and pressure solution (parallel processes)
  v = pb%cns_params%L*(pb%cns_params%lambda*(y_gr + y_ps) + (1-pb%cns_params%lambda)*y_ps_bulk)

  ! The nett rate of change of porosity is given by the relation
  ! phi_dot = -(1 - phi)*strain_rate, where the strain rate equals
  ! the compaction rate by pressure solution (e_ps) minus the dilatation
  ! rate by granular flow (tan_psi*v/w). Note that compaction is measured
  ! positive, dilatation negative
  dth_dt = (1 - phi)*(tan_psi*y_gr - e_ps)
  dth2_dt = -(1 - phi2)*e_ps_bulk

  ! Next, calculate the partial derivatives dv_dtau and dv_dphi for the
  ! pressure solution creep. The functional form depends on the rate-limiting
  ! mechanism (diffusion, dissolution, or precipitation)
  ! dv/dphi_d = 2 dv/dphi_s, but since pressure solution contributes
  ! negligibly to radiation damping, we ignore this factor 2

  dv_dtau_ps = y_ps/tau
  dv_dphi_ps = y_ps/(pb%cns_params%phi_c - phi)

  ! The partial derivatives dv_dtau and dv_dphi of the overall slip velocity
  ! are the sum of the partial derivatives of the pressure solution and
  ! granular flow velocities.

  dv_dtau = L_sb*(dv_dtau_ps + y_gr*(1-mu_tilde*tan_psi)*denom)
  dv_dphi = L_sb*(dv_dphi_ps + y_gr*2*pb%cns_params%H*(sigma + mu_tilde*tau)*denom )

  ! When thermal pressurisation is requested, update dV_dP
  if (pb%features%tp == 1) then
    dv_dP = pb%cns_params%L*y_gr*tau*(1 + tan_psi**2)*pb%cns_params%a_tilde*denom**2
  endif

end subroutine CNS_derivs

!===============================================================================
! SEISMIC: calculate the instantaneous fault slip velocity, resulting from
! parallel operation of granular flow and pressure solution creep
!===============================================================================
function compute_velocity(tau,sigma,phi,phi2,alpha,pb) result(v)

  use constants, only : PI

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: tau, sigma, phi
  double precision, dimension(pb%mesh%nn), intent(in) :: phi2, alpha
  double precision, dimension(pb%mesh%nn) :: tan_psi, denom, mu_star
  double precision, dimension(pb%mesh%nn) :: y_ps, y_gr, tau_gr
  double precision, dimension(pb%mesh%nn) :: y_ps_bulk, v

  ! Dilatation angle
  tan_psi = calc_tan_psi(phi, pb)

  ! Pre-compute the denominator for efficiency
  denom = 1.0/(pb%cns_params%a_tilde*(sigma+tau*tan_psi))
  mu_star = pb%cns_params%mu_tilde_star

  ! The granular flow strain rate [CS, Eqn. 33b]
  y_gr = pb%cns_params%y_gr_star*exp((tau*(1 - mu_star*tan_psi) - sigma*(mu_star + tan_psi))*denom)

  ! Shear strain rate [1/s] due to viscous flow (pressure solution)
  y_ps = calc_e_ps(tau, phi, .false., .false., pb)
  y_ps_bulk = 0
  if (pb%features%localisation == 1) then
    y_ps_bulk = calc_e_ps(tau, phi2, .false., .true., pb)
  endif

  ! The total slip velocity is the combined contribution of granular flow
  ! and pressure solution (parallel processes)
  v = pb%cns_params%L*(pb%cns_params%lambda*(y_gr + y_ps) + (1-pb%cns_params%lambda)*y_ps_bulk)

end function compute_velocity

!===============================================================================
! SEISMIC: calculate the dilatation angle, which is a geometric consequence
! of the instantaneous porosity (phi). See [NS], Eqn. 2
!===============================================================================
function calc_tan_psi(phi,pb) result(tan_psi)

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: phi
  double precision, dimension(pb%mesh%nn) :: tan_psi

  tan_psi = 2*pb%cns_params%H*(pb%cns_params%phi_c - phi)

end function calc_tan_psi

!===============================================================================
! SEISMIC: calculate the grain-boundary coefficient of friction, which is
! logarithmically rate-strengthening, assuming some atomic scale jump process
! See [CS], Section 3.5
!===============================================================================
function calc_mu_tilde(y_gr,pb) result(mu_tilde)

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: y_gr
  double precision, dimension(pb%mesh%nn) :: mu_tilde

  ! NOTE: the granular flow strain rate should be used here,
  ! not the total strain rate
  mu_tilde =  pb%cns_params%mu_tilde_star + pb%cns_params%a_tilde * &
              log(abs(y_gr)/pb%cns_params%y_gr_star)

end function calc_mu_tilde

!===============================================================================
! SEISMIC: calculate the rate of pressure solution as a function of effective
! normal or shear stress, and porosity. To prevent porosities from going
! negative compaction rates can be truncated by a cut-off porosity. The logical
! 'truncate' indicates if truncation is needed. This only applies to compaction
! in the normal direction, and not to shear deformation (viscous flow), hence
! truncate should be false when calculating shear strain rate
!===============================================================================
function calc_e_ps(sigma,phi,truncate,bulk,pb) result(e_ps_dot)

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: sigma, phi
  double precision, dimension(pb%mesh%nn) :: e_ps_dot_d, e_ps_dot_s, e_ps_dot
  double precision, dimension(pb%mesh%nn) :: Z_diff, Z_diss1, Z_diss2
  logical, intent(in) :: truncate, bulk

  ! Check if the pressure soltuion rate is calculated for the localised
  ! zone, or for the bulk. As the grain size may be different between these
  ! zones, different kinetics may be supplied
  if (bulk .eqv. .true.) then
    Z_diff = pb%cns_params%IPS_const_diff_bulk
    Z_diss1 = pb%cns_params%IPS_const_diss1_bulk
    Z_diss2 = pb%cns_params%IPS_const_diss2_bulk
  else
    Z_diff = pb%cns_params%IPS_const_diff
    Z_diss1 = pb%cns_params%IPS_const_diss1
    Z_diss2 = pb%cns_params%IPS_const_diss2
  endif

  e_ps_dot_d = Z_diff*sigma/(2*pb%cns_params%phi_c - 2*phi)**2
  e_ps_dot_s =  Z_diss1*Z_diss2*sigma*pb%cns_params%phi_c/(pb%cns_params%phi_c - phi)

  e_ps_dot = e_ps_dot_d + e_ps_dot_s

  ! The porosity cut-off should only apply to volumetric strain rates (i.e.
  ! compaction, not for shear). Phi0 can be taken as the percolation threshold
  ! See [VdE], Eqn. 6 (NOTE: update this when paper is accepted)
  if (truncate .eqv. .true.) then
    e_ps_dot = e_ps_dot * (phi - pb%cns_params%phi0)/pb%cns_params%phi_c
  endif

end function calc_e_ps

end module friction_cns
