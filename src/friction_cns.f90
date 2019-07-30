!===============================================================================
!
! SEISMIC: CNS microphysical model
!
! Rate-and-state friction law i_rns_law = 3 refers to the CNS (Chen-Niemeijer-
! Spiers) microphysical model.  In this interpretation, rate-and-state friction
! results from the interplay between compaction by one or more creep mechanisms
! (such as pressure solution) and dilation by granular flow. At low sliding
! velocities, the contribution from the creep mechanism(s) to the volumetric
! strain is relatively fast compared to granular dilatation, so that the gouge
! densifies (porosity reduces). By contrast, at high sliding velocities, the
! gouge dilates and 'softens', i.e. supports lower shear stress. At a given
! velocity, a steady-state porosity can be attained that corresponds to the
! porosity at which compaction balances dilatation, and hence results in a
! steady-state shear strength. Velocity-weakening behaviour emerges naturally
! from this model. This model is compatible with the thermal pressurisation
! model, and with localisation
!
! References:
!
!   [NS]  Niemeijer & Spiers (2007), J. Geophys. Res., doi:10.1029/2007JB00500
!   [CS]  Chen & Spiers (2016), JGR: Solid Earth, doi:10.1002/2016JB013470
!   [VdE] van den Ende et al. (2018), Tecto, doi:10.1016/j.tecto.2017.11.040
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

  public :: dphi_dt, compute_velocity, CNS_derivs

contains

!===============================================================================
! SEISMIC: calculate compaction/dilatation rate (dphi/dt)
!===============================================================================
subroutine dphi_dt(v,tau,sigma,phi,phi2,dth_dt,dth2_dt,pb)

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: v, tau, sigma
  double precision, dimension(pb%mesh%nn), intent(in) :: phi, phi2
  double precision, dimension(pb%mesh%nn) :: dth_dt, dth2_dt
  double precision, dimension(pb%mesh%nn) :: tan_psi, y_dot_gr
  double precision, dimension(pb%mesh%nn) :: e_dot_creep, y_dot_creep
  double precision, dimension(pb%mesh%nn) :: e_dot_creep_bulk, y_dot_creep_bulk

  tan_psi = calc_tan_psi(phi, pb)   ! Dilatation angle
  e_dot_creep = calc_creep_rate(sigma, phi, .true., .false., pb)   ! Compaction rate
  y_dot_creep = calc_creep_rate(tau, phi, .false., .false., pb)   ! Press. soln. shear rate

  e_dot_creep_bulk = 0
  y_dot_creep_bulk = 0

  ! If localisation is requested, calculate the values in the bulk
  if (pb%features%localisation == 1) then
    e_dot_creep_bulk = calc_creep_rate(sigma, phi2, .true., .true., pb)
    y_dot_creep_bulk = calc_creep_rate(tau, phi2, .false., .true., pb)
  endif

  ! Calculate granular flow rate (assumed to operate only within the localised
  ! zone) as y_dot_gr = y_t - y_dot_creep
  y_dot_gr =  (1.0/pb%cns_params%lambda)*(v/pb%cns_params%L - &
          (1-pb%cns_params%lambda)*y_dot_creep_bulk) - y_dot_creep

  ! The nett rate of change of porosity is given by the relation
  ! phi_dot = -(1 - phi)*strain_rate [CS], where the strain rate equals
  ! the compaction rate by creep (e_dot_creep) minus the dilatation
  ! rate by granular flow (tan_psi*y_dot_gr). Note that compaction is measured
  ! positive, dilatation negative
  dth_dt = (1 - phi)*(tan_psi*y_dot_gr - e_dot_creep)
  dth2_dt = -(1 - phi2)*e_dot_creep_bulk

end subroutine dphi_dt

!===============================================================================
! SEISMIC: the subroutine below was written to optimise/reduce function calls,
! and to improve code handling (less code duplication, hopefully fewer bugs).
!===============================================================================
subroutine CNS_derivs(  v, dth_dt, dth2_dt, dv_dtau, dv_dphi, dv_dP, &
                        tau, sigma, phi, phi2, pb )

  use constants, only : PI

  ! Input variables
  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: tau, sigma, phi, phi2

  ! Output variables
  double precision, dimension(pb%mesh%nn) :: v, dth_dt, dth2_dt
  double precision, dimension(pb%mesh%nn) :: dv_dtau, dv_dphi, dv_dP

  ! Internal variables
  double precision, dimension(pb%mesh%nn) :: tan_psi, denom, dummy_var
  double precision, dimension(pb%mesh%nn) :: mu_star, y_dot_gr, L_sb
  double precision, dimension(pb%mesh%nn) :: e_dot_creep, y_dot_creep
  double precision, dimension(pb%mesh%nn) :: e_dot_creep_bulk, y_dot_creep_bulk
  double precision, dimension(pb%mesh%nn) :: dv_dtau_creep, dv_dphi_creep
  integer :: N_creep

  N_creep = pb%cns_params%N_creep

  ! Pre-define some internal variables
  tan_psi = calc_tan_psi(phi, pb)               ! Dilatation angle
  mu_star = pb%cns_params%mu_tilde_star         ! Reference GB friction
  L_sb = pb%cns_params%lambda*pb%cns_params%L   ! Thickness of localised zone

  ! Pre-compute the denominator for efficiency
  denom = 1d0/(pb%cns_params%a_tilde*(sigma+tau*tan_psi))
  ! Pre-compute a dummy variable that we will re-use a few times
  dummy_var = (tau*(1 - mu_star*tan_psi) - sigma*(mu_star + tan_psi))*denom

  ! The granular flow strain rate [CS, Eqn. 33b]
  y_dot_gr = pb%cns_params%y_gr_star*exp(dummy_var)

  ! Compute strain rates and derivatives associated with N creep mechanisms
  call calc_creep_rates_derivs( sigma, tau, phi, e_dot_creep, y_dot_creep, &
                                dv_dtau_creep, dv_dphi_creep, pb)

  ! If localisation is requested, calculate the creep rates
  ! in the bulk zone, otherwise set to zero
  e_dot_creep_bulk = 0d0
  y_dot_creep_bulk = 0d0
  if (pb%features%localisation == 1) then
    e_dot_creep_bulk = calc_creep_rate(sigma, phi2, .true., .true., pb)
    y_dot_creep_bulk = calc_creep_rate(tau, phi2, .false., .true., pb)
  endif

  ! The total slip velocity is the combined contribution of granular flow
  ! and shear creep (parallel processes)
  v = pb%cns_params%L*(pb%cns_params%lambda*(y_dot_gr + y_dot_creep) + &
      (1-pb%cns_params%lambda)*y_dot_creep_bulk)

  ! The nett rate of change of porosity is given by the relation
  ! phi_dot = -(1 - phi)*strain_rate, where the strain rate equals
  ! the compaction rate by creep (e_dot_creep) minus the dilatation
  ! rate by granular flow (tan_psi*v/w). Note that compaction is measured
  ! positive, dilatation negative
  dth_dt = (1 - phi)*(tan_psi*y_dot_gr - e_dot_creep)
  dth2_dt = -(1 - phi2)*e_dot_creep_bulk

  ! The partial derivatives dv_dtau and dv_dphi of the overall slip velocity
  ! are the sum of the partial derivatives of the creep mechanisms and
  ! granular flow velocities.
  dv_dtau = L_sb*(dv_dtau_creep + &
            y_dot_gr*pb%cns_params%a_tilde*sigma*(1 + tan_psi**2)*denom**2)
  dv_dphi = L_sb*(dv_dphi_creep + &
            y_dot_gr*2*pb%cns_params%H*pb%cns_params%a_tilde*(sigma**2 + tau**2)*denom**2)

  ! When thermal pressurisation is requested, update dV_dP
  if (pb%features%tp == 1) then
    dv_dP = pb%cns_params%L*y_dot_gr*tau*(1 + tan_psi**2)*pb%cns_params%a_tilde*denom**2
  endif

end subroutine CNS_derivs

!===============================================================================
! TODO: update description
! SEISMIC: the subroutine below was written to optimise/reduce function calls,
! and to improve code handling (less code duplication, hopefully fewer bugs).
!===============================================================================
subroutine calc_creep_rates_derivs( sigma, tau, phi, e_dot_creep_all, &
                                    y_dot_creep_all, dv_dtau_creep_all, &
                                    dv_dphi_creep_all, pb)

  ! Input variables
  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: sigma, tau, phi

  ! Output variables
  double precision, dimension(pb%mesh%nn) :: e_dot_creep_all, y_dot_creep_all
  double precision, dimension(pb%mesh%nn) :: dv_dtau_creep_all, dv_dphi_creep_all

  ! Internal variables
  double precision, dimension(pb%mesh%nn) :: e_dot_creep, y_dot_creep
  double precision, dimension(pb%mesh%nn) :: dv_dtau_creep, dv_dphi_creep
  double precision, dimension(pb%mesh%nn) :: f_phi1, f_phi2, A, n, m
  double precision, dimension(pb%mesh%nn) :: inv_phi, inv_tau
  integer :: i, N_creep

  N_creep = pb%cns_params%N_creep
  y_dot_creep_all = 0d0
  e_dot_creep_all = 0d0
  dv_dtau_creep_all = 0d0
  dv_dphi_creep_all = 0d0

  ! Precompute expensive variables
  inv_phi = 1d0 / (pb%cns_params%phi_c - phi)
  inv_tau = 1d0 / tau

  ! The porosity function for shear creep
  f_phi1 = pb%cns_params%phi_c * inv_phi

  ! The porosity function for compaction creep
  ! The porosity cut-off should only apply to volumetric strain rates (i.e.
  ! compaction, not for shear). Phi0 can be taken as the percolation threshold
  ! See [VdE], Eqn. 6
  f_phi2 = (phi - pb%cns_params%phi0) * inv_phi

  do i=1,N_creep
    A = pb%cns_params%A(i::N_creep)
    n = pb%cns_params%n(i::N_creep)
    m = pb%cns_params%m(i::N_creep)

    y_dot_creep = A * (tau * f_phi1**m)**n
    e_dot_creep = A * (sigma * f_phi2**m)**n

    dv_dtau_creep = n * y_dot_creep * inv_tau
    dv_dphi_creep = m * n * y_dot_creep * inv_phi

    y_dot_creep_all = y_dot_creep_all + y_dot_creep
    e_dot_creep_all = e_dot_creep_all + e_dot_creep

    dv_dtau_creep_all = dv_dtau_creep_all + dv_dtau_creep
    dv_dphi_creep_all = dv_dphi_creep_all + dv_dphi_creep
  enddo

  ! TODO: add components of bulk?

end subroutine calc_creep_rates_derivs

!===============================================================================
! SEISMIC: calculate the instantaneous fault slip velocity, resulting from
! parallel operation of granular flow and N creep mechanisms
!===============================================================================
function compute_velocity(tau,sigma,phi,phi2,pb) result(v)

  use constants, only : PI

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: tau, sigma, phi
  double precision, dimension(pb%mesh%nn), intent(in) :: phi2
  double precision, dimension(pb%mesh%nn) :: tan_psi, denom, mu_star
  double precision, dimension(pb%mesh%nn) :: y_dot_creep, y_dot_gr
  double precision, dimension(pb%mesh%nn) :: y_dot_creep_bulk, v

  ! Dilatation angle
  tan_psi = calc_tan_psi(phi, pb)

  ! Pre-compute the denominator for efficiency
  denom = 1.0/(pb%cns_params%a_tilde*(sigma+tau*tan_psi))
  mu_star = pb%cns_params%mu_tilde_star

  ! The granular flow strain rate [CS, Eqn. 33b]
  y_dot_gr = pb%cns_params%y_gr_star*exp((tau*(1 - mu_star*tan_psi) - &
            sigma*(mu_star + tan_psi))*denom)

  ! Shear strain rate [1/s] due to shear creep
  y_dot_creep = calc_creep_rate(tau, phi, .false., .false., pb)
  y_dot_creep_bulk = 0d0
  if (pb%features%localisation == 1) then
    y_dot_creep_bulk = calc_creep_rate(tau, phi2, .false., .true., pb)
  endif

  ! The total slip velocity is the combined contribution of granular flow
  ! and creep (parallel processes)
  v = pb%cns_params%L*(pb%cns_params%lambda*(y_dot_gr + y_dot_creep) + &
      (1-pb%cns_params%lambda)*y_dot_creep_bulk)

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
! NOTE: this function is currently not used, but kept here for reference
!===============================================================================
! function calc_mu_tilde(y_dot_gr,pb) result(mu_tilde)
!
!   type(problem_type), intent(in) :: pb
!   double precision, dimension(pb%mesh%nn), intent(in) :: y_dot_gr
!   double precision, dimension(pb%mesh%nn) :: mu_tilde
!
!   ! NOTE: the granular flow strain rate should be used here,
!   ! not the total strain rate
!   mu_tilde =  pb%cns_params%mu_tilde_star + pb%cns_params%a_tilde * &
!               log(abs(y_dot_gr)/pb%cns_params%y_gr_star)
!
! end function calc_mu_tilde

!===============================================================================
! SEISMIC: calculate the rate of shear/normal creep as a function of effective
! normal or shear stress, and porosity. To prevent porosities from going
! negative compaction rates can be truncated by a cut-off porosity. The logical
! 'truncate' indicates if truncation is needed. This only applies to compaction
! in the normal direction, and not to shear deformation (viscous flow), hence
! truncate should be false when calculating shear strain rate
!===============================================================================
function calc_creep_rate(sigma,phi,truncate,bulk,pb) result(e_dot_creep)

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn) :: sigma, phi, f_phi
  double precision, dimension(pb%mesh%nn) :: e_dot_creep
  double precision, dimension(pb%mesh%nn) :: A, n, m
  logical, intent(in) :: truncate, bulk
  integer :: i, N_creep

  N_creep = pb%cns_params%N_creep
  e_dot_creep = 0d0

  if (truncate .eqv. .false.) then
    ! The porosity function for shear creep
    f_phi = pb%cns_params%phi_c / (pb%cns_params%phi_c - phi)
  else
    ! The porosity function for compaction creep
    ! The porosity cut-off should only apply to volumetric strain rates (i.e.
    ! compaction, not for shear). Phi0 can be taken as the percolation threshold
    ! See [VdE], Eqn. 6
    f_phi = (phi - pb%cns_params%phi0) / (pb%cns_params%phi_c - phi)
  endif

  ! Check if the creep rate is calculated for the localised zone, or for the
  ! bulk. As e.g. the grain size may be different between these zones, different
  ! flow law parameters may be supplied
  if (bulk .eqv. .false.) then
    do i=1,N_creep
      A = pb%cns_params%A(i::N_creep)
      n = pb%cns_params%n(i::N_creep)
      m = pb%cns_params%m(i::N_creep)
      e_dot_creep = e_dot_creep + A * (sigma * f_phi**m)**n
    enddo
  else
    do i=1,N_creep
      A = pb%cns_params%A_bulk(i::N_creep)
      n = pb%cns_params%n_bulk(i::N_creep)
      m = pb%cns_params%m_bulk(i::N_creep)
      e_dot_creep = e_dot_creep + A * (sigma * f_phi**m)**n
    enddo
  endif

end function calc_creep_rate

end module friction_cns
