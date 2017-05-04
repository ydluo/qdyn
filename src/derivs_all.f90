!Module for derivs

module derivs_all

  use problem_class
  implicit none

  type(problem_type), pointer :: odepb

  private

  public :: derivs, derivs_lsoda, derivs_rk45, jac_lsoda
  public :: odepb


contains


!====================Subroutine derivs===========================
!-------Compute theta, vslip and time deriv for one time step----
!-------Call rdft in compute_stress : real DFT-------------------
!----------------------------------------------------------------

subroutine derivs(time,yt,dydt,pb)

  use fault_stress, only : compute_stress
  use friction, only : dtheta_dt, dalpha_dt, dmu_dv_dtheta, compute_velocity, CNS_derivs

  type(problem_type), intent(inout) :: pb
  double precision, intent(in) :: time, yt(pb%neqs*pb%mesh%nn)
  double precision, intent(out) :: dydt(pb%neqs*pb%mesh%nn)

  double precision, dimension(pb%mesh%nn) :: theta, theta2, sigma, alpha, tau, v
  double precision, dimension(pb%mesh%nn) :: dsigma_dt, dtau_dt, dth_dt, dth2_dt
  double precision, dimension(pb%mesh%nn) :: dmu_dv, dmu_dtheta
  double precision :: dtau_per

  integer :: ind_stress_coupling, ind_cohesion, ind_localisation

  ! storage conventions:
  !
  ! theta = yt(1::pb%neqs)
  ! v = yt(2::pb%neqs)
  ! sigma = yt(ind_stress_coupling::pb%neqs)
  ! alpha = yt(ind_cohesion::pb%neqs)
  !
  ! dtheta/dt = dydt(1::pb%neqs)
  ! dv/dt = dydt(2::pb%neqs)
  ! dsigma/dt = dydt(ind_stress_coupling::pb%neqs)
  ! dalpha/dt = dydt(ind_cohesion::pb%neqs)

  ! SEISMIC: define the indices of yt and dydt based on which
  ! features are requested (defined in input file)
  ind_stress_coupling = 2 + pb%features%stress_coupling
  ind_cohesion = ind_stress_coupling + pb%features%cohesion
  ind_localisation = ind_cohesion + pb%features%localisation

  ! SEISMIC: start unpacking values of yt and dydt
  ! Note that the values of tau, sigma, etc. are only set if they're being
  ! used. They just serve as dummy variables otherwise
  theta = yt(1::pb%neqs)
  theta2 = 0

  if (pb%features%stress_coupling == 1) then
    sigma = yt(ind_stress_coupling::pb%neqs)
    dsigma_dt = dydt(ind_stress_coupling::pb%neqs)
  else
    sigma = pb%sigma
  endif

  if (pb%features%cohesion == 1) then
    alpha = yt(ind_cohesion::pb%neqs)
  endif

  if (pb%features%localisation == 1) then
    theta2 = yt(ind_localisation::pb%neqs)
  endif

  if (pb%i_rns_law == 3) then
    ! SEISMIC: the CNS model is solved for stress, not for velocity, so we
    ! compute the velocity and use that to compute dtau_dt, which is stored
    ! in dydt(2::pb%neqs). tau is stored as yt(2::pb%neqs). The system of
    ! ordinary differential equations is then solved in the same way as with
    ! the rate-and-state formulation.

    tau = yt(2::pb%neqs)

    ! The subroutine below calculates the slip velocity, time-derivatives of
    ! the porosity (theta), and partial derivatives required for radiation
    ! damping. These operations are combined into one subroutine for efficiency
    call CNS_derivs(v, dth_dt, dth2_dt, dmu_dv, dmu_dtheta, tau, sigma, theta, theta2, alpha, pb)
  else
    ! SEISMIC: for the classical rate-and-state model formulation, the slip
    ! velocity is stored in yt(2::pb%neqs), and tau is not used (set to zero)
    v = yt(2::pb%neqs)

    ! SEISMIC: calculate time-derivative of theta
    call dtheta_dt(v,tau,sigma,theta,theta2,dth_dt,dth2_dt,pb)
  endif

  ! compute shear stress rate from elastic interactions, for 0D, 1D & 2D
  ! SEISMIC: note the use of v instead of yt(2::pb%neqs)
  call compute_stress(dtau_dt,dsigma_dt,pb%kernel,v-pb%v_star)

  !YD we may want to modify this part later to be able to
  !impose more complicated loading/pertubation
  !functions involved: problem_class/problem_type; input/read_main
  !                    initialize/init_field;  derivs_all/derivs
  ! periodic loading
  dtau_per = pb%Omper * pb%Aper * cos(pb%Omper*time)

  ! SEISMIC: start repacking values of dydt

  ! state evolution law, dtheta/dt = f(v,theta)
  ! SEISMIC: in the CNS formulation, theta is the gouge porosity
  dydt(1::pb%neqs) = dth_dt

  ! SEISMIC: pack normal stress changes if needed
  if (pb%features%stress_coupling == 1) then
    dydt(ind_stress_coupling::pb%neqs) = dsigma_dt
  endif

  ! SEISMIC: pack rate of grain-boundary healing, if needed
  if (pb%features%cohesion == 1) then
    dydt(ind_cohesion::pb%neqs) = dalpha_dt(v,tau,sigma,theta,alpha,pb)
  endif

  ! SEISMIC: if localisation is included, pack porosity change of bulk
  if (pb%features%localisation == 1) then
    dydt(ind_localisation::pb%neqs) = dth2_dt
  endif

  ! SEISMIC: calculate radiation damping. For rate-and-state, dmu_dv and
  ! dmu_dtheta contain the partial derivatives of friction to V and theta,
  ! respectively. For the CNS model, these variables contain the partials of
  ! velocity to shear stress (dV/dtau) and velocity to porosity (dV/dtheta),
  ! respectively. For compatibility, the names of these variables are not
  ! changed.

  if (pb%i_rns_law == 3) then
    ! SEISMIC: the total dtau_dt results from the slip deficit and
    ! periodic loading, which is stored in dydt(2::pb%neqs)
    ! Damping is included from rewriting the following expression:
    ! dtau/dt = k(Vlp - Vs) - eta*(dV/dtau * dtau/dt + dV/dtheta * dtheta/dt)
    ! Rearrangement gives:
    ! dtau/dt = ( k[Vlp - Vs] - eta*dV/dtheta * dtheta/dt)/(1 + eta*dV/dtau)
    dydt(2::pb%neqs) = (dtau_dt + dtau_per - pb%zimpedance*dmu_dtheta*dth_dt)/(1 + pb%zimpedance*dmu_dv)
  else
    ! SEISMIC: the rate-and-state formulation computes the the time-derivative
    ! of velocity, rather than stress
    call dmu_dv_dtheta(dmu_dv,dmu_dtheta,v,tau,sigma,theta,theta2,pb)

    ! Time derivative of the elastic equilibrium equation
    !  dtau_load/dt + dtau_elastostatic/dt -impedance*dv/dt = sigma*( dmu/dv*dv/dt + dmu/dtheta*dtheta/dt )
    ! Rearranged in the following form:
    !  dv/dt = ( dtau_load/dt + dtau_elastostatic/dt - sigma*dmu/dtheta*dtheta/dt )/( sigma*dmu/dv + impedance )

    dydt(2::pb%neqs) = ( dtau_per + dtau_dt - sigma*dmu_dtheta*dydt(1::pb%neqs) ) &
                     / ( sigma*dmu_dv + pb%zimpedance )
   endif

end subroutine derivs

!--------------------------------------------------------------------------------------
! SEISMIC: the subroutine derivs_lsoda is a wrapper that interfaces between derivs
! and the LSODA solver routine
!--------------------------------------------------------------------------------------
subroutine derivs_lsoda(neq, time, yt, dydt)
  ! NOTE: neq = pb%neq * pb%mesh%nn
  integer :: neq, i
  double precision :: time
  double precision :: yt(neq)
  double precision :: dydt(neq)

  call derivs(time,yt,dydt,odepb)

end subroutine derivs_lsoda

!--------------------------------------------------------------------------------------
! SEISMIC: jacobian for lsoda
!--------------------------------------------------------------------------------------
subroutine jac_lsoda(neq, time, yt, ml, mu, pd, nrpd)
integer :: neq, ml, mu, nrpd
double precision :: time, yt(neq), pd(nrpd,neq)

end subroutine jac_lsoda


!--------------------------------------------------------------------------------------
! SEISMIC: the subroutine derivs_rk45 is a wrapper that interfaces between derivs
! and the Runge-Kutta-Fehlberg solver routine
!--------------------------------------------------------------------------------------
subroutine derivs_rk45(time, yt, dydt)
  double precision :: time
  double precision :: yt(*)
  double precision :: dydt(*)

  call derivs(time,yt,dydt,odepb)

end subroutine derivs_rk45

end module derivs_all
