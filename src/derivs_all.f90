!Module for derivs

module derivs_all

  implicit none

  private

  public :: derivs


contains


!====================Subroutine derivs===========================
!-------Compute theta, vslip and time deriv for one time step----
!-------Call rdft in compute_stress : real DFT-------------------
!----------------------------------------------------------------

subroutine derivs(time,yt,dydt,pb)

  use problem_class
  use fault_stress, only : compute_stress
  use friction, only : dtheta_dt, dalpha_dt, dmu_dv_dtheta, compute_velocity

  type(problem_type), intent(inout) :: pb
  double precision, intent(in) :: time, yt(pb%neqs*pb%mesh%nn)
  double precision, intent(out) :: dydt(pb%neqs*pb%mesh%nn)

  double precision, dimension(pb%mesh%nn) :: theta, sigma, alpha, tau, v
  double precision, dimension(pb%mesh%nn) :: dsigma_dt, dtau_dt
  double precision, dimension(pb%mesh%nn) :: dmu_dv, dmu_dtheta
  double precision :: dtau_per

  integer :: ind_stress_coupling, ind_cohesion

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

  ! SEISMIC: start unpacking values of yt and dydt
  ! Note that the values of tau, sigma, etc. are only set if they're being
  ! used. They just serve as dummy variables otherwise
  theta = yt(1::pb%neqs)

  if (pb%features%stress_coupling == 1) then
    sigma = yt(ind_stress_coupling::pb%neqs)
    dsigma_dt = dydt(ind_stress_coupling::pb%neqs)
  else
    sigma = pb%sigma
  endif

  if (pb%features%cohesion == 1) then
    alpha = yt(ind_cohesion::pb%neqs)
  endif

  if (pb%i_rns_law == 3) then
    ! SEISMIC: the CNS model is solved for stress, not for velocity, so we
    ! compute the velocity and use that to compute dtau_dt, which is stored
    ! in dydt(2::pb%neqs). tau is stored as yt(2::pb%neqs). The system of
    ! ordinary differential equations is then solved in the same way as with
    ! the rate-and-state formulation.
    ! NOTE: from what I understand from the comments in solver.f90, this method
    ! is a temporary solution. I put a comment in solver.f90 to refer to here
    ! in the case that this temporary fix receive a more permanent upgrade

    ! SEISMIC: if the CNS model is selected, compute the slip velocity as
    ! a function of the (current) state of stress and porosity
    tau = yt(2::pb%neqs)
    dtau_dt = dydt(2::pb%neqs)
    v = compute_velocity(tau, sigma, theta, alpha, pb)
  else
    ! SEISMIC: for the classical rate-and-state model formulation, the slip
    ! velocity is stored in yt(2::pb%neqs), and tau is not used (set to zero)
    v = yt(2::pb%neqs)
  endif

  ! compute shear stress rate from elastic interactions, for 0D, 1D & 2D
  ! SEISMIC: note the use of v instead of yt(2::pb%neqs)
  call compute_stress(dtau_dt,dsigma_dt,pb%kernel,v-pb%v_star)

  ! JPA Coulomb
  ! v = 0d0
  ! v(pb%rs_nodes) = yt(2::pb%neqs)
  !call compute_stress(pb%dtau_dt,pb%kernel,v-pb%v_star)

  !YD we may want to modify this part later to be able to
  !impose more complicated loading/pertubation
  !functions involved: problem_class/problem_type; input/read_main
  !                    initialize/init_field;  derivs_all/derivs
  ! periodic loading
  dtau_per = pb%Omper * pb%Aper * cos(pb%Omper*time)

  ! SEISMIC: start repacking values of dydt
  ! Question: I don't think the values of pb%dsigma_dt, pb%dtau_dt etc. need
  ! to be set here

  ! state evolution law, dtheta/dt = f(v,theta)
  ! SEISMIC: in the CNS formulation, theta is the gouge porosity
  dydt(1::pb%neqs) = dtheta_dt(v,tau,sigma,theta,pb)

  if (pb%features%stress_coupling == 1) then
    dydt(ind_stress_coupling::pb%neqs) = dsigma_dt
  endif

  if (pb%features%cohesion == 1) then
    dydt(ind_cohesion::pb%neqs) = dalpha_dt(v,tau,sigma,theta,alpha,pb)
  endif

  if (pb%i_rns_law == 3) then
    ! SEISMIC: the total dtau_dt results from the slip deficit and
    ! periodic loading, which is stored in dydt(2::pb%neqs)
    dydt(2::pb%neqs) = dtau_dt + dtau_per
  else
    ! SEISMIC: the rate-and-state formulation computes the the time-derivative
    ! of velocity, rather than stress

    ! Time derivative of the elastic equilibrium equation
    !  dtau_load/dt + dtau_elastostatic/dt -impedance*dv/dt = sigma*( dmu/dv*dv/dt + dmu/dtheta*dtheta/dt )
    ! Rearranged in the following form:
    !  dv/dt = ( dtau_load/dt + dtau_elastostatic/dt - sigma*dmu/dtheta*dtheta/dt )/( sigma*dmu/dv + impedance )
    ! SEISMIC NOTE: use of pb%sigma, even when pb%neqs == 3? Replaced with sigma (defined above)
    call dmu_dv_dtheta(dmu_dv,dmu_dtheta,v,theta,pb)
    dydt(2::pb%neqs) = ( dtau_per + dtau_dt - sigma*dmu_dtheta*dydt(1::pb%neqs) ) &
                     / ( sigma*dmu_dv + pb%zimpedance )
   endif

end subroutine derivs

end module derivs_all
