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
  use friction, only : dtheta_dt, dmu_dv_dtheta, compute_velocity

  type(problem_type), intent(inout) :: pb
  double precision, intent(in) :: time, yt(pb%neqs*pb%mesh%nn)
  double precision, intent(out) :: dydt(pb%neqs*pb%mesh%nn)

  double precision, dimension(pb%mesh%nn) :: dmu_dv, dmu_dtheta, sigma
  double precision :: dtau_per

  double precision, dimension(pb%mesh%nn) :: v

  ! storage conventions:
  ! v = yt(2::pb%neqs)
  ! theta = yt(1::pb%neqs)
  ! dv/dt = dydt(2::pb%neqs)
  ! dtheta/dt = dydt(1::pb%neqs)
  ! SEISMIC: dydt(3::pb%neqs) = dsigma/dt ? (normal stress)

  ! SEISMIC: Chen's law is solved for stress, not for velocity, so we
  ! compute the velocity and use that to compute dtau_dt, which is stored
  ! in dydt(2::pb%neqs). tau is stored as yt(2::pb%neqs). The system of
  ! ordinary differential equations is then solved in the same way as with
  ! the rate-and-state formulation.
  if (pb%i_rns_law == 3) then
    ! SEISMIC: in the case of normal stress coupling (whatever that may be)
    ! use the current value of sigma. Chen's model depends explicitly on sigma
    ! NOTE: from what I understand from the comments in solver.f90, this method
    ! is a temporary solution. I put a comment in solver.f90 to refer to here
    ! in the case that this temporary fix receive a more permanent upgrade
    if ( pb%neqs == 3) then
      sigma = yt(3::pb%neqs)
    else
      sigma = pb%sigma
    endif
    ! SEISMIC: if Chen's model is selected, compute the slip velocity as
    ! a function of the (current) state of stress and porosity
    v = compute_velocity(yt(2::pb%neqs), sigma, yt(1::pb%neqs), pb)
  else
    ! SEISMIC: for the classical rate-and-state model formulation, the slip
    ! velocity is stored in yt(2::pb%neqs)
    v = yt(2::pb%neqs)
  endif

  ! compute shear stress rate from elastic interactions, for 0D, 1D & 2D
  ! SEISMIC: note the use of v instead of yt(2::pb%neqs)
  call compute_stress(pb%dtau_dt,dydt(3::pb%neqs),pb%kernel,v-pb%v_star)

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

  ! state evolution law, dtheta/dt = f(v,theta)
  ! SEISMIC: in Chen's formulation, theta is the gouge porosity
  dydt(1::pb%neqs) = dtheta_dt(v,yt(1::pb%neqs),pb)

  if (pb%i_rns_law == 3) then
    ! SEISMIC: the total dtau_dt results from the slip deficit and
    ! periodic loading, which is stored in dydt(2::pb%neqs)
    dydt(2::pb%neqs) = pb%dtau_dt + dtau_per
  else
    ! SEISMIC: the rate-and-state formulation computes the the time-derivative
    ! of velocity, rather than stress

    ! Time derivative of the elastic equilibrium equation
    !  dtau_load/dt + dtau_elastostatic/dt -impedance*dv/dt = sigma*( dmu/dv*dv/dt + dmu/dtheta*dtheta/dt )
    ! Rearranged in the following form:
    !  dv/dt = ( dtau_load/dt + dtau_elastostatic/dt - sigma*dmu/dtheta*dtheta/dt )/( sigma*dmu/dv + impedance )
    call dmu_dv_dtheta(dmu_dv,dmu_dtheta,yt(2::pb%neqs),yt(1::pb%neqs),pb)
    dydt(2::pb%neqs) = ( dtau_per + pb%dtau_dt - pb%sigma*dmu_dtheta*dydt(1::pb%neqs) ) &
                     / ( pb%sigma*dmu_dv + pb%zimpedance )
   endif

end subroutine derivs

end module derivs_all
