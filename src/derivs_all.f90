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
  use friction, only : dtheta_dt, dmu_dv_dtheta

  type(problem_type), intent(inout) :: pb
  double precision, intent(in) :: time, yt(pb%neqs*pb%mesh%nn)
  double precision, intent(out) :: dydt(pb%neqs*pb%mesh%nn)

  double precision, dimension(pb%mesh%nn) :: dmu_dv, dmu_dtheta
  double precision :: dtau_per

  ! storage conventions:
  ! v = yt(2::pb%neqs)
  ! theta = yt(1::pb%neqs)
  ! dv/dt = dydt(2::pb%neqs)
  ! dtheta/dt = dydt(1::pb%neqs)

  ! compute shear stress rate from elastic interactions, for 0D, 1D & 2D
  call compute_stress(pb%dtau_dt,dydt(3::pb%neqs),pb%kernel,yt(2::pb%neqs)-pb%v_star)

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
  dydt(1::pb%neqs) = dtheta_dt(yt(2::pb%neqs),yt(1::pb%neqs),pb)

  ! Time derivative of the elastic equilibrium equation
  !  dtau_load/dt + dtau_elastostatic/dt -impedance*dv/dt = sigma*( dmu/dv*dv/dt + dmu/dtheta*dtheta/dt )
  ! Rearranged in the following form:
  !  dv/dt = ( dtau_load/dt + dtau_elastostatic/dt - sigma*dmu/dtheta*dtheta/dt )/( sigma*dmu/dv + impedance )

  ! SEISMIC: the acceleration of slip (dv/dt = dydt(2::pb%neqs)) has been added to the input parameters
  ! TODO: how does this affect integration? dv/dt is updated again after this call
  call dmu_dv_dtheta(dmu_dv,dmu_dtheta,yt(2::pb%neqs),dydt(2::pb%neqs),yt(1::pb%neqs),pb)
  dydt(2::pb%neqs) = ( dtau_per + pb%dtau_dt - pb%sigma*dmu_dtheta*dydt(1::pb%neqs) ) &
                   / ( pb%sigma*dmu_dv + pb%zimpedance )

end subroutine derivs

end module derivs_all
