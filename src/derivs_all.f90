!Module for derivs

module derivs_all

  implicit none
  
  private

  public :: derivs


contains


!====================Subroutine derivs===========================
 !-------Compute thata,vslip and time deriv for one time step-----
 !-------CALL rdft in compute_stress : real DFT-------------------
 !---------------------------------------------------------------|
 
subroutine derivs(time,yt,dydt,pb) 
   
  use problem_class
  use fault_stress, only : compute_stress
  
  type(problem_type), intent(inout) :: pb
  double precision , intent(in) :: time, yt(pb%neqs*pb%mesh%nn)
  double precision , intent(out) :: dydt(pb%neqs*pb%mesh%nn)

  double precision :: omega(pb%mesh%nn)
  double precision :: dtau_per

  ! compute shear stress rate from elastic interactions, for 0D, 1D & 2D
  call compute_stress(pb%dtau_dt,pb%kernel,yt(2::pb%neqs)-pb%v_star)

!YD we may want to modify this part later to be able to
!impose more complicated loading/pertubation
!functions involved: problem_class/problem_type; input/read_main 
!                    initialize/init_field;  derivs_all/derivs 
  ! periodic loading
  dtau_per = pb%Omper * pb%Aper * cos(pb%Omper*time)     

  ! state evolution law
  omega = yt(2::pb%neqs) *  yt(1::pb%neqs) / pb%dc
  select case (pb%itheta_law)
  case(1) ! "aging" law
    dydt(1::pb%neqs) = 1.d0-omega
  case(2) ! "slip" law
    dydt(1::pb%neqs) = -omega*log(omega)
  case(0) ! "aging" in the no-healing approximation
    dydt(1::pb%neqs) = -omega
  end select

  dydt(2::pb%neqs) = (dtau_per+pb%dtau_dt - pb%sigma*pb%b*pb%v2*dydt(1::pb%neqs)/   &
     (pb%v2*yt(1::pb%neqs)+pb%dc) )/  &
     ( pb%sigma*pb%a*(1.d0/yt(2::pb%neqs)-1.d0/(pb%v1+yt(2::pb%neqs))) + pb%zimpedance )
    
end subroutine derivs

end module derivs_all
