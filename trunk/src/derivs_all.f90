!Module for derivs

module derivs_all

  implicit none
  
  private

  public derivs


contains


!====================Subroutine derivs===========================
 !-------Compute thata,vslip and time deriv for one time step-----
 !-------CALL rdft in compute_stress : real DFT-------------------
 !---------------------------------------------------------------|
 !-------|     in:pb, yt, dydt                           |-------|
 !-------|    out:pb, yt, dydt                           |-------|
 !-------|                                               |-------|
 !---------------------------------------------------------------|
 
subroutine derivs(time,yt,dydt,pb) 
!JPA this subroutine should have yout with intent(out)
!    and put there the derivatives
   
  use problem_class
  use calc
  
  type(problem_type), intent(inout) :: pb
  double precision , intent(inout) :: time, yt(pb%neqs*pb%mesh%nn), dydt(pb%neqs*pb%mesh%nn)

  double precision :: omega(pb%mesh%nn)
  double precision :: dtau_per
  integer :: i

  ! compute shear stress rate from elastic interactions, for 0D, 1D & 2D
  call compute_stress(pb%dtau_dt,pb%kernel,yt(2::pb%neqs),pb%vstar)

!YD This part we may want to modify it later to be able to
!impose more complicated loading/pertubation
!functions involved: problem_class/problem_type; input/read_main 
!                    initialize/init_field;  derivs_all/derivs 
  
  ! periodic loading
  dtau_per = pb%Omper * pb%Aper * dcos(pb%Omper*time)     

  !--------State evolution law START--------------------------------
  omega = yt(2::pb%neqs) *  yt(1::pb%neqs) / pb%dc

  !      "aging" law
  if (pb%itheta_law == 1) then
    dydt(1::pb%neqs) = 1.d0-omega

  !      "slip" law
  elseif (pb%itheta_law== 2) then
    do i=1,pb%mesh%nn
      if (omega(i) > 0d0) then
        dydt(i*pb%neqs-1) = -omega(i)*dlog(omega(i))
      else
        dydt(i*pb%neqs-1) = 0d0
      endif
    enddo

  !      "aging" in the no-healing approximation
  elseif (pb%itheta_law == 0) then
    dydt(1::pb%neqs) = -omega

  endif
  !--------State evolution law END----------------------------------

  dydt(2::pb%neqs) = (dtau_per+pb%dtau_dt - pb%sigma*pb%b*pb%v2*dydt(1::pb%neqs)/   &
     (pb%v2*yt(1::pb%neqs)+pb%dc) )/  &
     ( pb%sigma*pb%a*(1.d0/yt(2::pb%neqs)-1.d0/(pb%v1+yt(2::pb%neqs))) + pb%zimpedance )
    
end subroutine derivs



end module derivs_all
