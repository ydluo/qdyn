!Module for derivs

module derivs_all

  implicit none
  
  private

  public derivs


contains


!C====================Subroutine derivs START===========================
 !-------Compute thata,vslip and time deriv for one time step------------
 !-------CALL rdft: real DFT---------------------------------------------
 !---------------------------------------------------------------|
 !-------|     in:pb                                     |-------|
 !-------|    out:pb                                     |-------|
 !-------|                                               |-------|
 !---------------------------------------------------------------|
 
subroutine derivs(pb)
   
  use problem_class
  use fftsg
  
  type(problem_type), intent(inout) :: pb

  double precision :: omega(pb%mesh%nn)
  double precision :: dtau_per
  integer :: i

  ! compute shear stress rate from elastic interactions
  !JPA: this should be replaced by "call compute_stress" which depends on dimension (0D, 1D, 2D fault)
  ! the rest of this subroutine does not depend on dimension
  if (pb%mesh%nn > 1) then
    pb%dtau_dt = pb%v_star-pb%v
    pb%dtau_dt( pb%mesh%nn+1 : pb%kernel%k2f%nnfft ) = 0d0 
    call rdft(pb%kernel%k2f%nnfft,1,pb%dtau_dt,pb%kernel%k2f%m_fft%iworkfft,pb%kernel%k2f%m_fft%rworkfft)
    pb%dtau_dt = pb%kernel%k2f%kernel * pb%dtau_dt
    call rdft(pb%kernel%k2f%nnfft,-1,pb%dtau_dt,pb%kernel%k2f%m_fft%iworkfft,pb%kernel%k2f%m_fft%rworkfft)
  else
    pb%dtau_dt(1) = pb%kernel%k2f%kernel(1)*( pb%v_star(1)-pb%v(1) )
  endif

  ! periodic loading
  dtau_per = pb%Omper * pb%Aper * dcos(pb%Omper*pb%time)     

  !--------State evolution law START--------------------------------
  omega = pb%v * pb%theta / pb%dc

  !      "aging" law
  if (pb%itheta_law == 1) then
    pb%dtheta_dt = 1.d0-omega

  !      "slip" law
  elseif (pb%itheta_law== 2) then
    do i=1,pb%mesh%nn
      if (omega > 0d0) then
        pb%dtheta_dt(i) = -omega*dlog(omega)
      else
        pb%dtheta_dt(i) = 0d0
      endif
    enddo

  !      "aging" in the no-healing approximation
  elseif (pb%itheta_law == 0) then
    pb%dtheta_dt = -omega

  endif
  !--------State evolution law END----------------------------------

  pb%dv_dt = (dtau_per+pb%dtau_dt - pb%sigma*pb%b*pb%v2*pb%dtheta_dt/(pb%v2*pb%theta+pb%dc) )    &
             /( pb%sigma*pb%a*(1.d0/pb%v-1.d0/(pb%v1+pb%v)) + pb%zimpedance )
    
end subroutine derivs



end module derivs_all
