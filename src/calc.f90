!Calculate stress etc.

module calc

  implicit none
  private

  public compute_stress

contains

subroutine compute_stress(pb)

  use fftsg, only : my_rdft
  use problem_class
  type(problem_type), intent(inout)  :: pb

  ! compute shear stress rate from elastic interactions
  ! compute_stress depends on dimension (0D, 1D, 2D fault)

  if (pb%mesh%kind == 0) then  ! 0D or 1D fault 
    if (pb%mesh%nn > 1) then   ! 1D fault
      pb%dtau_dt = pb%v_star-pb%v
      pb%dtau_dt( pb%mesh%nn+1 : pb%kernel%k2f%nnfft ) = 0d0 
      call my_rdft(1,pb%dtau_dt,pb%kernel%k2f%m_fft) 
      pb%dtau_dt = pb%kernel%k2f%kernel * pb%dtau_dt
      call my_rdft(-1,pb%dtau_dt,pb%kernel%k2f%m_fft)
    else
      pb%dtau_dt(1) = pb%kernel%k2f%kernel(1)*( pb%v_star(1)-pb%v(1) )   !0D fault
    endif
  endif
  
end subroutine compute_stress

end module calc
