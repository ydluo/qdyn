!Calculate stress etc.

module calc

  implicit none
  private

  public compute_stress

contains

! compute shear stress rate from elastic interactions
! compute_stress depends on dimension (0D, 1D, 2D fault)
!
subroutine compute_stress(dtau_dt,K,v,v_pl)

  use fftsg, only : my_rdft
  use problem_class, only :: kernel_type

  type(kernel_type), intent(inout)  :: K
  double precision , intent(out) :: dtau_dt(:)
  double precision , intent(in) :: v(:)

  double precision :: tmp(K%k2f%nnfft)
  integer :: nn

  nn = size(v)

  if (K%kind == 0) then  ! 0D or 1D fault 
    if (nn > 1) then   ! 1D fault
      tmp( 1 : nn ) = v_pl - v
      tmp( nn+1 : K%k2f%nnfft ) = 0d0  
      call my_rdft(1,tmp,K%k2f%m_fft) 
      tmp = K%k2f%kernel * tmp
      call my_rdft(-1,tmp,K%k2f%m_fft)
      dtau_dt = tmp(1:nn)
    else
      dtau_dt(1) = K%k2f%kernel(1)*( v_pl(1)-v(2) )   !0D fault
    endif
  endif
  
end subroutine compute_stress

end module calc
