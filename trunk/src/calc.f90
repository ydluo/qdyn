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

  use problem_class, only : kernel_type

  type(kernel_type), intent(inout)  :: K
  double precision , intent(out) :: dtau_dt(:)
  double precision , intent(in) :: v(:),v_pl(:)

  select case (K%kind)
    case(1); call compute_stress_1d(dtau_dt,K%k1,v,v_pl)
    case(2); call compute_stress_2d(dtau_dt,K%k2f,v,v_pl)
    case(3); call compute_stress_3d(dtau_dt,K%k3,v,v_pl)
  end select

end subroutine compute_stress

!--------------------------------------------------------
subroutine compute_stress_1d(dtau_dt,k1,v,v_pl)

  double precision , intent(out) :: dtau_dt(1)
  double precision , intent(in) :: k1,v(1),v_pl(1)

  dtau_dt = k1*( v_pl-v )

end subroutine compute_stress_1d

!--------------------------------------------------------
subroutine compute_stress_2d(dtau_dt,k2f,v,v_pl)

  use problem_class, only : kernel_2d_fft
  use fftsg, only : my_rdft
  
  type(kernel_2d_fft), intent(inout)  :: k2f
  double precision , intent(out) :: dtau_dt(:)
  double precision , intent(in) :: v(:),v_pl(:)

  double precision :: tmp(k2f%nnfft)
  integer :: nn

  nn = size(v)
  tmp( 1 : nn ) = v_pl - v
  tmp( nn+1 : k2f%nnfft ) = 0d0  
  call my_rdft(1,tmp,k2f%m_fft) 
  tmp = k2f%kernel * tmp
  call my_rdft(-1,tmp,k2f%m_fft)
  dtau_dt = tmp(1:nn)

end subroutine compute_stress_2d

end module calc
