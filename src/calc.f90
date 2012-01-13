!Calculate stress etc.

module calc

  implicit none
  private

  public compute_stress

contains

! compute shear stress rate from elastic interactions
!   tau = - K*slip 
!   dtau_dt = - K*slip_velocity 
!
! To account for steady plate velocity, the input velocity must be v-vpl:
!   tau = K*( slip - Vpl*t ) 
!   dtau_dt =  K*( vpl - v )
!
! compute_stress depends on dimension (0D, 1D, 2D fault)
!
subroutine compute_stress(tau,K,v)

  use problem_class, only : kernel_type

  type(kernel_type), intent(inout)  :: K
  double precision , intent(out) :: tau(:)
  double precision , intent(in) :: v(:)

  select case (K%kind)
    case(1); call compute_stress_1d(tau,K%k1,v)
    case(2); call compute_stress_2d(tau,K%k2f,v)
    case(3); call compute_stress_3d(tau,K%k3,v)
  end select

end subroutine compute_stress

!--------------------------------------------------------
subroutine compute_stress_1d(tau,k1,v)

  double precision , intent(out) :: tau(1)
  double precision , intent(in) :: k1,v(1)

  tau =  - k1*v

end subroutine compute_stress_1d

!--------------------------------------------------------
subroutine compute_stress_2d(tau,k2f,v)

  use problem_class, only : kernel_2d_fft
  use fftsg, only : my_rdft
  
  type(kernel_2d_fft), intent(inout)  :: k2f
  double precision , intent(out) :: tau(:)
  double precision , intent(in) :: v(:)

  double precision :: tmp(k2f%nnfft)
  integer :: nn

  nn = size(v)
  tmp( 1 : nn ) = v
  tmp( nn+1 : k2f%nnfft ) = 0d0  
  call my_rdft(1,tmp,k2f%m_fft) 
  tmp = - k2f%kernel * tmp
  call my_rdft(-1,tmp,k2f%m_fft)
  tau = tmp(1:nn)

end subroutine compute_stress_2d

!--------------------------------------------------------
subroutine compute_stress_3d(tau,k3,v)

  use problem_class, only : kernel_3D
  type(kernel_3D), intent(in)  :: k3
  double precision , intent(out) :: tau(1)
  double precision , intent(in) :: v(:)
  integer :: nn,nw,nx,i,j,ix,iw,jx  

  nn = size(v)
  nw = size(k3%kernel,1)
  nx = nn/nw
  
  tau = 0d0
  do i = 1,nn
    ix = mod((i-1),nx)+1          ! find column of obs
    iw = 1+(i-ix)/nx              ! find row of obs
    if (ix == 1)  then            ! obs at first column, directly stored in kernel (iw,nw*nx)
      do j = 1,nn
        tau(i) = tau(i) - k3%kernel(iw,j) * v(j)
      end do
    else                          ! obs at other column, calculate index to get kernel
      do j = 1,nn
        jx = mod((j-1),nx)+1        ! find column of source
        if (jx >= ix)  then         ! source on the right of ods, shift directly
          tau(i) = tau(i) - k3%kernel(iw,j+1-ix) * v(j)
        else                        ! source on the left, use symmetry
          tau(i) = tau(i) - k3%kernel(iw,j+1+ix-2*jx) * v(j)
        end if
      end do
    end if
  end do

end subroutine compute_stress_3d

end module calc