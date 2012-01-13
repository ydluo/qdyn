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

!--------------------------------------------------------
subroutine compute_stress_3d(dtau_dt,k3,v,v_pl)

  use problem_class, only : kernel_3D
  type(kernel_3D), intent(in)  :: k3
  double precision , intent(out) :: dtau_dt(1)
  double precision , intent(in) :: v(:),v_pl(:)
  integer :: nn,nw,nx,i,j,ix,iw,jx  

  nn = size(v)
  nw = size(k3%kernel)/nn
  nx = nn/nw
  
  dtau_dt = 0d0
  do i = 1,nn
    ix = mod((i-1),nx)+1          ! find column of obs
    iw = 1+(i-ix)/nx              ! find row of obs
    if (ix == 1)  then            ! obs at first column, directly stored in kernel (iw,nw*nx)
      do j = 1,nn
        dtau_dt(i) = dtau_dt(i) + k3%kernel(iw,j)*( v_pl(j)-v(j) )
      end do
    else                          ! obs at other column, calculate index to get kernel
      do j = 1,nn
        jx = mod((j-1),nx)+1        ! find column of source
        if (jx >= ix)  then         ! source on the right of ods, shift directly
          dtau_dt(i) = dtau_dt(i) + k3%kernel(iw,j+1-ix)*( v_pl(j)-v(j) )
        else                        ! source on the left, use symmetricity
          dtau_dt(i) = dtau_dt(i) + k3%kernel(iw,j+1+ix-2*jx)*( v_pl(j)-v(j) )
        end if
      end do
    end if
  end do

end subroutine compute_stress_3d

end module calc
