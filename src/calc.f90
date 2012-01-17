!Calculate stress etc.

module calc

  implicit none
  private

  public compute_stress

contains
! K is stiffness, different in sign with convention in Diet92
! compute shear stress rate from elastic interactions
!   tau = - K*slip 
!   dtau_dt = - K*slip_velocity 
!
! To account for steady plate velocity, the input velocity must be v-vpl:
!   tau = - K*( slip - Vpl*t ) 
!   dtau_dt = - K*( v - Vpl )
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
    case(3); call compute_stress_3d_fft(tau,K%k3f,v)
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
!subroutine compute_stress_3d(tau,k3,v)
!
!  use problem_class, only : kernel_3D
!
!  type(kernel_3D), intent(in)  :: k3
!  double precision , intent(out) :: tau(:)
!  double precision , intent(in) :: v(:)
!
!  integer :: nn,nw,nx,i,j,jj,ix,iw,jx,jw
!
!  nn = size(v)
!  nw = size(k3%kernel,1)
!  nx = nn/nw
! ! write(6,*) 'nn,nw,nx', nn, nw, nx
!
!  tau = 0d0
!  do i = 1,nn
!    ix = mod((i-1),nx)+1          ! find column of obs
!    iw = 1+(i-ix)/nx              ! find row of obs
!    if (ix == 1)  then            ! obs at first column, directly stored in kernel (iw,nw*nx)
!      do j = 1,nn
!        tau(i) = tau(i) - k3%kernel(iw,j) * v(j)
!  !      write(6,*) i,j,k3%kernel(iw,j)
!      end do
!    else                          ! obs at other column, calculate index to get kernel
!      do j = 1,nn
!        jx = mod((j-1),nx)+1        ! find column of source
!        if (jx >= ix)  then         ! source on the right of ods, shift directly
!          tau(i) = tau(i) - k3%kernel(iw,j+1-ix) * v(j)
!   !       write(6,*) i,j,k3%kernel(iw,j+1-ix)
!        else                        ! source on the left, use symmetry
!          tau(i) = tau(i) - k3%kernel(iw,j+1+ix-2*jx) * v(j)
!   !       write(6,*) i,j,k3%kernel(iw,j+1+ix-2*jx)
!        end if
!      end do
!    end if
!
!  end do
!
!  tau = 0d0
!  i=0
!  do iw=1,nw
!  do ix=1,nx
!    i = i+1
!    j = 0
!    do jw=1,nw
!    do jx=1,nx
!      j = j+1
!      idx = abs(jx-ix)  ! note: abs(x) assumes some symmetries in the kernel
!      jj = (jw-1)*nx + idx + 1 
!      tau(i) = tau(i) - k3%kernel(iw,jj) * v(j)
!    end do
!    end do
!  end do
!  end do
!    
!end subroutine compute_stress_3d

!--------------------------------------------------------
! JPA: try this subroutine instead
subroutine compute_stress_3d_fft(tau,k3f,v)

  use problem_class, only : kernel_3D_fft
  use fftsg, only : my_rdft

  type(kernel_3D_fft), intent(inout)  :: k3f
  double precision , intent(out) :: tau(:)
  double precision , intent(in) :: v(:)

  double precision :: tmpzk(k3f%nw,k3f%nxfft), tmp(k3f%nxfft)
  integer :: n,k

! version 2, does not need to be implemented
! It requires a different storage: 
!   v and tau: (nw,nx) reshaped into vector of length nw*nx (note that nw, along-dip, is first)
!   kernel(1:nw,1:nw,0:(nx-1))
!  tau = 0d0
!  ii = 0
!  do i = 1,nx
!    mm = 0
!    do m = 1,nx
!      k = abs(i-m)
!      tau(ii+1:ii+nw) = tau(ii+1:ii+nw) + matmul( k3%kernel(:,:,k) , v(mm+1:mm+nw) )
!      mm = mm + nw
!    enddo
!    ii = ii + nw
!  enddo

!JPA version 3, with FFT along-strike version.
! Assumes kernel has been FFT'd during initialization
! Assumes storage with along-strike index running faster than along-dip
! Note: to avoid periodic wrap-around, fft convolution requires twice longer arrays
!       and zero-padding (pre-processing) and chop-in-half (post-processing)

  do n = 1,k3f%nw
    tmp( 1 : k3f%nx ) = v( (n-1)*k3f%nx+1 : n*k3f%nx )
    tmp( k3f%nx+1 : k3f%nxfft ) = 0d0  ! convolution requires zero-padding
    call my_rdft(1,tmp,k3f%m_fft) 
    tmpzk(n,:) = tmp
  enddo
  do k = 1,k3f%nxfft
    tmpzk(:,k) = matmul( k3f%kernel(:,:,k), tmpzk(:,k) )
  enddo
  do n = 1,k3f%nw
    tmp = - tmpzk(n,:)
    call my_rdft(-1,tmp,k3f%m_fft)
    tau( (n-1)*k3f%nx+1 : n*k3f%nx ) = tmp(1:k3f%nx) ! take only first half of array
  enddo

end subroutine compute_stress_3d_fft

end module calc
