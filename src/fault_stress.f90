!Calculate stress etc.

module fault_stress

  use fftsg, only : OouraFFT_type

  implicit none
  private

  type kernel_2D_fft
    double precision, dimension(:), allocatable :: kernel
    integer :: nnfft, finite
    type (OouraFFT_type) :: m_fft
  end type kernel_2D_fft

  type kernel_3D
    double precision, dimension(:,:), allocatable :: kernel, kernel_n
    integer :: nw, nx, nwLocal, nwGlobal, nnLocal, nnGlobal
  end type kernel_3D

!  type fault_coord
!    double precision, dimension(:), allocatable :: xlocarray,ylocarray,zlocarray
!    double precision, dimension(:), allocatable :: xglobarray,yglobarray,zglobarray
!  end type fault_coord

  type kernel_3D_fft
    integer :: nxfft, nw, nx, nwLocal, nwGlobal, nnLocalfft, nnGlobalfft, nnLocal, nnGlobal
    double precision, dimension(:,:,:), allocatable :: kernel, kernel_n
    type (OouraFFT_type) :: m_fft
!    type (fault_coord) :: fault
  end type kernel_3D_fft

  type kernel_3D_fft2d
    integer :: nwfft, nxfft, nw, nx
    double precision, dimension(:,:), allocatable :: kernel
    type (OouraFFT_type) :: m_fft
  end type kernel_3D_fft2d

  type kernel_type
    integer :: kind = 0
    integer :: i_sigma_cpl = 0
    double precision :: k1
    type (kernel_2D_fft), pointer :: k2f
    type (kernel_3D), pointer :: k3
    type (kernel_3D_fft), pointer :: k3f
    type (kernel_3D_fft2d), pointer :: k3f2
  end type kernel_type


! For MPI parallel
  integer, allocatable, save :: nnLocalfft_perproc(:),nnoffset_perproc(:),&
                                nnLocal_perproc(:),nnoffset_glob_perproc(:),nwLocal_perproc(:),&
                                nwoffset_glob_perproc(:)

  public :: init_kernel, compute_stress, kernel_type, nnLocalfft_perproc,nnoffset_perproc,&
            nnLocal_perproc,nnoffset_glob_perproc,nwLocal_perproc,nwoffset_glob_perproc

contains
! K is stiffness, different in sign with convention in Dieterich (1992)
! compute shear stress rate from elastic interactions
!   tau = - K*slip
!   dtau_dt = - K*slip_velocity
!
! To account for steady plate velocity, the input velocity must be v-vpl:
!   tau = - K*( slip - Vpl*t )
!   dtau_dt = - K*( v - Vpl )

!=============================================================
subroutine init_kernel(lambda,mu,m,k)

  use mesh, only : mesh_type

  double precision, intent(in) :: lambda,mu
  type(mesh_type), intent(in) :: m
  type(kernel_type), intent(inout) :: k

  write(6,*) 'Intializing kernel: ...'

  select case (k%kind)
    case(1); call init_kernel_1D(k%k1,mu,m%Lfault)
    case(2); call init_kernel_2D(k%k2f,mu,m)
    case(3); call init_kernel_3D(k%k3,lambda,mu,m,k%i_sigma_cpl==1) ! 3D no fft
    case(4); call init_kernel_3D_fft(k%k3f,lambda,mu,m,k%i_sigma_cpl==1) ! 3D with FFT along-strike
    case(5); call init_kernel_3D_fft2d(k%k3f2,lambda,mu,m) ! 3D with 2DFFT
  end select

  write(6,*) 'Kernel intialized'

end subroutine init_kernel

!----------------------------------------------------------------------
subroutine init_kernel_1D(k,mu,L)

  double precision, intent(out) :: k
  double precision, intent(in) :: mu,L

  write(6,*) 'Single degree-of-freedom spring-block system'
  k = mu/L

end subroutine init_kernel_1D


!----------------------------------------------------------------------
subroutine init_kernel_2D(k,mu,m)

  use mesh, only : mesh_type
  use constants, only : PI, SRC_PATH

  type(kernel_2d_fft), intent(inout) :: k
  type(mesh_type), intent(in) :: m
  double precision, intent(in) :: mu

  double precision, allocatable :: kk(:)
  double precision :: tau_co
  integer :: i

  k%nnfft = (k%finite+1)*m%nn
  allocate (k%kernel(k%nnfft))

  write(6,*) 'FFT applied'

  if (k%finite == 0) then

   ! Kernel for a 1D in a 2D homogeneous elastic medium,
   ! assuming antiplane deformation (slip in the direction off the 2D plane)
   !   kernel(k) = 1/2*mu*|k|
   ! where
   !   mu = shear modulus
   !   k = wavenumber
   !
   ! Kernel for 1D fault surrounded by a damaged zone of uniform thickness and compliance
   ! derived from equations 47-49 of Ampuero et al (2002, http://onlinelibrary.wiley.com/doi/10.1029/2001JB000452/full)
   ! by setting s=0 (low frequency limit = static)
   !   kernel(k) = 1/2*mu*(1-D)*|k| * cotanh[ H*|k| + arctanh(1-D) ]
   ! where
   !   H = half-thickness of the fault damage zone
   !   D = damage level = 1 - (damaged shear modulus) / (intact shear modulus)

   ! Define wavenumber
    allocate( kk(k%nnfft) )
   !kk(1:2:m%nn-1) = 2d0*PI/m%Lfault*(/ i, i=0:m%nn/2-1 /)
    do i=0,m%nn/2-1
      kk(2*i+1) = 2d0*PI/m%Lfault*i
    enddo
    kk(2:2:m%nn) = kk(1:2:m%nn-1)
    kk(2) = PI*m%nn/m%Lfault ! Nyquist

   ! To mimic the width W of the fault in the dimension normal to the 2D plane,
   ! we introduce a "2.5D" approximation: we replace k by sqrt(k^2 + (2*pi/W)^2)
    kk = sqrt( kk*kk + (2*PI/m%W)**2 )

   ! Compute kernel for homogeneous medium
    k%kernel = 0.5d0*mu*kk

   ! Compute kernel for damaged medium
   ! TO DO : define D and H as inputs, then uncomment the lines below
   ! if (D>0 .and. H>0) k%kernel = k%kernel * (1-D) * cotanh( H*kk + arctanh(1-D) )

   ! factor 2/N from the inverse FFT convention
    k%kernel = k%kernel *2.d0/m%nn

 !- Read coefficient I(n) from pre-calculated file.
  elseif (k%finite == 1) then
    open(57,file=SRC_PATH//'/kernel_I.tab')
    do i=1,k%nnfft/2-1
      read(57,*,end=100) k%kernel(2*i+1)
    enddo
    read(57,*,end=100) k%kernel(2) ! Nyquist
    close(57)
    ! The factor 2/N comes from the inverse FFT convention
    tau_co = PI*mu / (2d0*m%Lfault) *2.d0/k%nnfft
    k%kernel(1) = 0d0
    k%kernel(2) = tau_co*dble(k%nnfft/2)*k%kernel(2)
    do i = 1,k%nnfft/2-1
      k%kernel(2*i+1) = tau_co*dble(i)*k%kernel(2*i+1)
      k%kernel(2*i+2) = k%kernel(2*i+1)
    enddo
  end if

  return

100 stop 'Kernel file src/kernel_I.tab is too short. Use src/TabKernelFiniteFlt.m to create a longer one.'

end subroutine init_kernel_2D

!----------------------------------------------------------------------
subroutine init_kernel_3D_fft(k,lambda,mu,m,sigma_coupling)

  use mesh, only : mesh_type
  use okada, only : compute_kernel
  use fftsg, only : my_rdft
  use utils, only : save_vector3
  use constants, only : FAULT_TYPE
  use my_mpi, only : is_mpi_parallel, is_mpi_master, my_mpi_rank

  type(kernel_3d_fft), intent(inout) :: k
  double precision, intent(in) :: lambda,mu
  type(mesh_type), intent(in) :: m
  logical, intent(in) :: sigma_coupling

  double precision :: tau,sigma_n, y_src, z_src, dip_src, dw_src, y_obs, z_obs,dip_obs
  double precision, allocatable :: tmp(:), tmp_n(:)   ! for FFT
  integer :: i, j, ii, jj, n, nn, IRET

  if (is_mpi_master()) then
    write(6,*) 'Generating 3D kernel...'
    write(6,*) 'OouraFFT applied along-strike'
  endif

  k%nx = m%nx
  k%nxfft = 2*m%nx ! fft convolution requires twice longer array

  if (is_MPI_parallel()) then
    write(6,*) 'MY_RANK:', my_mpi_rank()
    k%nwLocal = m%nw
    k%nnLocal=nnoffset_glob_perproc(my_mpi_rank())
    k%nwGlobal=sum(nwLocal_perproc)
    k%nnGlobal=k%nwGlobal*k%nx

  else
    k%nwGlobal = m%nw
    k%nnGlobal=k%nwGlobal*k%nx
    k%nwLocal = k%nwGlobal
    k%nnLocal=0
  endif

  write(6,*) 'nwGlobal:',k%nwGlobal
  write(6,*) 'nwLocal:',k%nwLocal
  write(6,*) 'nnLocal:',k%nnLocal
  k%nnLocalfft  = k%nwLocal*k%nxfft
  k%nnGlobalfft = k%nwGlobal*k%nxfft

  allocate(k%kernel(k%nwLocal,k%nwGlobal,k%nxfft))
  allocate(tmp(k%nxfft))
  if (sigma_coupling) allocate(k%kernel_n(k%nwLocal,k%nwGlobal,k%nxfft))
  allocate(tmp_n(k%nxfft))
  ! assumes faster index runs along-strike
  do n=1,k%nwGlobal
    nn = (n-1)*m%nx+1
    ! note that if serial (no MPI) the glob arrays point to the local arrays
    y_src = m%yglob(nn)
    z_src = m%zglob(nn)
    dip_src = m%dipglob(nn)
    dw_src = m%dwglob(n)
    do j=1,k%nwLocal
      jj = (j-1)*m%nx+1
      y_obs = m%y(jj)
      z_obs = m%z(jj)
      dip_obs = m%dip(jj)
      do i=-m%nx+1,m%nx
        call compute_kernel(lambda,mu, &
                i*m%dx, y_src, z_src, dip_src, m%dx, dw_src,   &
                0d0, y_obs, z_obs, dip_obs, &
                IRET,tau,sigma_n,FAULT_TYPE)
        ii = i+1
        ! wrap up the negative relative-x-positions in the second half of the
        ! array to comply with conventions of fft convolution
        if (i<0) ii = ii + k%nxfft
        tmp(ii) = tau
        tmp_n(ii) = sigma_n
      enddo
      call my_rdft(1,tmp,k%m_fft)
      k%kernel(j,n,:) = tmp / dble(m%nx)
      if (sigma_coupling) then
        call my_rdft(1,tmp_n,k%m_fft)
        k%kernel_n(j,n,:) = tmp_n / dble(m%nx)
      endif
    enddo
  enddo

end subroutine init_kernel_3D_fft

! -----------------------------------------------

subroutine init_kernel_3D_fft2d(k,lambda,mu,m)

  use mesh, only : mesh_type
  use fftsg, only : my_rdft2

  type(kernel_3d_fft2d), intent(inout) :: k
  double precision, intent(in) :: lambda, mu
  type(mesh_type), intent(in) :: m

  double precision, allocatable, dimension(:,:) :: Kij
  double precision :: Im, Ip, Jm, Jp, ImJp, ImJm, IpJp, IpJm, T1, T2, T3, T4, koef
  double precision, parameter :: PI = 3.14159265358979d0
  integer :: i, j

  write(6,*) 'Generating 3D kernel...'
  write(6,*) 'Ooura FFT2 applied for fault plane'
  k%nw = m%nw
  k%nx = m%nx
  k%nwfft = 2 * m%nw ! fft convolution requires twice longer array
  k%nxfft = 2 * m%nx
  allocate(k%kernel(k%nxfft,k%nwfft))
  allocate(Kij(k%nxfft,k%nwfft))
  write(6,*) 'Allocated for FFT2 Dimensions of ', k%nxfft, ' x ', k%nwfft

  koef = 2.0d0 * (lambda + mu) / (lambda + 2.0d0*mu)

  ! Construct the static kernel (modified from F. Gallovic code); note that slip
  ! is assumed to be in the along-strike direction.
  ! We mirror the kernel along-strike and along-dip for 2D convolution
  do i = 1,k%nxfft
    if (i <= m%nx) then
      Im = (dble(i-1) - 0.5d0) * m%dx
      Ip = (dble(i-1) + 0.5d0) * m%dx
    else
      Im = (dble(m%nx*2 - i + 1) - .5d0) * m%dx
      Ip = (dble(m%nx*2 - i + 1) + .5d0) * m%dx
    endif
    do j = 1,k%nwfft
      if (j <= m%nw) then
        Jm = (dble(j-1) - .5d0) * m%dw(1)
        Jp = (dble(j-1) + .5d0) * m%dw(1)
      else
        Jm = (dble(m%nw*2 - j + 1) - .5d0) * m%dw(1)
        Jp = (dble(m%nw*2 - j + 1) + .5d0) * m%dw(1)
      endif
      ImJp = sqrt(Im**2 + Jp**2)
      ImJm = sqrt(Im**2 + Jm**2)
      IpJp = sqrt(Ip**2 + Jp**2)
      IpJm = sqrt(Ip**2 + Jm**2)
      T1 = (Jp/ImJp - Jm/ImJm) / Im
      T2 = (Jp/IpJp - Jm/IpJm) / Ip
      T3 = (Ip/IpJm - Im/ImJm) / Jm
      T4 = (Ip/IpJp - Im/ImJp) / Jp
      Kij(i,j) = mu/(4.d0*PI) * ( koef*(T1 - T2) + T3 - T4 )
    enddo
  enddo

  ! Perform 2D FFT on kernel and pre-normalize
  call my_rdft2(1, Kij, k%m_fft)
  k%kernel = Kij * 2.0d0 / dble(k%nwfft * k%nxfft)

  deallocate(Kij)

end subroutine init_kernel_3D_fft2d

!----------------------------------------------------------------------
subroutine init_kernel_3D(k,lambda,mu,m,sigma_coupling)

  use mesh, only : mesh_type
  use okada, only : compute_kernel
  use constants, only : FAULT_TYPE

  type(kernel_3d), intent(inout) :: k
  double precision, intent(in) :: lambda,mu
  type(mesh_type), intent(in) :: m
  logical, intent(in) :: sigma_coupling

  double precision :: tau, sigma_n
  integer :: i, j, IRET

    write(6,*) 'Generating 3D kernel...'
    write(6,*) 'NO FFT applied'
    !kernel(i,j): response at i of source at j
    !because dx = constant, only need to calculate i at first column

    allocate (k%kernel(m%nw,m%nn))
    if (sigma_coupling) allocate (k%kernel_n(m%nw,m%nn))
    do i = 1,m%nw
      do j = 1,m%nn
        call compute_kernel(lambda,mu,m%x(j),m%y(j),m%z(j),  &
               m%dip(j),m%dx,m%dw((j-1)/m%nx+1),   &
               m%x(1+(i-1)*m%nx),m%y(1+(i-1)*m%nx),   &
               m%z(1+(i-1)*m%nx),m%dip(1+(i-1)*m%nx),IRET,tau,sigma_n,FAULT_TYPE)
        if (IRET == 0) then
          k%kernel(i,j) = tau
          if (sigma_coupling) k%kernel_n(i,j) = sigma_n
        else
          write(6,*) '!!WARNING!! : Kernel Singular, set value to 0,(i,j)',i,j
          k%kernel(i,j) = 0d0
          if (sigma_coupling) k%kernel_n(i,j) = 0d0
        end if
      end do
    end do
    do j = 1,m%nn
      write(99,*) k%kernel(1,j)
      if (sigma_coupling) write(99,*) k%kernel_n(1,j)
    end do

end subroutine init_kernel_3D

!=========================================================
subroutine compute_stress(tau,sigma_n,K,v)

  type(kernel_type), intent(inout)  :: K
  double precision , intent(out) :: tau(:), sigma_n(:)
  double precision , intent(in) :: v(:)

  ! depends on dimension (0D, 1D or 2D fault)
  select case (K%kind)
    case(1); call compute_stress_1d(tau,K%k1,v)
    case(2); call compute_stress_2d(tau,K%k2f,v)
    case(3)
      call compute_stress_3d(tau,sigma_n,K%k3,v)
    case(4); call compute_stress_3d_fft(tau,sigma_n,K%k3f,v)
    case(5); call compute_stress_3d_fft2d(tau,K%k3f2,v)
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

subroutine compute_stress_3d(tau,sigma_n,k3,v)

  type(kernel_3D), intent(in)  :: k3
  double precision, intent(inout) :: tau(:), sigma_n(:)
  double precision, intent(in) :: v(:)

  integer :: nnLocal,nnGlobal,nw,nxLocal,nxGlobal,k,iw,ix,j,jw,jx,idx,jj,ix0_proc
  double precision :: tsum

  nnLocal = size(tau)
  nw = size(k3%kernel,1)
  nxLocal = nnLocal/nw

  nnGlobal = size(v)
  nxGlobal = nnGlobal/nw

  ix0_proc = 0

  !$OMP PARALLEL PRIVATE(iw,ix,tsum,idx,jw,jx,jj,j,k)
  !$OMP DO SCHEDULE(STATIC)
   do k=1,nnLocal
     iw = (k-1)/nxLocal +1
     ix = k-(iw-1)*nxLocal + ix0_proc
       j = 0
       tsum = 0.0d0
       do jw=1,nw
         do jx=1,nxGlobal
           j = j+1
           idx = abs(jx-ix)  ! note: abs(x) assumes some symmetries in the kernel
           jj = (jw-1)*nxGlobal + idx + 1
           tsum = tsum - k3%kernel(iw,jj) * v(j)
           !NOTE: it could be more efficient to store kernel in (jj,iw) form
         end do
       end do
     tau(k) = tsum
   end do
  !$OMP END DO
  !$OMP END PARALLEL


  if (allocated(k3%kernel_n)) then

  !$OMP PARALLEL PRIVATE(iw,ix,tsum,idx,jw,jx,jj,j,k)
  !$OMP DO SCHEDULE(STATIC)
   do k=1,nnLocal
     iw = (k-1)/nxLocal +1
     ix = k-(iw-1)*nxLocal + ix0_proc
       j = 0
       tsum = 0.0d0
       do jw=1,nw
         do jx=1,nxGlobal
           j = j+1
           idx = abs(jx-ix)  ! note: abs(x) assumes some symmetries in the kernel
           jj = (jw-1)*nxGlobal + idx + 1
           tsum = tsum - k3%kernel_n(iw,jj) * v(j)
         end do
       end do
       sigma_n(k) = tsum
   end do
  !$OMP END DO
  !$OMP END PARALLEL

  end if

end subroutine compute_stress_3d

!--------------------------------------------------------
! version with FFT along-strike
! Assumes kernel has been FFT'd during initialization
! Assumes storage with along-strike index running faster than along-dip
! Note: to avoid periodic wrap-around, fft convolution requires twice longer arrays,
!       zero-padding (pre-processing) and chop-in-half (post-processing)
!
! MPI partitioning: each processor gets a range of depths

subroutine compute_stress_3d_fft(tau,sigma_n,k3f,v)

  use fftsg, only : my_rdft
  use utils, only : save_vector
  use my_mpi, only : is_MPI_parallel, gather_allvdouble

  type(kernel_3D_fft), intent(inout)  :: k3f
  double precision , intent(inout) :: tau(:), sigma_n(:) !PG: Collect tau and sigma_n in all processor.
  double precision , intent(in) :: v(:)
  double precision :: vzk(k3f%nwGlobal,k3f%nxfft), tmpzk(k3f%nwLocal,k3f%nxfft)
  double precision :: tmpx(k3f%nxfft)
  integer :: n,k
  integer :: iglobal,ilocal,iwlocal,iwglobal,ixlocal,ixglobal
  double precision :: tmpzkarray(k3f%nnLocalfft),vzkarray(k3f%nnGlobalfft)

!tmpx needs to be private to avoid superposition with the other threads.
!$OMP PARALLEL PRIVATE(tmpx)

!-- load velocity and apply FFT along strike

!$OMP DO SCHEDULE(STATIC)
  do n = 1,k3f%nwLocal
    tmpx( 1 : k3f%nx ) = v( (n-1)*k3f%nx+1 : n*k3f%nx )
    tmpx( k3f%nx+1 : k3f%nxfft ) = 0d0  ! convolution requires zero-padding
    call my_rdft(1,tmpx,k3f%m_fft)
    tmpzk(n,:) = tmpx
  enddo
!$OMP END DO

!$OMP SINGLE
  if (is_MPI_parallel()) then
   !Local array, convert from matrix to vector form
    ilocal=0
    do iwlocal=1,k3f%nwLocal
      do ixlocal=1,k3f%nxfft
        ilocal = ilocal+1
        tmpzkarray(ilocal) = tmpzk(iwlocal,ixlocal)
      enddo
    enddo

   ! gather the global vzk from the pieces in all processors
    call gather_allvdouble(tmpzkarray,k3f%nnLocalfft,vzkarray,nnLocalfft_perproc, &
                     nnoffset_perproc,k3f%nnGlobalfft)

   !Global array, convert from vector to matrix form
    iglobal=0
    do iwglobal=1,k3f%nwGlobal
      do ixglobal=1,k3f%nxfft
            iglobal=iglobal+1
         vzk(iwglobal,ixglobal)=vzkarray(iglobal)
      enddo
    enddo

  else
    vzk = tmpzk
  endif
!$OMP END SINGLE

  !-- compute shear stress

  ! convolution in Fourier domain is a product of complex numbers:
  ! K*V = (ReK + i*ImK)*(ReV+i*ImV)
  !     = ReK*ReV - ImK*ImV  + i*( ReK*ImV + ImK*ReV )
  !
  !$OMP SINGLE
  ! wavenumber = 0, real
  tmpzk(:,1) = matmul( k3f%kernel(:,:,1), vzk(:,1) )
  ! wavenumber = Nyquist, real
  tmpzk(:,2) = matmul( k3f%kernel(:,:,2), vzk(:,2) )
  !$OMP END SINGLE

  ! higher wavenumbers, complex
  !$OMP DO SCHEDULE(STATIC)
  do k = 3,k3f%nxfft-1,2
    ! real part = ReK*ReV - ImK*ImV
    tmpzk(:,k)   = matmul( k3f%kernel(:,:,k), vzk(:,k) )  &
                 - matmul( k3f%kernel(:,:,k+1), vzk(:,k+1) )
    ! imaginary part = ReK*ImV + ImK*ReV
    tmpzk(:,k+1) = matmul( k3f%kernel(:,:,k), vzk(:,k+1) )&
                 + matmul( k3f%kernel(:,:,k+1), vzk(:,k) )
  enddo
  !$OMP END DO

  !$OMP DO SCHEDULE(STATIC)
  do n = 1,k3f%nwLocal
    tmpx = - tmpzk(n,:)
    call my_rdft(-1,tmpx,k3f%m_fft)
    tau( (n-1)*k3f%nx+1 : n*k3f%nx ) = tmpx(1:k3f%nx) ! take only first half of array
  enddo
  !$OMP END DO

!-- compute normal stress

  if (allocated(k3f%kernel_n)) then
  ! Same steps as for shear stress

    !$OMP SINGLE
    tmpzk(:,1) = matmul( k3f%kernel_n(:,:,1), vzk(:,1) )
    tmpzk(:,2) = matmul( k3f%kernel_n(:,:,2), vzk(:,2) )
    !$OMP END SINGLE
    !$OMP DO SCHEDULE(STATIC)
    do k = 3,k3f%nxfft-1,2
      tmpzk(:,k)   = matmul( k3f%kernel_n(:,:,k), vzk(:,k) )  &
                   - matmul( k3f%kernel_n(:,:,k+1), vzk(:,k+1) )
      tmpzk(:,k+1) = matmul( k3f%kernel_n(:,:,k), vzk(:,k+1) )  &
                   + matmul( k3f%kernel_n(:,:,k+1), vzk(:,k) )
    enddo
    !$OMP END DO

    !$OMP DO SCHEDULE(STATIC)
    do n = 1,k3f%nwLocal
      tmpx = - tmpzk(n,:)
      call my_rdft(-1,tmpx,k3f%m_fft)
      sigma_n( (n-1)*k3f%nx+1 : n*k3f%nx ) = tmpx(1:k3f%nx)
    enddo
    !$OMP END DO

  endif

!$OMP END PARALLEL

end subroutine compute_stress_3d_fft

! ---------------------------------------------------------------------------
! Version to perform the 2D-convolution using 2D-FFTs of the velocity on the
! fault plane and a pre-computed 2D-FFT of an elastic kernel in an infinite medium.
! The kernel assumes slip is in the along-strike direction.

subroutine compute_stress_3d_fft2d(tau, k, v)

  use fftsg, only : my_rdft2

  double precision, intent(out) :: tau(:)
  type(kernel_3d_fft2d), intent(inout) :: k
  double precision, intent(in) :: v(:)

  double precision :: tmpx(k%nxfft, k%nwfft), tmpz(k%nxfft, k%nwfft)
  integer :: nw, nx, nwfft, nxfft, nw2, nw21

  ! Retrieve dimensions and half indices
  nw = k%nw
  nx = k%nx
  nwfft = k%nwfft
  nxfft = k%nxfft
  nw2 = nwfft / 2 + 1
  nw21 = nw2 + 1

  ! Store the velocity and zero-pad (faster index along-strike)
  tmpx = 0.0d0
  tmpx(1:nx,1:nw) = reshape(v, (/ nx, nw /))

  ! Compute the 2D-FFT on the velocity
  call my_rdft2(1, tmpx, k%m_fft)
  tmpz = tmpx

  ! Pointwise multiply with the FFT'd kernel
  ! Complex multiplication: do real parts first
  tmpx(1:2,1)   = k%kernel(1:2,1)   * tmpz(1:2,1)
  tmpx(1:2,nw2) = k%kernel(1:2,nw2) * tmpz(1:2,nw2)
  ! Now do higher wavenumbers
  tmpx(3::2,:) = k%kernel(3::2,:) * tmpz(3::2,:) &
               - k%kernel(4::2,:) * tmpz(4::2,:)
  tmpx(4::2,:) = k%kernel(3::2,:) * tmpz(4::2,:) &
               + k%kernel(4::2,:) * tmpz(3::2,:)
  tmpx(1,2:nw) = k%kernel(1,2:nw) * tmpz(1,2:nw) &
               - k%kernel(2,2:nw) * tmpz(2,2:nw)
  tmpx(2,2:nw) = k%kernel(1,2:nw) * tmpz(2,2:nw) &
               + k%kernel(2,2:nw) * tmpz(1,2:nw)
  ! Upper right is switched: first row is imaginary, second row is real
  tmpx(2,nw21:) = k%kernel(2,nw21:) * tmpz(2,nw21:) &
                - k%kernel(1,nw21:) * tmpz(1,nw21:)
  tmpx(1,nw21:) = k%kernel(2,nw21:) * tmpz(1,nw21:) &
                + k%kernel(1,nw21:) * tmpz(2,nw21:)

  ! Compute inverse 2D-FFT on complex product
  call my_rdft2(-1, tmpx, k%m_fft)

  ! Extract only the valid part
  tau = reshape(tmpx(1:nx,1:nw), (/ nx*nw /))

end subroutine compute_stress_3d_fft2d

end module fault_stress
