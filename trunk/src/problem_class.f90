!define problem 

module problem_class
 
  use fftsg, only : OouraFFT_type

  implicit none

  public

 ! timeseries outputs: at every time step, but only macroscopic quantities
  type ot_type
    double precision :: lcold,lcnew,llocnew,llocold
    integer :: unit,ic,ntout,ivmax
  end type ot_type

 ! snapshot outputs: at every fault point, but only at few selected times
  type ox_type
    integer :: count,unit,nxout
  end type ox_type

  type mesh_type
    integer :: dim = 0  ! dim = 1, 2 ,3 ~xD
    integer :: nx, nw, nn ! along-strike, along-dip, total grid number
    double precision :: dx !along-strike grid size(constant)  
    double precision :: Lfault, W, Z_CORNER ! fault length, width, lower-left corner z (follow Okada's convention)
    double precision, allocatable :: dw(:), DIP_W(:) !along-dip grid size and dip (adjustable), nw count
    double precision, allocatable :: x(:), y(:), z(:), dip(:) !coordinates and dip of every grid (nx*nw count)
  end type mesh_type

  type kernel_2D_fft
    integer :: kind = 0
    double precision, dimension(:), allocatable :: kernel
    integer :: nnfft, finite
    type (OouraFFT_type) :: m_fft
  end type kernel_2D_fft

  type kernel_3D
    double precision, dimension(:,:), allocatable :: kernel
  end type kernel_3D

  type kernel_3D_fft
    integer :: nxfft, nw, nx
    double precision, dimension(:,:,:), allocatable :: kernel
    type (OouraFFT_type) :: m_fft
  end type kernel_3D_fft

  type kernel_type
    integer :: kind = 0
    double precision :: k1
    type (kernel_2D_fft), pointer :: k2f
    type (kernel_3D), pointer :: k3
    type (kernel_3D_fft), pointer :: k3f
  end type kernel_type

  type problem_type
    double precision, dimension(:), allocatable :: &
      tau, dtau_dt, tau_init, &
      sigma, &
      slip, v, dv_dt, theta, dtheta_dt,  &
      a, b, dc, v1, v2, mu_star, v_star, theta_star
    double precision :: pot, pot_rate
    double precision :: beta=0d0, smu=0d0, lam=0d0, zimpedance=0d0

!YD This part we may want to modify it later to be able to
!impose more complicated loading/pertubation
!functions involved: problem_class/problem_type; input/read_main 
!                    initialize/init_field;  derivs_all/derivs 

    double precision :: Tper=0d0, Aper=0d0, Omper=0d0
    double precision :: time=0d0
    integer :: itheta_law, neqs

    double precision :: dt_try, dt_did, dt_next, dt_max=0d0, tmax, acc
    integer :: NSTOP,itstop,it

    type (mesh_type) :: mesh
    type (ot_type) :: ot
    type (ox_type) :: ox
    type (kernel_type) :: kernel
  end type problem_type

end module problem_class


