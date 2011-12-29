!define problem 

 module problem_class
 
  implicit none

  public


  type output_type
    integer :: ntout,nxout
    double precision :: lcold,lcnew,llocnew,llocold
  end type output_type

  type mesh_type
    integer :: kind = 0
    integer :: nn;
    double precision :: dx, Lfault, W 
  end type mesh_type

  type OouraFFT_type
    !     working arrays for Ooura's fft
    integer :: nwfft
    integer,allocatable :: iworkfft(:)
    double precision,allocatable :: rworkfft(:)
  end type OouraFFT_type

  type kernel_2D_fft
    integer :: kind = 0
    double precision, dimension(:), allocatable :: kernel
    integer :: nnfft, finite
    type (OouraFFT_type) :: m_fft
  end type kernel_2D_fft

  type kernel_type
    integer :: kind = 0
    type (kernel_2D_fft) :: k2f
  end type kernel_type




  type constant_type
    double precision :: day, week, month, year, pi   
  end type constant_type


  type problem_type
    double precision, dimension(:), allocatable :: &
      tau, dtau_dt, tau_init, &
      sigma, &
      slip, v, dv_dt, theta, dtheta_dt,  &
      a, b, dc, v1, v2, mu_star, v_star, theta_star
    double precision :: beta, smu, zimpedance
    double precision :: Tper, Aper, Omper
    double precision :: time
    integer :: itheta_law, neqs

    double precision :: dt_try, dt_did, dt_next, dt_max, tmax, acc
    integer :: NSTOP,itstop

    type (mesh_type) :: mesh
    type (output_type) :: output
    type (kernel_type) :: kernel
    type (constant_type) :: const
  end type problem_type

end module problem_class



