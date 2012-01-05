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
    integer :: kind = 0
    integer :: nn
    double precision :: dx, Lfault, W
    double precision, allocatable :: x(:) 
  end type mesh_type

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

  type problem_type
    double precision, dimension(:), allocatable :: &
      tau, dtau_dt, tau_init, &
      sigma, &
      slip, v, dv_dt, theta, dtheta_dt,  &
      a, b, dc, v1, v2, mu_star, v_star, theta_star
    double precision :: pot, pot_rate
    double precision :: beta, smu, zimpedance
    double precision :: Tper, Aper, Omper
    double precision :: time
    integer :: itheta_law, neqs

    double precision :: dt_try, dt_did, dt_next, dt_max, tmax, acc
    integer :: NSTOP,itstop,it

    type (mesh_type) :: mesh
    type (ot_type) :: ot
    type (ox_type) :: ox
    type (kernel_type) :: kernel
  end type problem_type

end module problem_class



