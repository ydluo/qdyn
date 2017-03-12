!define problem 

module problem_class
 
  use fault_stress, only : kernel_type
  use mesh, only : mesh_type

  implicit none

  public

 ! timeseries outputs: at every time step, but only macroscopic quantities
  type ot_type
    double precision :: lcold,lcnew,llocnew,llocold
    integer :: unit=0,ic=0,ntout=0,ivmax=0
    integer :: ivmaxglob=0 ! For MPI
  end type ot_type

 ! snapshot outputs: at every fault point, but only at few selected times
  type ox_type
    integer :: count,dyn_count2,unit,nxout,nxout_dyn,countglob,&
                i_ox_seq, i_ox_dyn, dyn_stat, dyn_stat2, dyn_count 
  end type ox_type

  type problem_type
    type (mesh_type) :: mesh
    type (kernel_type) :: kernel
   ! Basic variables
    double precision, dimension(:), allocatable :: tau, dtau_dt, sigma, slip, v, theta
   ! Boundary conditions
    integer :: i_sigma_cpl=0, finite=0
   ! Friction properties
    double precision, dimension(:), allocatable :: a, b, dc, v1, v2, mu_star, v_star, theta_star, coh
    integer :: itheta_law=1, i_rns_law=1, neqs=2
   ! Elastic properties
    double precision :: beta=0d0, smu=0d0, lam=0d0, D=0d0, H=0d0, zimpedance=0d0
   ! Periodic loading
    double precision :: Tper=0d0, Aper=0d0, Omper=0d0
   ! Time solver
    double precision :: time=0d0, dt_try=0d0, dt_did=0d0, dt_next=0d0, dt_max=0d0, tmax, acc
    integer :: NSTOP, itstop=-1, it=0
   ! For outputs
    double precision, dimension(:), allocatable :: v_pre, v_pre2, tau_max, t_rup, v_max, t_vmax
    integer, dimension(:), allocatable :: iot, iasp
    double precision :: vmaxglob, v_th, pot, pot_rate, pot_pre
    type (ot_type) :: ot
    type (ox_type) :: ox
   ! For MPI outputs
    double precision, dimension(:), allocatable :: &
      tau_glob, dtau_dt_glob, sigma_glob,&
      slip_glob, v_glob, theta_glob, &
      tau_max_glob, t_rup_glob, v_max_glob, t_vmax_glob
   ! QSB
    double precision :: DYN_M,DYN_th_on,DYN_th_off
    integer :: DYN_FLAG,DYN_SKIP
  end type problem_type

end module problem_class



