!define problem

module problem_class

  use fault_stress, only : kernel_type
  use mesh, only : mesh_type

  implicit none

  public

 ! timeseries outputs: at every time step, but only macroscopic quantities
  type ot_type
    double precision :: lcold,lcnew,llocnew,llocold, v_th, pot, pot_rate
    double precision, dimension(:), allocatable :: xsta, ysta, zsta, v_pre, v_pre2
    integer :: unit=0,ic=0,ntout=0,ivmax=0
    integer :: ivmaxglob=0 ! For MPI
    integer, dimension(:), allocatable :: iot, iasp
  end type ot_type

 ! snapshot outputs: at every fault point, but only at few selected times
  type ox_type
    integer :: count,dyn_count2,unit,nxout,nxout_dyn,countglob,&
                i_ox_seq, i_ox_dyn, dyn_stat, dyn_stat2, dyn_count
    double precision :: pot_pre
    double precision, dimension(:), allocatable :: tau_max, t_rup, v_max, t_vmax
  end type ox_type

  ! SEISMIC: definition of structure that holds the CNS model parameters
  ! See input.f90 for a description of the parameters
  type cns_type
    double precision, dimension(:), allocatable :: &
      a, mu_tilde_star, IPS_const_diff, IPS_const_diss1, IPS_const_diss2, H, w, y_gr_star, phi0, lambda, &
      IPS_const_diff_bulk, IPS_const_diss1_bulk, IPS_const_diss2_bulk
  end type cns_type
  ! End of the CNS model structure

  ! SEISMIC: definition of structure that holds the cohesion model parameters
  ! See input.f90 for a description of the parameters
  type cohesion_type
    double precision, dimension(:), allocatable :: &
      alpha0, alpha_c, compl, C_star, E_surf, NG_const
  end type cohesion_type
  ! End of cohesion model structure

  ! SEISMIC: requested features structure (normal stress coupling, cohesion)
  type features_type
    integer :: stress_coupling, cohesion, localisation
  end type features_type
  ! End of features structure

  type problem_type
    type (mesh_type) :: mesh
    type (kernel_type) :: kernel
   ! Basic variables
    double precision, dimension(:), allocatable :: tau, dtau_dt, sigma, slip, v, theta, theta2, alpha
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
    double precision :: vmaxglob
   ! For outputs
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

    ! SEISMIC: add structure that holds the CNS model parameters
    type (cns_type) :: cns_params
    type (cohesion_type) :: coh_params
    type (features_type) :: features
  end type problem_type

end module problem_class
