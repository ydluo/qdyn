!define problem

module problem_class

  use fault_stress, only : kernel_type
  use mesh, only : mesh_type, spectral_mesh_type

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
    integer :: ntout=0
    integer :: count,dyn_count2,unit,nxout,nxout_dyn,countglob,&
                i_ox_seq, i_ox_dyn, dyn_stat, dyn_stat2, dyn_count
    double precision :: pot_pre
    double precision, dimension(:), allocatable :: tau_max, t_rup, v_max, t_vmax
  end type ox_type

  ! SEISMIC: structure that holds the CNS model parameters
  ! See input.f90 for a description of the parameters
  type cns_type
    double precision, dimension(:), allocatable :: &
      a_tilde, mu_tilde_star, H, L, & ! Granular flow and microstructural params
      y_gr_star, phi_c, phi0, lambda, &
      A, n, m, A_bulk, n_bulk, m_bulk ! Creep mechanism parameters
    integer :: N_creep=-1 ! Number of active creep mechanisms
  end type cns_type
  ! End of the CNS model structure

  ! SEISMIC: structure that holds the thermal pressurisation (TP) model parameters
  ! See input.f90 for a description of the parameters
  ! Spectral mesh parameters (Dlogl, lw_max, Nl) are hard-coded in mesh.f90
  type tp_type
    type (spectral_mesh_type) :: mesh
    double precision, dimension(:), allocatable :: &
      rhoc, inv_rhoc, beta, eta, k_t, k_p, l, w, inv_w, P, T, P_a, T_a, &
      Pi, Theta, PiTheta, Omega, dP_dt, Theta_prev, PiTheta_prev, &
      alpha_th, alpha_hy, Lam, Lam_prime, Lam_T, phi_b, dilat_factor
    double precision :: t_prev=0d0
  end type tp_type
  ! End of the TP model structure

  ! SEISMIC: requested features structure
  ! stress_coupling: normal stress variations at subduction interfaces
  ! tp: thermal pressurisation
  ! cohesion: time-dependent cohesion (CNS model)
  ! localisation: localisation of deformation (CNS model)
  type features_type
    integer :: stress_coupling, tp, localisation
  end type features_type
  ! End of features structure

  ! SEISMIC: Runge-Kutta-Fehlberg solver parameters
  type rk45_type
    integer :: iflag
    double precision, dimension(:), allocatable :: work
    integer, dimension(:), allocatable :: iwork
  end type rk45_type
  ! End of features structure

  ! SEISMIC: Runge-Kutta-Fehlberg solver parameters
  type rk45_2_type
    double precision, dimension(5) :: t_coeffs
    double precision, dimension(5) :: error_coeffs
    double precision, dimension(19) :: coeffs
    double precision, dimension(:), allocatable :: k1, k2, k3, k4, k5, k6
    double precision, dimension(:), allocatable :: dy1, dy2, dy3, dy4, dy5, dy6
    double precision, dimension(:), allocatable :: y_new, h, e
  end type rk45_2_type
  ! End of features structure

  ! Structure for (unit)testing
  type test_type
    logical :: test_mode = .false.
    logical :: test_passed = .true.
  end type test_type
  ! End of testing structure

  type problem_type
    type (mesh_type) :: mesh
    type (kernel_type) :: kernel
   ! Basic variables
    double precision, dimension(:), allocatable :: tau, dtau_dt, dtheta_dt, dtheta2_dt
    double precision, dimension(:), allocatable :: sigma, slip, v, theta, theta2, alpha
   ! Boundary conditions
    integer :: i_sigma_cpl=0, finite=0
   ! Friction properties
    double precision, dimension(:), allocatable :: a, b, dc, v1, v2, mu_star, v_star, theta_star, coh, v_pl
    integer :: itheta_law=1, i_rns_law=1, neqs=2
   ! Elastic properties
    double precision :: beta=0d0, smu=0d0, lam=0d0, D=0d0, H=0d0, zimpedance=0d0
   ! Periodic loading
    double precision :: Tper=0d0, Aper=0d0, Omper=0d0
   ! Time solver
    double precision :: t_prev=0d0, time=0d0, dt_try=0d0, dt_did=0d0, dt_next=0d0, dt_max=0d0, tmax, acc
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
    ! Flag for unit testing
    logical :: test_mode = .false.

    ! SEISMIC: added structures
    type (cns_type) :: cns_params
    type (tp_type) :: tp
    type (features_type) :: features
    type (rk45_type) :: rk45
    type (rk45_2_type) :: rk45_2
    type (test_type) :: test
  end type problem_type

end module problem_class
