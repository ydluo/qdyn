!define problem

module problem_class

  use fault_stress, only : kernel_type
  use mesh, only : mesh_type

  implicit none

  public

 ! timeseries outputs: at every time step, but only macroscopic quantities
  type ot_type
    double precision :: lcold,lcnew,llocnew,llocold
    integer :: unit,ic,ntout,ivmax
!For MPI
    integer :: ivmaxglob
  end type ot_type

 ! snapshot outputs: at every fault point, but only at few selected times
  type ox_type
    integer :: count,dyn_count2,unit,nxout,nxout_dyn,countglob,&
                i_ox_seq, i_ox_dyn, dyn_stat, dyn_stat2, dyn_count
  end type ox_type

  ! SEISMIC: definition of structure that holds Chen's model parameters
  ! See input.f90 for a description of the parameters
  type chen_type
    double precision, dimension(:), allocatable :: &
      a, mu_tilde_star, IPS_const_diff, IPS_const_diss1, IPS_const_diss2, H, w, y_gr_star, phi0
  end type chen_type
  ! End of Chen's model structure

  type problem_type
    double precision, dimension(:), allocatable :: &
      tau, dtau_dt, tau_init, &
      sigma, dsigma_dt, v_pre, v_pre2, &
      tau_max, t_rup, v_max, t_vmax,  &
      slip, v, dv_dt, theta, dtheta_dt,  &
      a, b, dc, v1, v2, mu_star, v_star, &
      theta_star, iot, iasp, coh
!For MPI
    double precision, dimension(:), allocatable :: &
      tau_glob, dtau_dt_glob,&
      sigma_glob,&
      slip_glob, v_glob, dv_dt_glob, theta_glob, dtheta_dt_glob,&
      tau_max_glob, t_rup_glob, v_max_glob, t_vmax_glob
    double precision :: vmaxglob,sigma_vmaxglob

    double precision :: pot, pot_rate, pot_pre
    double precision :: beta=0d0, smu=0d0, lam=0d0, zimpedance=0d0, v_th
    logical :: station_found=.false. !Found station

!YD This part we may want to modify it later to be able to
!impose more complicated loading/pertubation
!functions involved: problem_class/problem_type; input/read_main
!                    initialize/init_field;  derivs_all/derivs

    double precision :: Tper=0d0, Aper=0d0, Omper=0d0
    double precision :: time=0d0
    integer :: itheta_law,i_rns_law, neqs

    double precision :: dt_try, dt_did, dt_next, dt_max=0d0, tmax, acc
    double precision :: DYN_M,DYN_th_on,DYN_th_off
    integer :: NSTOP,itstop,it,DYN_FLAG,DYN_SKIP

    type (mesh_type) :: mesh
    type (ot_type) :: ot
    type (ox_type) :: ox
    type (kernel_type) :: kernel
    ! SEISMIC: add structure that holds Chen model parameters
    type (chen_type) :: chen_params
  end type problem_type

end module problem_class
