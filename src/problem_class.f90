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
  end type ot_type

 ! snapshot outputs: at every fault point, but only at few selected times
  type ox_type
    integer :: count,dyn_count2,unit,nxout,nxout_dyn,        &
                i_ox_seq, i_ox_dyn, dyn_stat, dyn_stat2, dyn_count 
  end type ox_type

  type problem_type
    double precision, dimension(:), allocatable :: &
      tau, dtau_dt, tau_init, &
      sigma, dsigma_dt, v_pre, v_pre2, &
      tau_max, t_rup, v_max, t_vmax,  &
      slip, v, dv_dt, theta, dtheta_dt,  &
      a, b, dc, v1, v2, mu_star, v_star, theta_star, iot, iasp
    double precision :: pot, pot_rate, pot_pre
    double precision :: beta=0d0, smu=0d0, lam=0d0, zimpedance=0d0, v_th

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
  end type problem_type

end module problem_class



