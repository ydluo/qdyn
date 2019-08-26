module initialize

  implicit none
  private

  public :: init_all

contains

!=============================================================
subroutine init_all(pb)

  use problem_class
  use mesh, only : init_mesh, mesh_get_size
  use constants, only : PI, SOLVER_TYPE
  use my_mpi, only: is_MPI_master
  use fault_stress, only : init_kernel
  use output, only : ot_init, ox_init
  use friction, only : set_theta_star, friction_mu
  use friction_cns, only : compute_velocity, dphi_dt
  use solver, only : init_rk45
  use diffusion_solver, only: init_tp
  use ode_rk45_2, only : init_rk45_2
!!$  use omp_lib

  type(problem_type), intent(inout) :: pb

  integer :: n
!  integer :: TID, NTHREADS

  call init_mesh(pb%mesh)

 ! number of equations
  pb%neqs = 2 + pb%features%localisation
  if (pb%features%stress_coupling == 1 .and. pb%mesh%dim == 2) then
    pb%neqs = pb%neqs + 1
  endif

 ! dt_max & perturbation
 ! if periodic loading, set time step smaller than a fraction of loading period
  if (pb%Aper /= 0.d0 .and. pb%Tper > 0.d0) then
    if (pb%dt_max > 0.d0) then
      pb%dt_max = min(pb%dt_max,0.2d0*pb%Tper)
    else
      pb%dt_max = 0.2d0*pb%Tper
    endif
  endif
  if (pb%Tper > 0.d0) then
    pb%Omper = 2.d0*PI/pb%Tper
  else
    pb%Omper = 0.d0
  endif

 ! impedance
  if (pb%beta > 0d0) then
    pb%zimpedance = 0.5d0*pb%smu/pb%beta
  else
    pb%zimpedance = 0.d0
  endif
  if (is_mpi_master()) write(6,*) 'Impedance = ', pb%zimpedance
  !---------------------- impedance ------------------

  !---------------------- ref_value ------------------
  n = mesh_get_size(pb%mesh)
  ! SEISMIC: initialise pb%tau in input.f90 to be compatible with CNS model
  allocate ( pb%dtau_dt(n), pb%slip(n), pb%theta_star(n) )
  pb%slip = 0d0
  call set_theta_star(pb)

  ! SEISMIC: initialise thermal pressurisation model (diffusion_solver.f90)
  if (pb%features%tp == 1) then
    call init_tp(pb)
    write(6,*) "Spectral mesh initiated"
  else
    allocate(pb%tp%P(pb%mesh%nn))
    pb%tp%P = 0d0
  endif

  ! SEISMIC: the CNS model has the initial shear stress defined in the
  ! input file, so we can skip the initial computation of friction
  if (pb%i_rns_law /= 3) then
    pb%tau = (pb%sigma - pb%tp%P) * friction_mu(pb%v,pb%theta,pb) + pb%coh
  endif
  if (pb%i_rns_law == 3) then
    pb%v = compute_velocity(pb%tau, pb%sigma-pb%tp%P, pb%theta, pb%theta2, pb)
    if (pb%features%tp == 1) then
      call dphi_dt( pb%v, pb%tau, pb%sigma-pb%tp%P, pb%theta, pb%theta2, &
                    pb%dtheta_dt, pb%dtheta2_dt, pb)
    endif
  endif

  call init_kernel( pb%lam, pb%smu, pb%mesh, pb%kernel, pb%D, pb%H, &
                    pb%i_sigma_cpl, pb%finite, pb%test_mode)
  call ot_init(pb)
  call ox_init(pb)

  ! SEISMIC: initialise Runge-Kutta ODE solver, if selected
  if (SOLVER_TYPE == 2) then
    call init_rk45(pb)
  endif
  if (SOLVER_TYPE == 3) then
    call init_rk45_2(pb)
  endif

  if (is_mpi_master()) write(6,*) 'Initialization completed'

  ! Info about threads
!!$OMP PARALLEL PRIVATE(NTHREADS, TID)
!!$  TID = OMP_GET_THREAD_NUM()
!!$  write(6,*) 'Thread index = ', TID
!!$OMP BARRIER
!!$  if (TID == 0) then
!!$    NTHREADS = OMP_GET_NUM_THREADS()
!!$    write(6,*) 'Total number of threads = ', NTHREADS
!!$  end if
!!$OMP END PARALLEL

end subroutine init_all

end module initialize
