module unittests

  use constants, only : SOLVER
  use problem_class, only : problem_type
  use mesh, only : mesh_get_size
  use solver, only : init_rk45
  use friction
  ! use friction_cns
  use fault_stress, only : init_kernel

  use unittests_rsf
  ! use unittests_cns

  implicit none
  private

  public :: init_tests

contains

!===============================================================================
! Initiate the unit testing framework
!===============================================================================
subroutine init_tests(pb)
  type(problem_type) :: pb

  pb%test_mode = .true.

  write(6,*) ""
  write(6,*) "-----------------------------------------------------------------"
  write(6,*) "Initiating QDYN unit test suite"
  write(6,*) ""

  ! Define assert_close function

  ! Spring-block configuration
  call initiate_springblock(pb)

  ! Allocate mesh variables/parameters
  call allocate_mesh(pb)
  call initiate_parameters(pb)
  call init_kernel(pb%lam, pb%smu, pb%mesh, pb%kernel, &
                   pb%D, pb%H, pb%i_sigma_cpl, pb%finite)

  write(6,*) ""
  write(6,*) "Test suite initialised successfully"

  ! Test RSF
  call test_rsf_friction(pb)

  ! Test CNS (analytical solutions Chen?)
  ! call test_cns_friction(pb)

  ! Test TP (analytical soln. Noda)

  ! Test stress calculations (ask analytical soln. Pablo)

  ! Test MPI

  write(6,*) ""
  write(6,*) "Finalised QDYN unit test suite"
  write(6,*) "-----------------------------------------------------------------"
end subroutine init_tests

!===============================================================================
! Set the parameters corresponding to a 0-D spring-block configuration
!===============================================================================
subroutine initiate_springblock(pb)
  type(problem_type) :: pb
  ! Spring-block configuration
  pb%mesh%dim = 0
  pb%mesh%nn = 1
  pb%mesh%Lfault = 1
  pb%mesh%W = 1

  write(6,*) " * Spring-block configuration initialised"
end subroutine initiate_springblock

!===============================================================================
! Allocate mesh parameters for both RSF and CNS
!===============================================================================
subroutine allocate_mesh(pb)
  type(problem_type) :: pb
  integer :: n

  n = mesh_get_size(pb%mesh) ! number of nodes in this processor

  ! Allocate general parameters
  allocate ( pb%tau(n), pb%sigma(n), pb%v(n), pb%theta(n), pb%theta2(n),  &
             pb%v_star(n), pb%ot%iot(n), pb%ot%iasp(n), &
             pb%dc(n), pb%coh(n), pb%v_pl(n) )

  ! Allocate RSF parameters
  allocate ( pb%a(n), pb%b(n), pb%v1(n), pb%v2(n), pb%mu_star(n))

  ! Allocate CNS parameters
  ! allocate( pb%cns_params%a_tilde(n), pb%cns_params%mu_tilde_star(n), &
  !           pb%cns_params%y_gr_star(n), pb%cns_params%H(n), &
  !           pb%cns_params%phi_c(n), pb%cns_params%phi0(n), &
  !           pb%cns_params%IPS_const_diff(n), pb%cns_params%IPS_const_diss1(n), &
  !           pb%cns_params%IPS_const_diss2(n), pb%cns_params%L(n), &
  !           pb%cns_params%lambda(n) )

  allocate ( pb%dtau_dt(n), pb%slip(n), pb%theta_star(n) )

  write(6,*) " * Mesh variables/parameters allocated"
end subroutine allocate_mesh

!===============================================================================
! Allocate simulation settings and mesh parameters
!===============================================================================
subroutine initiate_parameters(pb)
  type(problem_type) :: pb

  ! Solver settings
  pb%dt_try = 1e-1
  pb%dt_max = 0.d0
  pb%acc = 1e-7

  ! Simulation features
  pb%features%localisation = 0
  pb%features%stress_coupling = 0
  pb%i_sigma_cpl = 0
  pb%features%tp = 0
  pb%neqs = 2
  pb%finite = 1

  ! Periodic loading
  pb%Omper = 0.d0
  pb%Aper = 0.d0
  pb%Tper = 31536000.d0

  ! Damage model
  pb%D = 0.d0
  pb%H = 0.d0

  ! Elastic moduli and radiation damping
  pb%lam = 30e9
  pb%smu = 30e9
  pb%beta = 3e3
  pb%zimpedance = 0.5d0*pb%smu/pb%beta

  !
  pb%coh = 0.d0
  pb%time = 0.d0

  write(6,*) " * Global parameters initialised"

end subroutine initiate_parameters

end module unittests
