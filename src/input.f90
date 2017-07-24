!qdyn input all parameters

module input

  implicit none
  private

  public :: read_main

contains
!=====================================================================
! Read input file qdyn.in
! If MPI, there is one qdyn*.in file per processor
!
subroutine read_main(pb)

  use problem_class
  use mesh, only : read_mesh_parameters, read_mesh_nodes, mesh_get_size
  use output, only : ot_read_stations
  use my_mpi, only: my_mpi_tag, is_MPI_parallel

  type(problem_type), intent(inout)  :: pb

  integer :: i,n

  write(6,*) 'Start reading input ...'

  open(unit=15, file='qdyn'//trim(my_mpi_tag())//'.in', action='read')

  call read_mesh_parameters(15,pb%mesh)
  write(6,*) '   Mesh input complete'

  if (pb%mesh%dim==1) read(15,*) pb%finite
  read(15,*) pb%itheta_law
  read(15,*) pb%i_rns_law
  read(15,*) pb%i_sigma_cpl
  ! SEISMIC: various simulation features can be turned on (1) or off (0)
  read(15,*)  pb%features%stress_coupling, pb%features%tp, &
              pb%features%cohesion, pb%features%localisation
  read(15,*)pb%ot%ntout, pb%ot%ic, pb%ox%nxout, pb%ox%nxout_dyn,    &
            pb%ox%i_ox_seq, pb%ox%i_ox_dyn
  read(15,*) pb%beta, pb%smu, pb%lam, pb%D, pb%H, pb%ot%v_th
  read(15,*)pb%Tper, pb%Aper
  read(15,*)pb%dt_try, pb%dt_max,pb%tmax, pb%acc
  read(15,*)pb%NSTOP
  read(15,*)pb%DYN_FLAG,pb%DYN_SKIP
  read(15,*)pb%DYN_M,pb%DYN_th_on,pb%DYN_th_off
  write(6,*) '  Flags input complete'

  n = mesh_get_size(pb%mesh) ! number of nodes in this processor
  allocate ( pb%tau(n), pb%sigma(n), pb%v(n), pb%theta(n),  &
             pb%v_star(n), pb%ot%iot(n), pb%ot%iasp(n), &
             pb%dc(n), pb%coh(n))

  !JPA if MPI, read only the nodes of this processor
  ! <SEISMIC>
  ! Read input parameters for the CNS model. These parameters are (in order):
  !   sigma:            applied normal stress (effective if no TP)
  !   tau:              initial shear stress
  !   theta:            initial porosity (phi in CNS formulation)
  !   v_star:           load-point velocity
  !   a_tilde:          coefficient of logarithmic rate dependence
  !   mu_tilde_star:    reference friction coefficient at y_gr_star
  !   y_gr_star:        reference granular fow strain rate
  !   H:                dilatancy geometric factor
  !   phi_c:            critical state porosity
  !   phi0:             lower cut-off porosity
  !   IPS_const_diff:   pressure solution (temperature-dependent) constant
  !                     for diffusion controlled pressure solution creep
  !   IPS_const_diss1:  pressure solution (temperature-dependent) constant
  !                     for dissolution controlled pressure solution creep
  !   IPS_const_diss2:  pressure solution (temperature-dependent) constant
  !                     for dissolution controlled pressure solution creep
  !   w:                total thickness of the fault zone
  !
  ! Note that these parameters are material (gouge) properties, and are
  ! generally not spatically uniform, and hence are allocatable
  ! See friction.f90 for a description of and references to the CNS model
  ! See user manual for detailed definitions of the above parameters (TODO)
  ! </SEISMIC>

  ! If the CNS model is selected
  if (pb%i_rns_law == 3) then
    allocate( pb%cns_params%a_tilde(n), pb%cns_params%mu_tilde_star(n), &
              pb%cns_params%y_gr_star(n), pb%cns_params%H(n), &
              pb%cns_params%phi_c(n), pb%cns_params%phi0(n), &
              pb%cns_params%IPS_const_diff(n), pb%cns_params%IPS_const_diss1(n), &
              pb%cns_params%IPS_const_diss2(n), pb%cns_params%L(n) )
    do i=1,n
      read(15,*)pb%sigma(i), pb%tau(i), pb%theta(i), pb%v_star(i),  &
                pb%cns_params%a_tilde(i), pb%cns_params%mu_tilde_star(i), &
                pb%cns_params%y_gr_star(i), pb%cns_params%H(i), &
                pb%cns_params%phi_c(i), pb%cns_params%phi0(i), &
                pb%cns_params%IPS_const_diff(i), pb%cns_params%IPS_const_diss1(i), &
                pb%cns_params%IPS_const_diss2(i), pb%cns_params%L(i), &
                pb%ot%iot(i), pb%ot%iasp(i)
      pb%dc(i) = 1.0
    end do
  ! Else, the RSF model is selected
  else
    allocate ( pb%a(n), pb%b(n), pb%v1(n), &
               pb%v2(n), pb%mu_star(n))
    do i=1,n
      read(15,*)pb%sigma(i), pb%v(i), pb%theta(i),  &
                pb%a(i), pb%b(i), pb%dc(i), pb%v1(i), &
                pb%v2(i), pb%mu_star(i), pb%v_star(i), &
                pb%ot%iot(i), pb%ot%iasp(i), pb%coh(i)
    end do
  endif

  ! <SEISMIC>
  ! Read input parameters for the time-dependent cohesion model (CNS only).
  ! These parameters and corresponding units are (in order):
  !
  if (pb%features%cohesion == 1) then
    allocate (  pb%alpha(n), pb%coh_params%alpha0(n), &
                pb%coh_params%alpha_c(n), pb%coh_params%compl(n), &
                pb%coh_params%C_star(n), pb%coh_params%E_surf(n), &
                pb%coh_params%NG_const(n))
    do i=1,n
      read(15,*)pb%alpha(i), pb%coh_params%alpha0(i), &
                pb%coh_params%alpha_c(i), pb%coh_params%compl(i), &
                pb%coh_params%C_star(i), pb%coh_params%E_surf(i), &
                pb%coh_params%NG_const(i)
    end do
  else
    allocate(pb%alpha(n))
    do i=1,n
      pb%alpha(i) = 0
    end do
  endif
  ! End reading cohesion model parameters
  ! </SEISMIC>

  ! <SEISMIC>
  ! Read input parameters for the localisation model (CNS only).
  ! These parameters and corresponding units are (in order):
  !
  if (pb%features%localisation == 1) then
    allocate (  pb%cns_params%lambda(n), pb%theta2(n), &
                pb%cns_params%IPS_const_diff_bulk(n), pb%cns_params%IPS_const_diss1_bulk(n), &
                pb%cns_params%IPS_const_diss2_bulk(n) )
    do i=1,n
      read(15,*)pb%cns_params%lambda(i), pb%theta2(i), &
                pb%cns_params%IPS_const_diff_bulk(i), pb%cns_params%IPS_const_diss1_bulk(i), &
                pb%cns_params%IPS_const_diss2_bulk(i)
    end do
  else
    allocate(pb%cns_params%lambda(n), pb%theta2(n))
    do i=1,n
      pb%cns_params%lambda(i) = 1
      pb%theta2(i) = 0
    end do
  endif
  ! End reading localisation model parameters
  ! </SEISMIC>

  ! <SEISMIC>
  ! Read input parameters for the thermal pressurisation (TP) model.
  ! These parameters and corresponding units are (in order):
  !
  !   rhoc:   density times specific heat capacity of host rock [J/K/m^3]
  !   beta:   bulk compressibility [1/Pa]
  !   eta:    dynamic viscocity [Pa s]
  !   w:      half-width of shear distribution profile [m]
  !   k_t:    thermal conductivity [J/s/K/m]
  !   k_p:    intrinsic hydraulic permeability [m^2]
  !   l:      nett thermal expansion coefficient [1/K]
  !   P_a:    ambient fluid pressure [Pa]
  !   T_a:    ambient temperature [K]
  !
  if (pb%features%tp == 1) then
    allocate (  pb%tp%rhoc(n), pb%tp%beta(n), pb%tp%eta(n), pb%tp%w(n), &
                pb%tp%k_t(n), pb%tp%k_p(n), pb%tp%l(n), &
                pb%tp%P_a(n), pb%tp%T_a(n) )
    do i=1,n
      read(15,*)pb%tp%rhoc(i), pb%tp%beta(i), pb%tp%eta(i), pb%tp%w(i), &
                pb%tp%k_t(i), pb%tp%k_p(i), pb%tp%l(i), &
                pb%tp%P_a(i), pb%tp%T_a(i)
    end do
  endif
  ! End reading TP model parameters
  ! </SEISMIC>

  if (is_MPI_parallel()) then
    call read_mesh_nodes(15,pb%mesh)
    call ot_read_stations(pb%ot)
  endif

  close(15)
  write(6,*) 'Input complete'

end subroutine read_main


end module input
