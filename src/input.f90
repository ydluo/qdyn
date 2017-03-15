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
  read(15,*) i_sigma_cpl
  read(15,*) pb%features%stress_coupling, pb%features%cohesion, pb%features%localisation ! SEISMIC
  pb%neqs = 2 + pb%features%stress_coupling + pb%features%cohesion + pb%features%localisation
  pb%kernel%has_sigma_coupling = (i_sigma_cpl == 1)
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
  allocate ( pb%sigma(n), pb%v(n), pb%theta(n),  &
             pb%a(n), pb%b(n), pb%dc(n), pb%v1(n), &
             pb%v2(n), pb%mu_star(n), pb%v_star(n), &
             pb%ot%iot(n),pb%ot%iasp(n),pb%coh(n))

  !JPA if MPI, read only the nodes of this processor
  ! SEISMIC: the CNS model allows for an initial state of stress that is
  ! read from the input file
  if (pb%i_rns_law == 3) then
    allocate ( pb%tau(n) )
    do i=1,n
      read(15,*)pb%sigma(i), pb%tau(i), pb%v(i), pb%theta(i),  &
                pb%a(i), pb%b(i), pb%dc(i), pb%v1(i), &
                pb%v2(i), pb%mu_star(i), pb%v_star(i), &
                pb%iot(i), pb%iasp(i), pb%coh(i)
    end do
  else
    do i=1,n
      read(15,*)pb%sigma(i), pb%v(i), pb%theta(i),  &
                pb%a(i), pb%b(i), pb%dc(i), pb%v1(i), &
                pb%v2(i), pb%mu_star(i), pb%v_star(i), &
                pb%iot(i), pb%iasp(i), pb%coh(i)
    end do
  endif

  ! <SEISMIC>
  ! Read input parameters for the CNS model. These parameters are (in order):
  !   a:                coefficient of logarithmic rate dependence
  !   mu_tilde_star:    reference friction coefficient at y_gr_star
  !   y_gr_star:        reference granular fow strain rate
  !   H:                dilatancy geometric factor
  !   phi0:             initial (reference) porosity
  !   IPS_const_diff:   pressure solution (temperature-dependent) constant
  !                     for diffusion controlled pressure solution creep
  !   IPS_const_diss1:  pressure solution (temperature-dependent) constant
  !                     for dissolution controlled pressure solution creep
  !   IPS_const_diss2:  pressure solution (temperature-dependent) constant
  !                     for dissolution controlled pressure solution creep
  !   w:                total thickness of the fault zone
  !   lambda:           relative thickness of localised zone
  !                     (0 = non-existent shear band, 1 = entire fault zone)
  !
  ! Note that these parameters are material (gouge) properties, and are
  ! generally not spatically uniform, and hence are allocatable
  ! See friction.f90 for a description of and references to the CNS model
  ! See user manual for detailed definitions of the above parameters (TODO)

  if (pb%i_rns_law == 3) then
    if (pb%itheta_law == 3) then
      allocate( pb%cns_params%a(n), pb%cns_params%mu_tilde_star(n), &
                pb%cns_params%y_gr_star(n), pb%cns_params%H(n), &
                pb%cns_params%phi0(n), pb%cns_params%IPS_const_diff(n), &
                pb%cns_params%IPS_const_diss1(n), pb%cns_params%IPS_const_diss2(n), &
                pb%cns_params%w(n) )

      do i=1,n
        read(15,*)pb%cns_params%a(i), pb%cns_params%mu_tilde_star(i), &
                  pb%cns_params%y_gr_star(i), pb%cns_params%H(i), &
                  pb%cns_params%phi0(i), pb%cns_params%IPS_const_diff(i), &
                  pb%cns_params%IPS_const_diss1(i), pb%cns_params%IPS_const_diss2(i), &
                  pb%cns_params%w(i)
        if (pb%cns_params%IPS_const_diff(i)*pb%cns_params%IPS_const_diss1(i) /= 0) then
          write(6,*) "input.f90: ambiguous rate-controlling mechanism for itheta_law = 3"
          write(6,*) "For each fault element, either IPS_const_diff or IPS_const_diss1 must be zero"
          stop
        endif
      end do
    else
      write(6,*) "input.f90: incompatible theta law with chosen rate-and-state law"
      write(6,*) "If the CNS friction law is chosen, theta law must be 3"
      stop
    endif
  endif

  ! End reading CNS model parameters
  ! </SEISMIC>

  if (pb%features%cohesion == 1) then
    allocate (  pb%alpha(n), pb%coh_params%alpha0(n), &
                pb%coh_params%alpha_c(n), pb%coh_params%compl(n), &
                pb%coh_params%C_star(n), pb%coh_params%E_surf(n), &
                pb%coh_params%NG_const(n), pb%dalpha_dt(n))

    do i=1,n
      read(15,*)pb%alpha(i), pb%coh_params%alpha0(i), &
                pb%coh_params%alpha_c(i), pb%coh_params%compl(i), &
                pb%coh_params%C_star(i), pb%coh_params%E_surf(i), &
                pb%coh_params%NG_const(i)
    end do
  endif

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
    allocate(pb%cns_params%lambda(n))
    do i=1,n
      pb%cns_params%lambda = 1
    end do
  endif

  if (is_MPI_parallel()) then
    call read_mesh_nodes(15,pb%mesh)
    call ot_read_stations(pb%ot)
  endif

  if (pb%mesh%dim==1) read(15,*) pb%finite
  read(15,*) pb%itheta_law
  close(15)
  write(6,*) 'Input complete'

end subroutine read_main


end module input
