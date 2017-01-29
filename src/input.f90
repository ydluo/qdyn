!qdyn input all parameters

module input

  implicit none
  private

  public :: read_main

contains
!=====================================================================
! read in all parameters
!
subroutine read_main(pb)

  use problem_class
  use mesh, only : read_mesh, mesh_get_size
  use constants, only : FFT_TYPE
  use my_mpi, only: my_mpi_tag, is_MPI_parallel

  type(problem_type), intent(inout)  :: pb

  integer :: i,n,nsta,ista,ik
  double precision :: xsta, ysta, zsta, dmin, d

  write(6,*) 'Start reading input: ...'

  !PG, if MPI then read one input file per processor
  open(unit=15, FILE='qdyn'//trim(my_mpi_tag())//'.in')
  call read_mesh(15,pb%mesh)
  write(6,*) '   Mesh input complete'

  if (pb%mesh%dim==2) then
    pb%kernel%kind =3+FFT_TYPE
    if (pb%mesh%nx < 4) then
      write(6,*) 'nx < 4, FFT disabled'
      pb%kernel%kind = 3
    endif
  else
      pb%kernel%kind = pb%mesh%dim+1
  endif

  select case (pb%kernel%kind)
    ! SEISMIC: in case of spring-block, allocate kernel%k2f or it may
    ! trigger a segmentation fault in output.f90::screen_init()
    case (1)
      allocate(pb%kernel%k2f)
      pb%kernel%k2f%finite = -1
    case(2)
      allocate(pb%kernel%k2f)
      read(15,*) pb%kernel%k2f%finite
    case(3)
      allocate(pb%kernel%k3)
    case(4)
      allocate(pb%kernel%k3f)
    case(5)
      allocate(pb%kernel%k3f2)
  end select

  read(15,*) pb%itheta_law
  read(15,*) pb%i_rns_law
  read(15,*) pb%kernel%i_sigma_cpl ! SEISMIC: this value is overridden later on?
  !read(15,*) pb%neqs ! Replace this with feature flags
  read(15,*) pb%features%stress_coupling, pb%features%cohesion ! SEISMIC
  pb%neqs = 2 + pb%features%stress_coupling + pb%features%cohesion

!JPA neqs should not be setup explicitly by the user
!    It should be inferred from the type of problem:
!    neq=2 if problem in homogeneous medium without free surface
!    neq=3 if bimaterial problem, or with free surface (for which normal stress changes
!          are coupled to slip, the 3rd variable is the normal stress)
!YD we may want to determine the neqs and other value in matlab script, then read-in
! as far as it will not make conflict with other variables/parameters
! because matlab is using more human-like language
!However, it will be safer to deal with variable/parameters here
!--?? Leave AS IS till we complete benchmark this 2D version ??--

  read(15,*)pb%ot%ntout, pb%ot%ic, pb%ox%nxout, pb%ox%nxout_dyn,    &
            pb%ox%i_ox_seq, pb%ox%i_ox_dyn
  read(15,*)pb%beta, pb%smu, pb%lam, pb%v_th

!YD This part we may want to modify it later to be able to
!impose more complicated loading/pertubation
!functions involved: problem_class/problem_type; input/read_main
!                    initialize/init_field;  derivs_all/derivs

  read(15,*)pb%Tper, pb%Aper
  read(15,*)pb%dt_try, pb%dt_max,pb%tmax, pb%acc
  read(15,*)pb%NSTOP
  read(15,*)pb%DYN_FLAG,pb%DYN_SKIP
  read(15,*)pb%DYN_M,pb%DYN_th_on,pb%DYN_th_off
  write(6,*) '  Flags input complete'

!JPA some of these arrays should be allocated in initialize.f90,
!    unless it's better to do it here to optimize memory access

  n = mesh_get_size(pb%mesh) ! number of nodes in this processor
  allocate ( pb%tau(n),     &
             pb%dtau_dt(n), pb%dsigma_dt(n), &
             pb%tau_init(n), pb%sigma(n), &
             pb%slip(n), pb%v(n), pb%dv_dt(n), &
             pb%theta(n),  pb%dtheta_dt(n),  &
             pb%a(n), pb%b(n), pb%dc(n),   &
             pb%v_pre(n), pb%v_pre2(n), &
             pb%t_rup(n), pb%tau_max(n),   &
             pb%v_max(n), pb%t_vmax(n),   &
             pb%v1(n), pb%v2(n), pb%mu_star(n),&
             pb%v_star(n), pb%theta_star(n),   &
             pb%iot(n),pb%iasp(n),pb%coh(n))

  !JPA if MPI, read only the nodes of this processor
  ! SEISMIC: Chen's model allows for an initial state of stress that is
  ! read from the input file
  if (pb%i_rns_law == 3) then
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
  ! Read input parameters for Chen's model. These parameters are (in order):
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
  !   w:                thickness of the (localised) fault zone
  !
  ! Note that these parameters are material (gouge) properties, and are
  ! generally not spatically uniform, and hence are allocatable
  ! See friction.f90 for a description of and references to Chen's model
  ! See user manual for detailed definitions of the above parameters (TODO)

  if (pb%i_rns_law == 3) then
    if (pb%itheta_law == 3) then
      allocate( pb%chen_params%a(n), pb%chen_params%mu_tilde_star(n), &
                pb%chen_params%y_gr_star(n), pb%chen_params%H(n), &
                pb%chen_params%phi0(n), pb%chen_params%IPS_const_diff(n), &
                pb%chen_params%w(n) )

      do i=1,n
        read(15,*)pb%chen_params%a(i), pb%chen_params%mu_tilde_star(i), &
                  pb%chen_params%y_gr_star(i), pb%chen_params%H(i), &
                  pb%chen_params%phi0(i), pb%chen_params%IPS_const_diff(i), &
                  pb%chen_params%w(i)
      end do
    elseif (pb%itheta_law == 4) then
      allocate( pb%chen_params%a(n), pb%chen_params%mu_tilde_star(n), &
                pb%chen_params%y_gr_star(n), pb%chen_params%H(n), &
                pb%chen_params%phi0(n), pb%chen_params%IPS_const_diss1(n), &
                pb%chen_params%IPS_const_diss2(n), pb%chen_params%w(n) )

      do i=1,n
        read(15,*)pb%chen_params%a(i), pb%chen_params%mu_tilde_star(i), &
                  pb%chen_params%y_gr_star(i), pb%chen_params%H(i), &
                  pb%chen_params%phi0(i), pb%chen_params%IPS_const_diss1(i), &
                  pb%chen_params%IPS_const_diss2(i), pb%chen_params%w(i)
      end do
    else
      write(6,*) "input.f90: incompatible theta law with chosen rate-and-state law"
      write(6,*) "If Chen's friction law is chosen, theta law must be either 3 or 4"
      write(6,*) "Pressure solution rate-limiting mechanism: 3 = diffusion, 4 = dissolution"
      stop
    endif
  endif

  ! End reading Chen's model parameters
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

  if (is_MPI_parallel()) then
    allocate(pb%mesh%x(pb%mesh%nn), pb%mesh%y(pb%mesh%nn),&
             pb%mesh%z(pb%mesh%nn), pb%mesh%dip(pb%mesh%nn))
    pb%mesh%dx = pb%mesh%Lfault/pb%mesh%nx !Check if dx is needed
!Reading mesh chunck
    do i=1,n
      read(15,*) pb%mesh%x(i),pb%mesh%y(i),pb%mesh%z(i),pb%mesh%dip(i)
    enddo

  !Finding stations in this processor
    dmin = huge(dmin)
    if (.not.(pb%ot%ic==1)) then !Reading stations, pb%ot%ic==1 is default
     open(unit=200,file='stations.dat',action='read',status='unknown')
     read(200,*) nsta
     do ista=1,nsta
       read(200,*) xsta, ysta, zsta
       sta_loop: do
        do ik=1,pb%mesh%nn
         d=sqrt((pb%mesh%x(ik)-xsta)**2+(pb%mesh%y(ik)-ysta)**2+(pb%mesh%z(ik)-zsta)**2)
         if (d<=dmin) then
           pb%ot%ic=ik
           write(6,*) 'processor: ',my_mpi_tag(),' Station found, index:',ik
           pb%station_found=.true.
           exit sta_loop
         endif
         if (ik==pb%mesh%nn) exit sta_loop
        enddo
       enddo sta_loop
      close(200)
     enddo
    endif

  endif

  close(15)
  write(6,*) 'Input complete'

end subroutine read_main


end module input
