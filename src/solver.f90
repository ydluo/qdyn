! Solve_master

module solver

  use problem_class, only : problem_type

  implicit none
  private

  integer(kind=8), save :: iktotal

  public  :: solve, init_lsoda, init_rk45

contains

!=====================================================================
! Master Solver
!
subroutine solve(pb)

  use output, only : screen_init, screen_write, ox_write, ot_write
  use my_mpi, only : is_MPI_parallel, is_mpi_master, finalize_mpi

  type(problem_type), intent(inout)  :: pb

  if (is_mpi_master()) call screen_init(pb)
  call update_field(pb)
  call screen_write(pb)
  call ox_write(pb)

  iktotal=0
  ! Time loop
  do while (pb%it /= pb%itstop)
!  do while (pb%it+1 /= pb%itstop)
    pb%it = pb%it + 1
!   if (is_mpi_master()) write(6,*) 'it:',pb%it
    call do_bsstep(pb)
! if stress exceeds yield call Coulomb_solver ! JPA Coulomb quick and dirty
!                         or (cleaner version) do linear adjustment of
!                         timestep then redo bsstep
!                         or (cleanest version) iterate tiemstep adjustment and
!                         bsstep until stress is exactly equal to yield
    call update_field(pb)
    ! SEISMIC: shouldn't there be a specific output time step for ot_write
    ! as well? I.e. to have something like ntout for screen_write and one
    ! for ot_write
    call ot_write(pb)
    call check_stop(pb)   ! here itstop will change
!--------Output onestep to screen and ox file(snap_shot)
! if(mod(pb%it-1,pb%ot%ntout) == 0 .or. pb%it == pb%itstop) then
    if(mod(pb%it,pb%ot%ntout) == 0 .or. pb%it == pb%itstop) then
!      if (is_mpi_master()) write(6,*) 'it:',pb%it,'iktotal=',iktotal,'pb%time=',pb%time
      call screen_write(pb)
    endif
! if (is_mpi_master()) call ox_write(pb)
    call ox_write(pb)
  enddo

  if (is_MPI_parallel()) call finalize_mpi()

end subroutine solve



!=====================================================================
! pack, do bs_step and unpack
!
! IMPORTANT NOTE : between pack/unpack pb%v & pb%theta are not up-to-date
! SEISMIC IMPORTANT NOTE: when the CNS model is used, pb%tau is not up-to-date
!
subroutine do_bsstep(pb)

  use derivs_all
  use ode_bs
  use ode_lsoda_main, only: dlsoda
  use ode_rk45, only: rkf45_d
  use output, only : screen_write, ox_write, ot_write
  use constants, only : USE_RK_SOLVER

  type(problem_type), intent(inout) :: pb

  double precision :: t0
  double precision, dimension(pb%neqs*pb%mesh%nn) :: yt, dydt, yt_scale
  integer :: ik, ind_stress_coupling, ind_localisation

  ! Pack v, theta into yt
  ! yt(2::pb%neqs) = pb%v(pb%rs_nodes) ! JPA Coulomb

  ! SEISMIC: define the indices of yt and dydt based on which
  ! features are requested (defined in input file)
  ind_stress_coupling = 2 + pb%features%stress_coupling
  ind_localisation = ind_stress_coupling + pb%features%localisation

  ! SEISMIC: in the case of the CNS model, solve for tau and not v
  if (pb%i_rns_law == 3) then   ! SEISMIC: CNS model
    yt(2::pb%neqs) = pb%tau
  else  ! SEISMIC: not CNS model (i.e. rate-and-state)
    yt(2::pb%neqs) = pb%v
  endif
  yt(1::pb%neqs) = pb%theta
  ! SEISMIC NOTE/WARNING: I don't know how permanent this temporary solution is,
  ! but in case it gets fixed more permanently, derivs_all.f90 needs adjustment
  if (pb%features%stress_coupling == 1) then           ! Temp solution for normal stress coupling
    yt(ind_stress_coupling::pb%neqs) = pb%sigma
  endif
  if (pb%features%localisation == 1) then
    yt(ind_localisation::pb%neqs) = pb%theta2
  endif

  ! SEISMIC: user-defined switch (constants.f90) to use either the Runge-Kutta
  ! method, or Bulirsch-Stoer
  if (USE_RK_SOLVER .eqv. .true.) then

    ! Set-up Runge-Kutta solver
    pb%rk45%iflag = -2 ! Reset to one-step mode each call
    pb%t_prev = pb%time

    ! Call Runge-Kutta solver routine
    call rkf45_d( derivs_rk45, pb%neqs*pb%mesh%nn, yt, pb%time, 2*pb%tmax, &
                  pb%acc, 0d0, pb%rk45%iflag, pb%rk45%work, pb%rk45%iwork)

    ! Set time step
    pb%dt_did = pb%time - pb%t_prev

    ! Basic error checking. See description of rkf45_d in ode_rk45.f90 for details
    select case (pb%rk45%iflag)
    case (3)
      write (6,*) "RK45 error [3]: relative error tolerance too small"
      stop
    case (4)
      !write (6,*) "RK45 warning [4]: integration took more than 3000 derivative evaluations"
    case (5)
      write (6,*) "RK45 error [5]: solution vanished, relative error test is not possible"
      stop
    case (6)
      write (6,*) "RK45 error [6]: requested accuracy could not be achieved"
      stop
    case (8)
      call ot_write(pb)
      call screen_write(pb)
      call ox_write(pb)
      write (6,*) "RK45 error [8]: invalid input parameters"
      stop
    end select
  else
    ! Use Bulirsch-Stoer method

    ! this update of derivatives is only needed to set up the scaling (yt_scale)
    call derivs(pb%time,yt,dydt,pb)
    yt_scale=dabs(yt)+dabs(pb%dt_try*dydt)
    ! One step
    call bsstep(yt,dydt,pb%neqs*pb%mesh%nn,pb%time,pb%dt_try,pb%acc,yt_scale,pb%dt_did,pb%dt_next,pb,ik)

    ! SEISMIC NOTE: what is happening here?
    if (pb%dt_max >  0.d0) then
      pb%dt_try = min(pb%dt_next,pb%dt_max)
    else
      pb%dt_try = pb%dt_next
    endif
  endif

  ! SEISMIC: below here is some code that calls the LSODA solver method
  ! In development, do not use!
  ! pb%lsoda%istate = 2
  ! call dlsoda(  derivs_lsoda, pb%lsoda%neq, yt, pb%time, pb%lsoda%tout, &
  !               pb%lsoda%itol, pb%lsoda%rtol, pb%lsoda%atol, pb%lsoda%itask, &
  !               pb%lsoda%istate, pb%lsoda%iopt, pb%lsoda%rwork, pb%lsoda%lrw, &
  !               pb%lsoda%iwork, pb%lsoda%liw, jac_lsoda, pb%lsoda%jt)
  !
  ! if (pb%lsoda%istate /= 2) then
  !   write (6,*) "Next iteration by LSODA solver failed!"
  !   write (6,*) "Current value of istate: ", pb%lsoda%istate
  !   write (6,*) "Will now terminate..."
  !   stop
  ! endif


  iktotal=ik+iktotal
  !  if (MY_RANK==0) write(6,*) 'iktotal=',iktotal,'pb%time=',pb%time
  ! Unpack yt into v, theta
  !  pb%v(pb%rs_nodes) = yt(2::pb%neqs) ! JPA Coulomb

  ! SEISMIC: retrieve the solution for tau in the case of the CNS model, else
  ! retreive the solution for slip velocity
  if (pb%i_rns_law == 3) then
    pb%tau = yt(2::pb%neqs)
    pb%dtau_dt = yt(2::pb%neqs)
  else
    pb%v = yt(2::pb%neqs)
    pb%dtau_dt = 0d0
  endif

  pb%theta = yt(1::pb%neqs)
  ! SEISMIC NOTE/WARNING: I don't know how permanent this temporary solution is,
  ! but in case it gets fixed more permanently, derivs_all.f90 needs adjustment
  if (pb%features%stress_coupling == 1) then           ! Temp solution for normal stress coupling
    pb%sigma = yt(ind_stress_coupling::pb%neqs)
  endif

  if (pb%features%localisation == 1) then
    pb%theta2 = yt(ind_localisation::pb%neqs)
  endif

end subroutine do_bsstep


!=====================================================================
! Update field: slip, tau, potency potency rate, crack,
!
subroutine update_field(pb)

  use output, only : crack_size
  use friction, only : friction_mu, dtheta_dt
  use friction_cns, only : compute_velocity
  use my_mpi, only: max_allproc, is_MPI_parallel
  use diffusion_solver, only : update_PT_final

  type(problem_type), intent(inout) :: pb

  double precision, dimension(pb%mesh%nn) :: P_prev
  integer :: i,ix,iw

  ! SEISMIC: obtain P at the previous time step
  P_prev = 0d0
  if (pb%features%tp == 1) P_prev = pb%tp%P_prev

  ! SEISMIC: in case of the CNS model, re-compute the slip velocity with
  ! the final value of tau, sigma, and porosity. Otherwise, use the standard
  ! rate-and-state expression to calculate tau as a function of velocity
  if (pb%i_rns_law == 3) then
    pb%v = compute_velocity(pb%tau, pb%sigma-P_prev, pb%theta, pb%theta2, pb)
  else
    pb%tau = (pb%sigma-P_prev) * friction_mu(pb%v,pb%theta,pb) + pb%coh
  endif
  ! Update slip
  ! SEISMIC NOTE: slip needs to be calculated after velocity!
  pb%slip = pb%slip + pb%v*pb%dt_did

  ! SEISMIC: if thermal pressurisation is requested, update P and T
  if (pb%features%tp == 1) then
    if (pb%i_rns_law == 3) then
      call dtheta_dt( pb%v, pb%tau, pb%sigma-P_prev, pb%theta, pb%theta2, &
                      pb%dtheta_dt, pb%dtheta2_dt, pb)
    endif
    call update_PT_final(pb%dt_did, pb)
  endif

  ! update potency and potency rate

  if (pb%mesh%dim == 0 .or. pb%mesh%dim == 1) then
    pb%ot%pot = sum(pb%slip*pb%mesh%dx(1))
    pb%ot%pot_rate = sum(pb%v*pb%mesh%dx(1))
  else
    pb%ot%pot=0d0;
    pb%ot%pot_rate=0d0;
    do iw=1,pb%mesh%nw
      do ix=1,pb%mesh%nx
        i=(iw-1)*pb%mesh%nx+ix
        pb%ot%pot = pb%ot%pot + pb%slip(i) * pb%mesh%dx(i) * pb%mesh%dw(iw)
        pb%ot%pot_rate = pb%ot%pot_rate + pb%v(i) * pb%mesh%dx(i) * pb%mesh%dw(iw)
      end do
    end do
  endif
!PG: the crack size only work in serial.
  ! update crack size
  pb%ot%lcold = pb%ot%lcnew
  pb%ot%lcnew = crack_size(pb%slip,pb%mesh%nn)
  pb%ot%llocold = pb%ot%llocnew
  pb%ot%llocnew = crack_size(pb%dtau_dt,pb%mesh%nn)
  ! Output time series at max(v) location
  pb%ot%ivmax = maxloc(pb%v,1)
  if (is_MPI_parallel()) then
    call max_allproc(pb%v(pb%ot%ivmax),pb%vmaxglob)
  else
    pb%vmaxglob = pb%v(pb%ot%ivmax)
  endif

end subroutine update_field

!=====================================================================
! check stop:
!
subroutine check_stop(pb)

  use output, only : time_write
  use my_mpi, only: is_MPI_parallel, is_mpi_master, finalize_mpi

  type(problem_type), intent(inout) :: pb

  double precision, save :: vmax_old = 0d0, vmax_older = 0d0

  if (pb%itstop>0) return

  select case (pb%NSTOP)

   ! STOP if time > tmax
    case (0)
      if (is_mpi_master()) call time_write(pb)
      if (pb%time > pb%tmax) pb%itstop = pb%it

   ! STOP soon after end of slip localization
    case (1)
      if (is_MPI_parallel()) then !JPA WARNING in progress
        print *, 'Stop criterion 1 not implemented yet for MPI'
        call finalize_mpi()
      endif
      if (pb%ot%llocnew > pb%ot%llocold) pb%itstop=pb%it+2*pb%ot%ntout

   ! STOP soon after maximum slip rate
    case (2)
      if (pb%it > 2 .and. vmax_old > vmax_older .and. pb%vmaxglob < vmax_old)  &
        pb%itstop = pb%it+10*pb%ot%ntout
      vmax_older = vmax_old
      vmax_old = pb%vmaxglob

   ! STOP at a slip rate threshold
    case (3)
      if (pb%vmaxglob > pb%tmax) pb%itstop = pb%it    !here tmax is threshhold velocity

    case default
      print *, 'Stop criterion ',pb%NSTOP,' not implemented'
      pb%itstop = pb%it

  end select

end subroutine check_stop


subroutine init_lsoda(pb)

  use problem_class
  use derivs_all
  use ode_lsoda_main, only: dlsoda

  type(problem_type), intent(inout) :: pb
  double precision, dimension(pb%neqs*pb%mesh%nn) :: yt, dydt, yt_scale
  integer :: LRN, LRS
  integer :: ind_stress_coupling, ind_localisation

  ! Basic parameters
  pb%lsoda%neq(1) = pb%neqs * pb%mesh%nn   ! number of equations
  pb%lsoda%tout = pb%dt_try                     ! ignored for itask = {2,5} (one-step)
  pb%lsoda%itol = 1                     ! use scalar rtol and atol
  pb%lsoda%rtol(1) = pb%acc                ! relative tolerance
  pb%lsoda%atol(1) = 0                     ! absolute tolerance
  pb%lsoda%jt = 2                       ! 1 = user supplied full, 2 = internal full jacobian
  pb%lsoda%iopt = 0                     ! no optional parameters
  pb%lsoda%itask = 2                    ! one-step mode, no t_crit

  ! Set-up rwork vector
  LRN = 20 + 16*pb%lsoda%neq(1)            ! size of rwork for non-stiff equations
  LRS = 22 + 9*pb%lsoda%neq(1) + pb%lsoda%neq(1)**2   ! size of rwork for stiff equations
  pb%lsoda%lrw = max(LRN, LRS)                  ! size of rwork for either case
  allocate(pb%lsoda%rwork(pb%lsoda%lrw))        ! allocate rwork
  pb%lsoda%rwork(5) = pb%dt_try         ! set initial step size

  ! Set-up iwork vector
  pb%lsoda%liw = 20+pb%lsoda%neq(1)                   ! size of iwork
  allocate(pb%lsoda%iwork(pb%lsoda%liw))  ! allocate iwork
  pb%lsoda%iwork(6) = 500                 ! max number of internal steps for convergence, default 500

  ! Check if max time step is limited by user
  if (pb%dt_max >  0.d0) then
    pb%lsoda%rwork(5:10) = 0.0            ! reset rwork
    pb%lsoda%iwork(5:10) = 0              ! reset iwork
    pb%lsoda%rwork(6) = pb%dt_max         ! max time step, do not set for infinite (default)
    pb%lsoda%iopt = 1                     ! using optional input
  endif

  pb%lsoda%istate = 1                     ! first call, do sanity checks etc.

  ! Prepare yt for first call
  ind_stress_coupling = 2 + pb%features%stress_coupling
  ind_localisation = ind_stress_coupling + pb%features%localisation

  ! SEISMIC: in the case of the CNS model, solve for tau and not v
  if (pb%i_rns_law == 3) then   ! SEISMIC: CNS model
    yt(2::pb%neqs) = pb%tau
  else  ! SEISMIC: not CNS model (i.e. rate-and-state)
    yt(2::pb%neqs) = pb%v
  endif
  yt(1::pb%neqs) = pb%theta
  ! SEISMIC NOTE/WARNING: I don't know how permanent this temporary solution is,
  ! but in case it gets fixed more permanently, derivs_all.f90 needs adjustment
  if (pb%features%stress_coupling == 1) then           ! Temp solution for normal stress coupling
    yt(ind_stress_coupling::pb%neqs) = pb%sigma
  endif
  if (pb%features%localisation == 1) then
    yt(ind_localisation::pb%neqs) = pb%theta2
  endif

  call dlsoda(  derivs_lsoda, pb%lsoda%neq, yt, pb%time, pb%lsoda%tout, &
                pb%lsoda%itol, pb%lsoda%rtol, pb%lsoda%atol, pb%lsoda%itask, &
                pb%lsoda%istate, pb%lsoda%iopt, pb%lsoda%rwork, pb%lsoda%lrw, &
                pb%lsoda%iwork, pb%lsoda%liw, jac_lsoda, pb%lsoda%jt)

  ! istate output values:
  !  1: nothing was done (t = tout)
  !  2: great success!
  ! -1: excessive amount of work done
  ! -2: too much accuracy requested for machine precision
  ! -3: illegal output detected before integrating
  ! -4: repeated error test failures
  ! -5: repeated convergence test failures
  ! -6: EWT(i) become zero (error became zero)
  ! -7: length of RWORK/IWORK was too small to proceed

  if (pb%lsoda%istate /= 2) then
    write (6,*) "Initiation of LSODA solver failed!"
    write (6,*) "Current value of istate: ", pb%lsoda%istate
    write (6,*) "Will now terminate..."
    stop
  endif
  ! We've made the first call
  ! No changes in parameters except tout and itask are allowed
  pb%lsoda%istate = 2

end subroutine init_lsoda


subroutine init_rk45(pb)

  use problem_class
  use derivs_all
  use ode_rk45, only: rkf45_d

  type(problem_type), intent(inout) :: pb
  double precision, dimension(pb%neqs*pb%mesh%nn) :: yt, dydt
  integer :: nwork
  integer :: ind_stress_coupling, ind_localisation

  write (6,*) "Initialising RK45 solver"

  nwork = 3 + 6*pb%neqs*pb%mesh%nn
  pb%rk45%iflag = -1
  allocate(pb%rk45%work(nwork))
  allocate(pb%rk45%iwork(5))

  ! Prepare yt for first call
  ind_stress_coupling = 2 + pb%features%stress_coupling
  ind_localisation = ind_stress_coupling + pb%features%localisation

  ! SEISMIC: in the case of the CNS model, solve for tau and not v
  if (pb%i_rns_law == 3) then   ! SEISMIC: CNS model
    yt(2::pb%neqs) = pb%tau
  else  ! SEISMIC: not CNS model (i.e. rate-and-state)
    yt(2::pb%neqs) = pb%v
  endif
  yt(1::pb%neqs) = pb%theta
  ! SEISMIC NOTE/WARNING: I don't know how permanent this temporary solution is,
  ! but in case it gets fixed more permanently, derivs_all.f90 needs adjustment
  if (pb%features%stress_coupling == 1) then           ! Temp solution for normal stress coupling
    yt(ind_stress_coupling::pb%neqs) = pb%sigma
  endif
  if (pb%features%localisation == 1) then
    yt(ind_localisation::pb%neqs) = pb%theta2
  endif

  call rkf45_d( derivs_rk45, pb%neqs*pb%mesh%nn, yt, pb%time, pb%time, &
                pb%acc, 0d0, pb%rk45%iflag, pb%rk45%work, pb%rk45%iwork)

  select case (pb%rk45%iflag)
  case (3)
    write (6,*) "RK45 error [3]: relative error tolerance too small"
    stop
  case (4)
    write (6,*) "RK45 warning [4]: integration took more than 3000 derivative evaluations"
  case (5)
    write (6,*) "RK45 error [5]: solution vanished, relative error test is not possible"
    stop
  case (6)
    write (6,*) "RK45 error [6]: requested accuracy could not be achieved"
    stop
  case (8)
    write (6,*) "RK45 error [8]: invalid input parameters"
    stop
  end select

write (6,*) "Finished initialising RK45 solver"

end subroutine init_rk45


end module solver
