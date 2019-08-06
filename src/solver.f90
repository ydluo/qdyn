! Solve_master

module solver

  use problem_class, only : problem_type
  use utils, only : pack, unpack

  implicit none
  private

  integer(kind=8), save :: iktotal

  public  :: solve, init_rk45

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
    call check_stop(pb)   ! here itstop will change
!--------Output onestep to screen and ox file(snap_shot)
! if(mod(pb%it-1,pb%ot%ntout) == 0 .or. pb%it == pb%itstop) then
    if(mod(pb%it,pb%ox%ntout) == 0) then
!      if (is_mpi_master()) write(6,*) 'it:',pb%it,'iktotal=',iktotal,'pb%time=',pb%time
      call screen_write(pb)
    endif

    if (pb%it /= pb%itstop) then
      call ox_write(pb)
    endif
! if (is_mpi_master()) call ox_write(pb)
    if (mod(pb%it, pb%ot%ntout) == 0 .and. pb%it /= pb%itstop) then
      call ot_write(pb)
    endif
  enddo

  ! Write data for last step
  call screen_write(pb)
  call ox_write(pb)
  call ot_write(pb)

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
  use ode_rk45, only: rkf45_d
  use ode_rk45_2, only: rkf45_d2
  use output, only : screen_write, ox_write, ot_write
  use constants, only : SOLVER_TYPE
  use diffusion_solver, only : update_PT_final

  type(problem_type), intent(inout) :: pb

  double precision, dimension(pb%neqs*pb%mesh%nn) :: yt, dydt, yt_scale
  double precision, dimension(pb%mesh%nn) :: main_var
  double precision :: t_out
  integer :: ik, neqs

  neqs = pb%neqs * pb%mesh%nn

  ! SEISMIC: in the case of the CNS model, solve for tau and not v
  if (pb%i_rns_law == 3) then   ! SEISMIC: CNS model
    main_var = pb%tau
  else  ! SEISMIC: not CNS model (i.e. rate-and-state)
    main_var = pb%v
  endif

  call pack(yt, pb%theta, main_var, pb%sigma, pb%theta2, pb)

  ! SEISMIC: user-defined switch to use either (1) the Bulirsch-Stoer method, or
  ! the (2) Runge-Kutta-Fehlberg method
  if (SOLVER_TYPE == 0) then
    ! Default value of SOLVER_TYPE has not been altered
    write (6,*) "The default solver type (0) has not been altered, and no solver was picked"
    write( 6,*) "Check the input script and define a solver type > 0"
    stop

  elseif (SOLVER_TYPE == 1) then
    ! Use Bulirsch-Stoer method

    ! this update of derivatives is only needed to set up the scaling (yt_scale)
    call derivs(pb%time,yt,dydt,pb)
    yt_scale=dabs(yt)+dabs(pb%dt_try*dydt)
    ! One step
    call bsstep(yt,dydt,neqs,pb%time,pb%dt_try,pb%acc,yt_scale,pb%dt_did,pb%dt_next,pb,ik)

    ! SEISMIC NOTE: what is happening here?
    if (pb%dt_max >  0.d0) then
      pb%dt_try = min(pb%dt_next,pb%dt_max)
    else
      pb%dt_try = pb%dt_next
    endif

  elseif (SOLVER_TYPE == 2) then
    ! Set-up Runge-Kutta solver

    pb%rk45%iflag = -2 ! Reset to one-step mode each call
    pb%t_prev = pb%time

    ! Call Runge-Kutta solver routine
    call rkf45_d( derivs_rk45, neqs, yt, pb%time, pb%tmax, &
                  pb%acc, 0d0, pb%rk45%iflag, pb%rk45%work, pb%rk45%iwork)

    ! Set time step
    pb%dt_did = pb%time - pb%t_prev

    ! Basic error checking. See description of rkf45_d in ode_rk45.f90 for details
    select case (pb%rk45%iflag)
    case (3)
      write (6,*) "RK45 error [3]: relative error tolerance too small"
      stop
    case (4)
      ! write (6,*) "RK45 warning [4]: integration took more than 3000 derivative evaluations"
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

  elseif (SOLVER_TYPE == 3) then
    ! Set-up Runge-Kutta solver
    pb%t_prev = pb%time
    ! Call Runge-Kutta solver routine
    call rkf45_d2(derivs, yt, pb%time, pb%dt_max, pb%acc, 0d0, pb)
    ! Set time step
    pb%dt_did = pb%time - pb%t_prev
  else
    ! Unknown solver type
    write (6,*) "Solver type", SOLVER_TYPE, "not recognised"
    stop
  endif

  iktotal=ik+iktotal

  call unpack(yt, pb%theta, main_var, pb%sigma, pb%theta2, pb)
  if (pb%features%tp == 1) call update_PT_final(pb%dt_did, pb)

  ! SEISMIC: retrieve the solution for tau in the case of the CNS model, else
  ! retreive the solution for slip velocity
  if (pb%i_rns_law == 3) then
    pb%tau = main_var
  else
    pb%v = main_var
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

  double precision, dimension(pb%mesh%nn) :: P
  integer :: i,ix,iw

  ! SEISMIC: obtain P at the previous time step
  P = 0d0
  if (pb%features%tp == 1) P = pb%tp%P

  ! SEISMIC: in case of the CNS model, re-compute the slip velocity with
  ! the final value of tau, sigma, and porosity. Otherwise, use the standard
  ! rate-and-state expression to calculate tau as a function of velocity
  if (pb%i_rns_law == 3) then
    pb%v = compute_velocity(pb%tau, pb%sigma-P, pb%theta, pb%theta2, pb)
  else
    pb%tau = (pb%sigma-P) * friction_mu(pb%v,pb%theta,pb) + pb%coh
  endif
  ! Update slip
  ! SEISMIC NOTE: slip needs to be calculated after velocity!
  ! NOTE 2: include slip in solver routine to get higher order accuracy
  pb%slip = pb%slip + pb%v*pb%dt_did

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
        pb%ot%pot = pb%ot%pot + pb%slip(i) * pb%mesh%dx(1) * pb%mesh%dw(iw)
        pb%ot%pot_rate = pb%ot%pot_rate + pb%v(i) * pb%mesh%dx(1) * pb%mesh%dw(iw)
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
      if (pb%time >= pb%tmax) pb%itstop = pb%it

   ! STOP soon after end of slip localization
    case (1)
      if (is_MPI_parallel()) then !JPA WARNING in progress
        print *, 'Stop criterion 1 not implemented yet for MPI'
        call finalize_mpi()
      endif
      if (pb%ot%llocnew > pb%ot%llocold) pb%itstop=pb%it+2*pb%ox%ntout

   ! STOP soon after maximum slip rate
    case (2)
      if (pb%it > 2 .and. vmax_old > vmax_older .and. pb%vmaxglob < vmax_old)  &
        pb%itstop = pb%it+10*pb%ox%ntout
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


subroutine init_rk45(pb)

  use problem_class
  use derivs_all
  use ode_rk45, only: rkf45_d

  type(problem_type), intent(inout) :: pb
  double precision, dimension(pb%neqs*pb%mesh%nn) :: yt, dydt
  double precision, dimension(pb%mesh%nn) :: main_var
  integer :: nwork

  write (6,*) "Initialising RK45 solver"

  nwork = 3 + 6*pb%neqs*pb%mesh%nn
  pb%rk45%iflag = -1
  allocate(pb%rk45%work(nwork))
  allocate(pb%rk45%iwork(5))

  if (pb%i_rns_law == 3) then   ! SEISMIC: CNS model
    main_var = pb%tau
  else  ! SEISMIC: not CNS model (i.e. rate-and-state)
    main_var = pb%v
  endif

  call pack(yt, pb%theta, main_var, pb%sigma, pb%theta2, pb)

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
