module output

! OUTPUT: This module manages outputs

  implicit none
  private

 ! id in output file names "fort.id"
  integer, parameter :: FID_OT = 18, FID_VMAX = 22, FID_IOT_0 = 10000

  public :: screen_init, ot_read_stations, ot_init, ox_init, &
            screen_write, ot_write, ox_write,  &
            time_write, crack_size

contains

!=====================================================================
!output initilized field to screen
subroutine screen_init(pb)

  use problem_class
  use constants, only : PI

  type (problem_type), intent(inout) :: pb

  double precision :: K

  ! SEISMIC: skip calculating critical stiffness for CNS model
  ! Is not very useful
  if (pb%i_rns_law /= 3) then

  if (pb%ot%ic<1) return

    write(6,*) 'Values at selected point of the fault:'
    K = pb%mesh%Lfault
    if (pb%mesh%dim == 1) then
      if (.not. pb%kernel%k2f%finite) K = pb%mesh%W
    endif
    K = PI*pb%smu/K
    if (pb%mesh%dim < 2) then
      write(6,*) 'K/Kc = ',K/(pb%sigma(pb%ot%ic)*(pb%b(pb%ot%ic)-pb%a(pb%ot%ic))/pb%dc(pb%ot%ic))
      write(6,*) 'K/Kb = ',K/(pb%sigma(pb%ot%ic)*pb%b(pb%ot%ic)/pb%dc(pb%ot%ic))
    end if

  endif

  write(6,*)
  write(6,*) '    it,  dt (secs), time (yrs), v_max (m/s), sigma_max (MPa)'

end subroutine screen_init



!=====================================================================
!output one step to screen
subroutine screen_write(pb)

  use constants, only : YEAR
  use problem_class
  use my_mpi, only : is_MPI_parallel, is_mpi_master, max_allproc

  type (problem_type), intent(in) :: pb
  double precision :: sigma_max, sigma_max_glob

  sigma_max = maxval(pb%sigma)

  if (is_MPI_parallel()) then
    call max_allproc(sigma_max,sigma_max_glob)
    if (is_mpi_master()) write(6,'(i7,x,4(e11.3,x),i5)') pb%it, pb%dt_did, pb%time/YEAR,&
                              pb%vmaxglob, sigma_max_glob/1.0D6
  else
    write(6,'(i7,x,4(e11.3,x),i5)') pb%it, pb%dt_did, pb%time/YEAR,    &
                            pb%vmaxglob, sigma_max/1.0D6
  endif

end subroutine screen_write


!=====================================================================
! write time of every step
subroutine time_write(pb)

  use problem_class
  use my_mpi, only: is_mpi_master
  type (problem_type), intent(inout) :: pb

  if (is_mpi_master()) write(121,*) pb%time, pb%dt_did

end subroutine time_write


!=====================================================================
subroutine ot_read_stations(ot)

  use problem_class, only : ot_type

  type (ot_type), intent(inout) :: ot

  integer :: nsta,ista

  open(unit=200,file='stations.dat',action='read')
  read(200,*) nsta
  allocate(ot%xsta(nsta),ot%ysta(nsta),ot%zsta(nsta))
  do ista=1,nsta
    read(200,*) ot%xsta(ista), ot%ysta(ista), ot%zsta(ista)
  enddo
  close(200)

end subroutine ot_read_stations

!=====================================================================
! write ot file header
subroutine ot_init(pb)

  use problem_class
  use constants, only: OUT_MASTER, BIN_OUTPUT
  use my_mpi, only : is_MPI_parallel, is_mpi_master, my_mpi_tag
  use mesh, only : mesh_get_size

  type (problem_type), intent(inout) :: pb
  integer :: i,j,n,nsta,ista,ik
  double precision :: dmin2, d2

  n = mesh_get_size(pb%mesh)
  pb%ot%lcnew = dble(n)
  pb%ot%llocnew = dble(n)

  if (is_MPI_parallel()) then

   ! find stations
   !PG, temporal fix, only works for one station.
    ! pb%ot%ic=0 !Setting all to zero
    ! dmin2 = 0.01d0*minval(pb%mesh%dx) ! distance tolerance = 1% grid size
    ! dmin2 = dmin2*dmin2
    ! nsta = 1 !NOTE: currently only one station implemented
    !
    ! do ista=1,nsta
    !   do ik=1,n
    !     d2 = (pb%mesh%x(ik)-pb%ot%xsta(ista))**2 &
    !        + (pb%mesh%y(ik)-pb%ot%ysta(ista))**2 &
    !        + (pb%mesh%z(ik)-pb%ot%zsta(ista))**2
    !     if (d2 < dmin2) then
    !       pb%ot%ic=ik
    !       write(6,*) 'Processor: ',my_mpi_tag(),', station ',ista, &
    !                  ' found, distance mismatch = ',d2
    !       exit
    !     endif
    !   enddo
    ! enddo

    if (OUT_MASTER .and. is_mpi_master() ) then
      pb%ot%unit = FID_OT
      write(pb%ot%unit,'(a)')'# macroscopic values:'
      write(pb%ot%unit,'(a)')'# 1=t'
      write(pb%ot%unit,'(a)')'# values at selected point:'
      write(pb%ot%unit,'(a)')'# 2=V, 3=theta, 4=V*theta/dc, 5=tau, 6=slip'
      if (pb%features%tp == 1) write(pb%ot%unit,'(a)')'# 7=P, 8=T'
      close(pb%ot%unit)
    endif

   !JPA WARNING VMAX and IOT outputs not implemented yet in parallel
   ! TODO: why are these not implemented?

  else

    pb%ot%unit = FID_OT
    if (.not.BIN_OUTPUT) then
      write(pb%ot%unit,'(a)')'# macroscopic values:'
      write(pb%ot%unit,'(a)')'# 1=t,2=loc_size,3=crack_size,4=potcy,5=pot_rate'
      write(pb%ot%unit,'(a)')'# values at selected point:'
      write(pb%ot%unit,'(a)')'# 6=V, 7=theta, 8=V*theta/dc, 9=tau, 10=slip'
      write(pb%ot%unit,'(a)')'# values at max(V) location:'
      write(pb%ot%unit,'(a)')'# 11=x, 12=V, 13=theta, 14=omeg, 15=tau, 16=slip, 17=sigma'
      if (pb%features%tp == 1) then
        write(pb%ot%unit,'(a)')'# 18=P, 19=T (selected point)'
        write(pb%ot%unit,'(a)')'# 20=P, 21=T (at max(V))'
      endif
    else
      open(pb%ot%unit,form='unformatted',access='stream')
    endif

    write(FID_VMAX,'(a)')'# Seismicity record:'
    write(FID_VMAX,'(a)')'# 1=loc, 2=t, 3=v'

    j = FID_IOT_0
    do i=1,n
      if (pb%ot%iot(i) == 1) then
        j = j+1
        write(j,'(a,i10)') '# nx= ', i
        write(j,'(a)') '# 1=t, 2=V, 3=theta, 4=tau, 5=slip, 6=sigma'
      endif
    enddo

  endif

  allocate ( pb%ot%v_pre(n), pb%ot%v_pre2(n) )
  pb%ot%v_pre = 0.d0
  pb%ot%v_pre2 = 0.d0

end subroutine ot_init


!=====================================================================
! write ox file header
subroutine ox_init(pb)

  use problem_class
  use my_mpi, only : is_MPI_parallel, is_mpi_master
  use mesh, only : mesh_get_size

  type (problem_type), intent(inout) :: pb

  integer :: i, n

  ! Number of mesh elements
  n = pb%mesh%nn
  ! Number of ox elements
  pb%ox%count = ceiling(n / float(pb%ox%nxout))
  ! Initial potency
  pb%ox%pot_pre = 0.d0

  ! ---------------------------------------------------------------------------
  ! Prepare data structure and headers for ox, ox_dyn, and dyn output

  ! Number of ox output quantities
  pb%ox%nox = 10
  pb%ox%nrup = 2
  ! If thermal pressurisation is requested, add 2 more
  if (pb%features%tp == 1) then
    pb%ox%nox = pb%ox%nox + 2
  endif
  ! Allocate space in array of pointers
  allocate(pb%ox%objects_dyn(pb%ox%nox))
  allocate(pb%ox%objects_glob(pb%ox%nox))
  allocate(pb%ox%objects_loc(pb%ox%nox))
  allocate(pb%ox%objects_ox(pb%ox%nox))
  allocate(pb%ox%objects_rup(pb%ox%nrup))
  allocate(pb%ox%fmt(pb%ox%nox))
  ! Default output format
  pb%ox%fmt = "(e15.7)"
  ! Time needs higher precision
  pb%ox%fmt(4) = "(e24.14)"

  ! Initialise MPI gather
  call init_pb_global(pb)

  ! Set pointers to specific output quantities

  ! Proc local variables
  ! Spatial coordinates
  pb%ox%objects_loc(1)%p => pb%mesh%x
  pb%ox%objects_loc(2)%p => pb%mesh%y
  pb%ox%objects_loc(3)%p => pb%mesh%z
  ! Time array
  pb%ox%objects_loc(4)%p => pb%mesh%time
  ! Mechanical quantities
  pb%ox%objects_loc(5)%p => pb%v
  pb%ox%objects_loc(6)%p => pb%theta
  pb%ox%objects_loc(7)%p => pb%tau
  pb%ox%objects_loc(8)%p => pb%dtau_dt
  pb%ox%objects_loc(9)%p => pb%slip
  pb%ox%objects_loc(10)%p => pb%sigma

  ! Global objects (required for ox_dyn and QSB output)
  ! Spatial coordinates
  pb%ox%objects_glob(1)%p => pb%mesh%xglob
  pb%ox%objects_glob(2)%p => pb%mesh%yglob
  pb%ox%objects_glob(3)%p => pb%mesh%zglob
  ! Time array
  pb%ox%objects_glob(4)%p => pb%mesh%time
  ! Mechanical quantities
  pb%ox%objects_glob(5)%p => pb%v_glob
  pb%ox%objects_glob(6)%p => pb%theta_glob
  pb%ox%objects_glob(7)%p => pb%tau_glob
  pb%ox%objects_glob(8)%p => pb%dtau_dt_glob
  pb%ox%objects_glob(9)%p => pb%slip_glob
  pb%ox%objects_glob(10)%p => pb%sigma_glob

  ! If thermal pressurisation is requested, add P and T
  if (pb%features%tp == 1) then
    ! Local
    pb%ox%objects_loc(11)%p => pb%tp%P
    pb%ox%objects_loc(12)%p => pb%tp%T
    ! Global
    pb%ox%objects_glob(11)%p => pb%P_glob
    pb%ox%objects_glob(12)%p => pb%T_glob
  endif

  ! Case 1: serial execution. All quantities are local
  if (.not. is_MPI_parallel()) then
    ! Assign local objects to ox and dynamic objects
    do i=1,pb%ox%nox
      pb%ox%objects_dyn(i)%p => pb%ox%objects_loc(i)%p
      pb%ox%objects_ox(i)%p => pb%ox%objects_loc(i)%p
    enddo
    ! Assign local tau/v to rupture object
    pb%ox%objects_rup(1)%p => pb%tau
    pb%ox%objects_rup(2)%p => pb%v

  ! Case 2: parallel execution with output master. All quantities are global
  elseif (is_MPI_parallel() .and. is_MPI_master()) then
    ! Assign global objects to ox and dynamic objects
    do i=1,pb%ox%nox
      pb%ox%objects_dyn(i)%p => pb%ox%objects_glob(i)%p
      pb%ox%objects_ox(i)%p => pb%ox%objects_glob(i)%p
    enddo
    ! Assign global tau/v to rupture object
    pb%ox%objects_rup(1)%p => pb%tau_glob
    pb%ox%objects_rup(2)%p => pb%v_glob

  endif

  ! Define headers
  pb%ox%header = '# x y z t v theta dtau tau_dot slip sigma'
  if (pb%features%tp == 1) then
    pb%ox%header = '# x y z t v theta dtau tau_dot slip sigma P T'
  endif

  ! Define output units

  ! Standard ox output unit
  pb%ox%unit = 19

  ! If ox_dyn is requested
  if (pb%ox%i_ox_dyn == 1) then
    ! ox_dyn unit
    pb%ox%ox_dyn_unit = 20000
    ! Allocate rupture quantities
    allocate ( pb%ox%t_rup(n), pb%ox%tau_max(n), pb%ox%v_max(n), pb%ox%t_vmax(n) )
  endif

  if (pb%DYN_FLAG == 1) then
    ! Define output units
    pb%ox%QSB_unit_pre = 100
    pb%ox%QSB_unit_post = 101
  endif

  ! If dynamic output is requested (ox_dyn or QSB)
  if ((pb%ox%i_ox_dyn == 1) .or. (pb%DYN_FLAG == 1)) then
    ! Number of dynamic events (start at 0)
    pb%ox%dyn_count = 0
    ! Dynamic event status (0 = no, 1 = yes)
    pb%ox%dyn_stat = 0
  endif

  ! End preparing data structures/headers
  ! ---------------------------------------------------------------------------

end subroutine ox_init

!=====================================================================
! Export timeseries
subroutine ot_write(pb)

  use problem_class
  use constants, only : OCTAVE_OUTPUT, BIN_OUTPUT
  use my_mpi, only : is_MPI_parallel, my_mpi_tag

  type (problem_type), intent(inout) :: pb
  integer :: i,j,ios
  character(30) :: ot_fmt

  if (is_MPI_parallel()) then
   ! if "ic" station is in this processor
    if (pb%ot%ic>0) then

      open(pb%ot%unit,access='APPEND',status='old',iostat=ios)
      if (ios>0) stop 'Fatal error: ot_write: Error opening a fort.18 file'
     !JPA add test for the first time we try to open this file but it does not exist yet
     !JPA add test to prevent appending data to a file from a previous simulation
      if (OCTAVE_OUTPUT) then
        ot_fmt = '(g0.16,5(",",g0.6))'
        if (pb%features%tp == 1) ot_fmt = '(g0.16,7(",",g0.6))'
      else
        ot_fmt = '(e24.16,5e14.6)'
        if (pb%features%tp == 1) ot_fmt = '(e24.16,7e14.6)'
      endif

      ! If thermal pressurisation is requested, output P and T
      if (pb%features%tp == 1) then
        write(pb%ot%unit,ot_fmt) pb%time, pb%v(pb%ot%ic), pb%theta(pb%ot%ic), &
            pb%v(pb%ot%ic)*pb%theta(pb%ot%ic)/pb%dc(pb%ot%ic), &
            pb%tau(pb%ot%ic), pb%slip(pb%ot%ic), &
            pb%tp%P(pb%ot%ic), pb%tp%T(pb%ot%ic)
      else
        write(pb%ot%unit,ot_fmt) pb%time, pb%v(pb%ot%ic), pb%theta(pb%ot%ic), &
            pb%v(pb%ot%ic)*pb%theta(pb%ot%ic)/pb%dc(pb%ot%ic), &
            pb%tau(pb%ot%ic), pb%slip(pb%ot%ic)
      endif
      close(pb%ot%unit)
    endif
   !JPA warning: ivmax outputs not implemented in parallel yet

  else

    if (.not.BIN_OUTPUT) then

      if (OCTAVE_OUTPUT) then
        ! for Octave: comma as field delimiter and no spaces
        ot_fmt = '(g0.16,16(",",g0.6))'
        if (pb%features%tp == 1) ot_fmt = '(g0.16,20(",",g0.6))'
      else
        ot_fmt = '(e24.16,16e14.6)'
        if (pb%features%tp == 1) ot_fmt = '(e24.16,20e14.6)'
      endif
      if (pb%features%tp == 1) then
        write(pb%ot%unit,ot_fmt) pb%time, pb%ot%llocnew*pb%mesh%dx(1),  &
          pb%ot%lcnew*pb%mesh%dx(1), pb%ot%pot, pb%ot%pot_rate,    &
          pb%v(pb%ot%ic), pb%theta(pb%ot%ic),  &
          pb%v(pb%ot%ic)*pb%theta(pb%ot%ic)/pb%dc(pb%ot%ic), &
          pb%tau(pb%ot%ic), pb%slip(pb%ot%ic),    &
          ! for ivmax
          pb%mesh%x(pb%ot%ivmax), pb%v(pb%ot%ivmax), pb%theta(pb%ot%ivmax),   &
          pb%v(pb%ot%ivmax)*pb%theta(pb%ot%ivmax)/pb%dc(pb%ot%ivmax),    &
          pb%tau(pb%ot%ivmax), pb%slip(pb%ot%ivmax), &
          pb%sigma(pb%ot%ivmax)-pb%tp%P(pb%ot%ivmax), &
          ! write P, T
          pb%tp%P(pb%ot%ic), pb%tp%T(pb%ot%ic), pb%tp%P(pb%ot%ivmax), pb%tp%T(pb%ot%ivmax)
      else
        write(pb%ot%unit,ot_fmt) pb%time, pb%ot%llocnew*pb%mesh%dx(1),  &
          pb%ot%lcnew*pb%mesh%dx(1), pb%ot%pot, pb%ot%pot_rate,    &
          pb%v(pb%ot%ic), pb%theta(pb%ot%ic),  &
          pb%v(pb%ot%ic)*pb%theta(pb%ot%ic)/pb%dc(pb%ot%ic), &
          pb%tau(pb%ot%ic), pb%slip(pb%ot%ic),    &
          ! for ivmax
          pb%mesh%x(pb%ot%ivmax), pb%v(pb%ot%ivmax), pb%theta(pb%ot%ivmax),   &
          pb%v(pb%ot%ivmax)*pb%theta(pb%ot%ivmax)/pb%dc(pb%ot%ivmax),    &
          pb%tau(pb%ot%ivmax), pb%slip(pb%ot%ivmax), pb%sigma(pb%ot%ivmax)
      endif

    else
      ! macroscopic values:
      !   1=t,2=loc_size,3=crack_size,4=potcy,5=pot_rate
      ! values at selected point:
      !   6=V, 7=theta, 8=V*theta/dc, 9=tau, 10=slip
      ! values at max(V) location:
      !   11=x, 12=V, 13=theta, 14=omeg, 15=tau, 16=slip, 17=sigma
      ! If thermal pressurisation: 18=P, 19=T, 20=P_max, 21=T_max

      if (pb%features%tp == 1) then
        write(pb%ot%unit) pb%time, pb%ot%llocnew*pb%mesh%dx(1),  &
          pb%ot%lcnew*pb%mesh%dx(1), pb%ot%pot, pb%ot%pot_rate,    &
          pb%v(pb%ot%ic), pb%theta(pb%ot%ic),  &
          pb%v(pb%ot%ic)*pb%theta(pb%ot%ic)/pb%dc(pb%ot%ic), &
          pb%tau(pb%ot%ic), pb%slip(pb%ot%ic),    &
          ! for ivmax
          pb%mesh%x(pb%ot%ivmax), pb%v(pb%ot%ivmax), pb%theta(pb%ot%ivmax),   &
          pb%v(pb%ot%ivmax)*pb%theta(pb%ot%ivmax)/pb%dc(pb%ot%ivmax),    &
          pb%tau(pb%ot%ivmax), pb%slip(pb%ot%ivmax), &
          pb%sigma(pb%ot%ivmax)-pb%tp%P(pb%ot%ivmax), &
          ! write P, T
          pb%tp%P(pb%ot%ic), pb%tp%T(pb%ot%ic), pb%tp%P(pb%ot%ivmax), pb%tp%T(pb%ot%ivmax)
      else
        write(pb%ot%unit) pb%time, pb%ot%llocnew*pb%mesh%dx(1),  &
          pb%ot%lcnew*pb%mesh%dx(1), pb%ot%pot, pb%ot%pot_rate,    &
          pb%v(pb%ot%ic), pb%theta(pb%ot%ic),  &
          pb%v(pb%ot%ic)*pb%theta(pb%ot%ic)/pb%dc(pb%ot%ic), &
          pb%tau(pb%ot%ic), pb%slip(pb%ot%ic),    &
          ! for ivmax
          pb%mesh%x(pb%ot%ivmax), pb%v(pb%ot%ivmax), pb%theta(pb%ot%ivmax),   &
          pb%v(pb%ot%ivmax)*pb%theta(pb%ot%ivmax)/pb%dc(pb%ot%ivmax),    &
          pb%tau(pb%ot%ivmax), pb%slip(pb%ot%ivmax), pb%sigma(pb%ot%ivmax)
      endif
      if (pb%it+1 == pb%itstop) close(pb%ot%unit)
    endif

  endif

  if (is_MPI_parallel()) return
 !JPA warning: the outputs below are not yet implemented in parallel

 ! output slip velocity maxima at selected nodes
  do i=1,pb%mesh%nn
    if ((pb%ot%iasp(i) == 1) .and. (pb%ot%v_pre(i) >= pb%ot%v_th) .and.      &
        (pb%v(i) < pb%ot%v_pre(i)) .and. (pb%ot%v_pre(i) >= pb%ot%v_pre2(i))) then
      write(FID_VMAX,'(i10,2e24.16)') i, pb%time, pb%ot%v_pre(i)
    endif
  enddo
  pb%ot%v_pre2=pb%ot%v_pre
  pb%ot%v_pre=pb%v

 ! output time series at selected nodes
  j = FID_IOT_0
  do i=1,pb%mesh%nn
    if (pb%ot%iot(i) == 1) then
      j = j+1
      write(j,'(e24.16,5e14.6)') pb%time, pb%v(i), pb%theta(i), pb%tau(i), pb%slip(i), pb%sigma(i)
    endif
  enddo

end subroutine ot_write

!=====================================================================
! Export snapshots
subroutine ox_write(pb)

  use problem_class
  use my_mpi, only: is_MPI_parallel, is_mpi_master

  type (problem_type), intent(inout) :: pb

  integer :: ixout, nox, nxout_dyn, unit
  logical ::  call_gather, close_unit, dynamic, falling_edge, last_call, &
              MPI_master, rising_edge, write_ox, write_ox_dyn, write_QSB
  double precision, dimension(pb%mesh%nnglob) :: tau, v

  nox = pb%ox%nox
  nxout_dyn = pb%ox%nxout_dyn

  !---------------------------------------------------------------------------
  ! Perform series of checks to see which operations need to be executed

  ! Is this the last call (last simulation step)?
  ! TODO: check if this is correct, maybe pb%itstop - 1?
  last_call = (pb%it == pb%itstop)

  ! Is this proc MPI master? (default .true. for serial)
  MPI_master = is_MPI_master()

  ! Check if we're crossing dynamic thresholds from below (rising edge)
  ! or from above (falling edge)
  rising_edge = ((pb%ox%dyn_stat == 0) .and. &
                (pb%v(pb%ot%ivmax) >= pb%DYN_th_on))
  falling_edge =  ((pb%ox%dyn_stat == 1) .and. &
                  (pb%v(pb%ot%ivmax) <= pb%DYN_th_off))
  ! If we're crossing a dynamic threshold, update global dynamic state
  if (rising_edge .and. MPI_master) then
    ! Set global dynamic state to 1
    pb%ox%dyn_stat = 1
  endif
  ! We're dynamic as long as dyn_stat == 1
  dynamic = (pb%ox%dyn_stat == 1)

  ! Should this proc write ox, ox_dyn, or QSB output?
  write_ox = (mod(pb%it, pb%ox%ntout) == 0 .or. last_call) .and. MPI_master
  write_ox_dyn = ((pb%ox%i_ox_dyn == 1) .and. dynamic) .and. MPI_master
  write_QSB = ((pb%DYN_FLAG == 1) .and. dynamic) .and. MPI_master

  ! Call an MPI gather when:
  !  1. Output is requested (either regular ox or other flavour)
  !  2. AND parallel execution
  !  3. AND this proc is MPI master
  call_gather = (write_ox .or. write_ox_dyn .or. write_QSB) .and. &
                is_MPI_parallel() .and. MPI_master

  ! Does the output unit need to be closed?
  close_unit = last_call

  ! End checks
  !---------------------------------------------------------------------------
  ! Perform an MPI global gather

  if (call_gather) then
    call pb_global(pb)
  endif

  ! End MPI global gather
  !---------------------------------------------------------------------------
  ! Write regular ox output

  if (write_ox) then

    ! Write ox data
    call write_ox_lines(pb%ox%unit, pb%ox%fmt, pb%ox%objects_ox, pb)

    ! Check if unit needs to be closed
    if (close_unit) then
      close(pb%ox%unit)
    endif

  endif

  ! End regular ox output
  !---------------------------------------------------------------------------
  ! Write ox_dyn output

  if (write_ox_dyn) then

    ! If we cross the dynamic threshold from below
    if (rising_edge) then

      ! Collect output quantities
      ! TODO: check slice is same as loop?
      pb%tau_max_glob = pb%tau_glob(::nxout_dyn)
      pb%v_max_glob = pb%v_glob(::nxout_dyn)
      pb%t_rup_glob = pb%time
      pb%t_vmax_glob = pb%time
      ! Define unit for rising edge ox output
      unit = pb%ox%ox_dyn_unit + 3 * pb%ox%dyn_count + 1
      ! Write ox data
      call write_ox_lines(unit, pb%ox%fmt, pb%ox%objects_dyn, pb)
      ! Close output unit
      close(unit)

    ! If we cross the dynamic threshold from above
    elseif (falling_edge) then

      ! First: falling edge ox output
      ! Define unit for falling edge ox output
      unit = pb%ox%ox_dyn_unit + 3 * pb%ox%dyn_count + 2
      ! Write ox data
      call write_ox_lines(unit, pb%ox%fmt, pb%ox%objects_dyn, pb)
      ! Close output unit
      close(unit)

      ! Second: max rupture stats
      ! TODO: move this to write_ox_lines subroutine?
      ! Define unit for max stats output
      unit = pb%ox%ox_dyn_unit + 3 * pb%ox%dyn_count + 3
      write(unit,'(a)') '# x y z t_rup tau_max t_vmax vmax'
      do ixout=1,pb%mesh%nnglob,nxout_dyn
        write(unit, '(3e15.7,4e28.20)') &
          pb%mesh%xglob(ixout), pb%mesh%yglob(ixout), pb%mesh%zglob(ixout), &
          pb%t_rup_glob(ixout), pb%tau_max_glob(ixout), pb%t_vmax_glob(ixout), &
          pb%v_max_glob(ixout)
      enddo
      ! Close output unit
      close(unit)

    ! No threshold crossing, but still dynamic
    else
      tau = pb%ox%objects_rup(1)%p
      v = pb%ox%objects_rup(2)%p
      ! Update max stress/slip rate for each ox element
      do ixout=1,pb%mesh%nnglob,nxout_dyn
        ! Mark location of rupture front (max tau at given location)
        if (tau(ixout) > pb%tau_max_glob(ixout)) then
          pb%tau_max_glob(ixout) = tau(ixout)
          pb%t_rup_glob(ixout) = pb%time
        endif
        ! Mark location of max slip rate (at given location)
        if (v(ixout) > pb%v_max_glob(ixout)) then
          pb%v_max_glob(ixout) = v(ixout)
          pb%t_vmax_glob(ixout) = pb%time
        endif
      enddo
    endif

  endif

  ! End ox_dyn output
  !---------------------------------------------------------------------------
  ! Write QSB output

  if (write_QSB) then

    if (rising_edge) then

      unit = pb%ox%QSB_unit_pre
      open(unit, file='DYN_PRE.txt', status='REPLACE')
      call write_ox_lines(unit, pb%ox%fmt, pb%ox%objects_dyn, pb)
      pb%ox%pot_pre = pb%ot%pot
      close(unit)

    else if (falling_edge) then

      unit = pb%ox%QSB_unit_post
      open(unit, file='DYN_POST.txt', status='REPLACE')
      call write_ox_lines(unit, pb%ox%fmt, pb%ox%objects_dyn, pb)
      close(unit)

      ! Check for seismic moment
      ! TODO: Why? What is this? Purpose? Relevance?
      if ((pb%ot%pot - pb%ox%pot_pre) * pb%smu >= pb%DYN_M) then
        pb%ox%dyn_count = pb%ox%dyn_count + 1
        if (pb%ox%dyn_count > pb%DYN_SKIP) then
          pb%itstop = pb%it
          write(222,'(3i10,e24.14)') pb%it,pb%ot%ivmax,pb%ox%count,pb%time
        endif
      endif
    endif

  endif

  ! End QSB output
  !---------------------------------------------------------------------------
  ! Check dynamic state

  if (falling_edge .and. MPI_master) then
    ! Set global dynamic state back to 0
    pb%ox%dyn_stat = 0
    ! Increment event counter
    pb%ox%dyn_count = pb%ox%dyn_count + 1
  endif

  ! End check dynamic state
  !---------------------------------------------------------------------------

  ! All done

end subroutine ox_write

!=====================================================================
! Write ox data to file
subroutine write_ox_lines(unit, fmt, objects, pb)

  use problem_class
  type (problem_type), intent(inout) :: pb
  integer :: unit, ixout, iox, nn
  character(len=16), dimension(pb%ox%nox) :: fmt
  type(oxptr), dimension(pb%ox%nox) :: objects

  ! Number of nodes to loop over (either local or global)
  nn = size(objects(1)%p)

  ! Write header
  write(unit,'(a)') trim(pb%ox%header)
  ! Loop over all ox elements
  do ixout=1,nn,pb%ox%nxout
    ! Loop over all ox output quantities (except last one)
    do iox=1,pb%ox%nox-1
      ! Write ox output quantity, do not advance to next line
      write(unit, fmt(iox), advance='no') objects(iox)%p(ixout)
    enddo
    ! Write last ox output quantity, advance to next line
    write(unit, fmt(pb%ox%nox)) objects(pb%ox%nox)%p(ixout)
  enddo

end subroutine write_ox_lines

!=====================================================================
! distance between largest peak on the left half
! and largest peak on the right half
function crack_size(s,n)

  integer :: n
  double precision ::  s(n)
  double precision ::  stemp,smax,s1,s2,xL,xR
  integer :: i,iL,iR,imin
  double precision :: crack_size

  if (n < 2) then
    crack_size = 0d0
    return
  endif

  ! assuming minimum stressing rate at nucleation point
  !      smin=s(1)
  !      do i=1,n
  !        if (s(i).le.smin) then
  !          imin=i
  !          smin=s(i)
  !        endif
  !      enddo

  ! WARNING: assumes one peak on x<0 and one peak on x>0
  imin=n/2

  iL=1
  stemp=0d0;
  do i=1,imin
    if (s(i) > stemp) then
      iL = i
      stemp = s(i)
    end if
  end do

  iR=imin
  stemp=0d0
  do i=imin,n
    if (s(i) >= stemp) then
      iR = i
      stemp = s(i)
    end if
  end do


  !iL = maxloc(s(1:imin))
  smax=s(iL)
  xL = dble(iL);
  if (iL > 1) then
    s1 = 0.5d0*( s(iL+1)-s(iL-1) )
    s2 = s(iL-1)-2d0*s(iL)+s(iL+1)
    if (s2 /= 0d0) xL = xL-s1/s2
  endif

  !iR = maxloc(s(imin:n))
  smax=s(iR)
  xR = dble(iR);
  if (iR < n) then
    s1 = 0.5d0*( s(iR+1)-s(iR-1) )
    s2 = s(iR-1)-2d0*s(iR)+s(iR+1)
    if (s2 /= 0d0) xR = xR-s1/s2
  endif

  crack_size = xR-xL

end function crack_size

!=====================================================================
! Initiate MPI global gather

subroutine init_pb_global(pb)

  use problem_class

  type(problem_type), intent(inout) :: pb
  integer :: n

  n = pb%mesh%nnglob

  ! Check if global variables are already allocated (should not be!)
  if (.not. pb%allocated_glob) then

    ! Base quantities
    allocate( pb%v_glob(n), pb%theta_glob(n), pb%tau_glob(n), &
              pb%dtau_dt_glob(n), pb%slip_glob(n), pb%sigma_glob(n))
    ! If thermal pressurisation is requested, allocate P and T
    if (pb%features%tp == 1) then
      allocate(pb%P_glob(n), pb%T_glob(n))
    endif

    ! Allocate rupture max stats
    allocate( pb%tau_max_glob(n), pb%t_rup_glob(n), &
              pb%v_max_glob(n), pb%t_vmax_glob(n))

    pb%tau_max_glob = 0d0
    pb%t_rup_glob = 0d0
    pb%v_max_glob = 0d0
    pb%t_vmax_glob = 0d0

    ! Update flag
    pb%allocated_glob = .true.

  ! If global quantities are already allocated, something is wrong!
  else
    write(6, *) "Error in output.f90::init_pb_global"
    write(6, *) "Global quantities are allocated while they shouldn't be"
    stop "Terminating..."
  endif

end subroutine init_pb_global

!=====================================================================
! Collect global fault nodes to master processor for outputs
subroutine pb_global(pb)

  use problem_class
  use mesh, only: nnLocal_perproc, nnoffset_glob_perproc
  use my_mpi, only: gather_allvdouble_root

  type(problem_type), intent(inout) :: pb
  integer :: i, nnLocal, nnGlobal

  nnLocal = pb%mesh%nn
  nnGlobal = pb%mesh%nnglob

  ! Loop over proc-specific quantities
  ! Skip: x, y, z, time
  do i=5,pb%ox%nox
    ! Initialise to zero
    write(6, *) i, size(pb%ox%objects_glob(i)%p)
    pb%ox%objects_glob(i)%p = 0d0
    ! Call MPI gather
    write(6, *) i, size(pb%ox%objects_glob(i)%p)
    call gather_allvdouble_root(  pb%ox%objects_loc(i)%p, nnLocal, &
                                  pb%ox%objects_glob(i)%p, nnLocal_perproc, &
                                  nnoffset_glob_perproc, nnGlobal)
    write(6, *) i, size(pb%ox%objects_glob(i)%p)
  enddo

end subroutine pb_global
!=====================================================================


end module output
