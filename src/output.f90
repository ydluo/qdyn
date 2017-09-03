module output

! OUTPUT: This module manages outputs

  implicit none
  private
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
  type (problem_type), intent(inout) :: pb

  write(121,*) pb%time, pb%tmax

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
  integer :: i,n,nsta,ista,ik
  double precision :: dmin2, d2

  n = mesh_get_size(pb%mesh)
  pb%ot%lcnew = dble(n)
  pb%ot%llocnew = dble(n)

  if (is_MPI_parallel()) then

   ! find stations
    dmin2 = 0.01d0*minval(pb%mesh%dx) ! distance tolerance = 1% grid size
    dmin2 = dmin2*dmin2
    nsta = 1 !NOTE: currently only one station implemented
    do ista=1,nsta
      do ik=1,n
        d2 = (pb%mesh%x(ik)-pb%ot%xsta(ista))**2 &
           + (pb%mesh%y(ik)-pb%ot%ysta(ista))**2 &
           + (pb%mesh%z(ik)-pb%ot%zsta(ista))**2
        if (d2 < dmin2) then
          pb%ot%ic=ik
          write(6,*) 'Processor: ',my_mpi_tag(),', station ',ista, &
                     ' found, distance mismatch = ',d2
          exit
        endif
      enddo
    enddo

    if (OUT_MASTER .and. is_mpi_master() ) then
      pb%ot%unit = 18
      write(pb%ot%unit,'(a)')'# macroscopic values:'
      write(pb%ot%unit,'(a)')'# 1=t'
      write(pb%ot%unit,'(a)')'# values at selected point:'
      write(pb%ot%unit,'(a)')'# 2=V, 3=theta, 4=V*theta/dc, 5=tau, 6=slip'
      if (pb%features%tp == 1) write(pb%ot%unit,'(a)')'# 7=P, 8=T'
      close(pb%ot%unit)
    endif
!JPA WARNING Implementation in progress
!    pb%ot%unit = 22
!    write(pb%ot%unit,'(a)')'# Seismicity record:'
!    write(pb%ot%unit,'(a)')'# 1=loc, 2=t, 3=v'
!  pb%ot%unit = 10000
!  do i=1,n
!    if (pb%ot%iot(i) == 1) then
!      pb%ot%unit = pb%ot%unit+1
!      write(pb%ot%unit,'(a,i10)')'# nx= ', i
!      write(pb%ot%unit,'(a)')'# 1=t, 2=V, 3=theta, 4=tau, 5=slip, 6=sigma'
!    endif
!  enddo

else

  pb%ot%unit = 18
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

  pb%ot%unit = 22
  write(pb%ot%unit,'(a)')'# Seismicity record:'
  write(pb%ot%unit,'(a)')'# 1=loc, 2=t, 3=v'

  pb%ot%unit = 10000
  do i=1,n
    if (pb%ot%iot(i) == 1) then
      pb%ot%unit = pb%ot%unit+1
      write(pb%ot%unit,'(a,i10)')'# nx= ', i
      write(pb%ot%unit,'(a)')'# 1=t, 2=V, 3=theta, 4=tau, 5=slip, 6=sigma'
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
  use constants, only : OUT_MASTER, BIN_OUTPUT
  use my_mpi, only : is_MPI_parallel, is_mpi_master
  use mesh, only : mesh_get_size

  type (problem_type), intent(inout) :: pb

  integer :: i,n

  n = mesh_get_size(pb%mesh)
  allocate ( pb%ox%t_rup(n), pb%ox%tau_max(n), pb%ox%v_max(n), pb%ox%t_vmax(n) )

  if (pb%ox%i_ox_seq == 0) then
    pb%ox%unit = 19
  else
!    pb%ox%unit = 1000
    pb%ox%unit = 999
  endif

  pb%ox%count=0
  pb%ox%dyn_count=0
  pb%ox%dyn_count2=0
  do i=1,pb%mesh%nn,pb%ox%nxout
    pb%ox%count = pb%ox%count+1
  enddo

  if (is_MPI_parallel()) then
    if (OUT_MASTER .and. is_mpi_master() ) then
      pb%ox%countglob=0
      do i=1,pb%mesh%nnglob, pb%ox%nxout
        pb%ox%countglob = pb%ox%countglob+1
      enddo
     if (pb%ox%unit==19) write(pb%ox%unit,'(a,i10)')'# nx= ',pb%ox%countglob
    endif

  else
     if (.not.BIN_OUTPUT) then
         write(pb%ox%unit,'(a,i10)') '# nx= ',pb%ox%count
     else
         open(pb%ox%unit,form='unformatted',access='stream')
         write(pb%ox%unit) pb%ox%count
         do i=1,pb%mesh%nn,pb%ox%nxout
           write(pb%ox%unit) pb%mesh%x(i)
         enddo
     endif
  endif

  pb%ox%dyn_stat = 0
  pb%ox%dyn_stat2 = 0

  pb%ox%pot_pre = 0.d0

end subroutine ox_init

!=====================================================================
! Export timeseries
subroutine ot_write(pb)

  use problem_class
  use constants, only : OCTAVE_OUTPUT, BIN_OUTPUT
  use my_mpi, only : is_MPI_parallel

  type (problem_type), intent(inout) :: pb
  integer :: i,ios
  character(30) :: ot_fmt

  pb%ot%unit = 18
  if (is_MPI_parallel()) then
   ! if "ic" station is in this processor
    if (pb%ot%ic>0) then
      open(pb%ot%unit,access='APPEND',status='old',iostat=ios)
      if (ios>0) stop 'Fatal error: ot_write: Error opening a fort.18 file'
     !JPA add test for the first time we try to open this file but it does not exist yet
     !JPA add test to prevent appenaingd data to a file from a previous simulation
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
      if (pb%it+1 == pb%itstop) close(pb%ot%unit)
    endif

  endif

  if (is_MPI_parallel()) return
 !JPA warning: the ot outputs below are not yet implemented in parallel

 ! output slip velocity maxima at selected nodes
  pb%ot%unit = 22
  do i=1,pb%mesh%nn
    if ((pb%ot%iasp(i) == 1) .and. (pb%ot%v_pre(i) >= pb%ot%v_th) .and.      &
        (pb%v(i) < pb%ot%v_pre(i)) .and. (pb%ot%v_pre(i) >= pb%ot%v_pre2(i))) then
      write(pb%ot%unit,'(i10,2e24.16)') i, pb%time, pb%ot%v_pre(i)
    endif
  enddo
  pb%ot%v_pre2=pb%ot%v_pre
  pb%ot%v_pre=pb%v

  pb%ot%unit = 10000

 ! output time series at selected nodes
  do i=1,pb%mesh%nn
    if (pb%ot%iot(i) == 1) then
      pb%ot%unit = pb%ot%unit+1
      write(pb%ot%unit,'(e24.16,5e14.6)') pb%time, pb%v(i), &
        pb%theta(i), pb%tau(i), pb%slip(i), pb%sigma(i)
    endif
  enddo

end subroutine ot_write

!=====================================================================
! Export snapshots
subroutine ox_write(pb)

  use problem_class
  use constants, only: OUT_MASTER, BIN_OUTPUT
  use my_mpi, only: is_MPI_parallel, is_mpi_master, my_mpi_tag, synchronize_all

  type (problem_type), intent(inout) :: pb

  integer :: ixout
  character(len=256) :: fileproc

if (is_MPI_parallel()) then
! In progress
! if (mod(pb%it-1,pb%ot%ntout) == 0 .or. pb%it == pb%itstop) then
 if (mod(pb%it,pb%ot%ntout) == 0 .or. pb%it == pb%itstop) then
  if (OUT_MASTER) then
  ! Collecting global nodes
    call pb_global(pb)
    if (is_mpi_master()) then
      pb%ot%ivmaxglob = maxloc(pb%v_glob,1)
      ! Writing fault points in single file fort.19
      if (pb%ox%i_ox_seq == 0) then
        write(pb%ox%unit,'(a,2i8,e14.6)') '# x y z t v theta dtau tau_dot slip sigma ',&
                                          pb%it,pb%ot%ivmaxglob,pb%time
        do ixout=1,pb%mesh%nnglob,pb%ox%nxout
          write(pb%ox%unit,'(3e15.7,e24.16,6e15.7)') pb%mesh%xglob(ixout),pb%mesh%yglob(ixout),&
          pb%mesh%zglob(ixout),pb%time,pb%v_glob(ixout),pb%theta_glob(ixout),&
          pb%tau_glob(ixout),pb%dtau_dt_glob(ixout),pb%slip_glob(ixout), pb%sigma_glob(ixout)
        enddo
      else
      !Writing in File output in snapshots. fort.XXXXXX
      pb%ox%unit = pb%ox%unit + 1
      write(pb%ox%unit,'(3i10,e24.14)') pb%it,pb%ot%ivmaxglob,pb%ox%countglob,pb%time
      write(pb%ox%unit,'(a)') '#  x  y  z  t  v  theta  dtau  tau_dot  slip '
      do ixout=1,pb%mesh%nnglob,pb%ox%nxout
        write(pb%ox%unit,'(3e15.7,e24.14,6e15.7)')       &
          pb%mesh%xglob(ixout),pb%mesh%yglob(ixout),pb%mesh%zglob(ixout),pb%time,     &
          pb%v_glob(ixout),pb%theta_glob(ixout),pb%tau_glob(ixout),   &
          pb%dtau_dt_glob(ixout),pb%slip_glob(ixout), pb%sigma_glob(ixout)
      enddo
      close(pb%ox%unit) !Closing snapshot
      endif
    endif
  else
  !local
      !Each processor writes an output file.
      pb%ox%unit = pb%ox%unit + 1
      write(fileproc,'(a,i6.6,a,a)') 'fort.',pb%ox%unit,'_proc',my_mpi_tag()
      open(pb%ox%unit,file=fileproc(1:len_trim(fileproc)),status='replace',form='formatted',action='write')
      write(pb%ox%unit,'(3i10,e24.14)') pb%it,pb%ot%ivmax,pb%ox%count,pb%time
      write(pb%ox%unit,'(a)') '#  x  y  z  t  v  theta  dtau  tau_dot  slip '
      do ixout=1,pb%mesh%nn,pb%ox%nxout
        write(pb%ox%unit,'(3e15.7,e24.14,6e15.7)')       &
          pb%mesh%x(ixout),pb%mesh%y(ixout),pb%mesh%z(ixout),pb%time,     &
          pb%v(ixout),pb%theta(ixout),pb%tau(ixout),   &
          pb%dtau_dt(ixout),pb%slip(ixout), pb%sigma(ixout)
      enddo
      close(pb%ox%unit) !Closing snapshot
  endif
!  close(pb%ox%unit)
 endif

 if (pb%ox%i_ox_dyn == 1) then
   if (OUT_MASTER) then
     call pb_global(pb)
     if (is_mpi_master()) then
      if (pb%ox%dyn_stat2 == 0 .and. pb%vmaxglob >= pb%DYN_th_on ) then
        pb%ox%dyn_stat2 = 1
        do ixout=1,pb%mesh%nnglob,pb%ox%nxout_dyn
          pb%tau_max_glob(ixout) = pb%tau_glob(ixout)
          pb%t_rup_glob(ixout) = pb%time
          pb%v_max_glob(ixout) = pb%v_glob(ixout)
          pb%t_vmax_glob(ixout) = pb%time
        enddo
        write(20001+3*pb%ox%dyn_count2,'(3i10,e24.14)')   &
              pb%it,pb%ot%ivmaxglob,pb%ox%countglob,pb%time
        write(20001+3*pb%ox%dyn_count2,'(a)') '#  x  y  z  t  v  theta  dtau  tau_dot  slip sigma'
        do ixout=1,pb%mesh%nnglob,pb%ox%nxout_dyn
           write(20001+3*pb%ox%dyn_count2,'(3e15.7,e24.14,6e15.7)')       &
           pb%mesh%xglob(ixout),pb%mesh%yglob(ixout),pb%mesh%zglob(ixout),pb%time,     &
           pb%v_glob(ixout),pb%theta_glob(ixout),pb%tau_glob(ixout),   &
           pb%dtau_dt_glob(ixout),pb%slip_glob(ixout), pb%sigma_glob(ixout)
        enddo
        close(20001+3*pb%ox%dyn_count2)
      endif

      if (pb%ox%dyn_stat2 == 1) then
        do ixout=1,pb%mesh%nnglob,pb%ox%nxout_dyn
          if (pb%tau_glob(ixout) > pb%tau_max_glob(ixout)) then
             pb%tau_max_glob(ixout) = pb%tau_glob(ixout)
             pb%t_rup_glob(ixout) = pb%time
          endif
          if (pb%v_glob(ixout) > pb%v_max_glob(ixout)) then
             pb%v_max_glob(ixout) = pb%v_glob(ixout)
             pb%t_vmax_glob(ixout) = pb%time
          endif
        enddo
      endif

      if (pb%ox%dyn_stat2 == 1 .and. pb%vmaxglob <= pb%DYN_th_off ) then
        pb%ox%dyn_stat2 = 0
        write(20002+3*pb%ox%dyn_count2,'(3i10,e24.14)')   &
              pb%it,pb%ot%ivmaxglob,pb%ox%countglob,pb%time
        write(20002+3*pb%ox%dyn_count2,'(a)') '#  x  y  z  t  v  theta  dtau  tau_dot  slip sigma'
        do ixout=1,pb%mesh%nnglob,pb%ox%nxout_dyn
           write(20002+3*pb%ox%dyn_count2,'(3e15.7,e24.14,6e15.7)')       &
           pb%mesh%xglob(ixout),pb%mesh%yglob(ixout),pb%mesh%zglob(ixout),pb%time,     &
           pb%v_glob(ixout),pb%theta_glob(ixout),pb%tau_glob(ixout),   &
           pb%dtau_dt_glob(ixout),pb%slip_glob(ixout),pb%sigma_glob(ixout)
        enddo
        close(20002+3*pb%ox%dyn_count2)

        write(20003+3*pb%ox%dyn_count2,'(a)') '#  x  y  z  t_rup tau_max t_vmax vmax'
        do ixout=1,pb%mesh%nnglob,pb%ox%nxout_dyn
           write(20003+3*pb%ox%dyn_count2,'(3e15.7,4e28.20)')       &
           pb%mesh%xglob(ixout),pb%mesh%yglob(ixout),pb%mesh%zglob(ixout),       &
           pb%t_rup_glob(ixout),pb%tau_max_glob(ixout),pb%t_vmax_glob(ixout),pb%v_max_glob(ixout)
        enddo
        close(20003+3*pb%ox%dyn_count2)
        pb%ox%dyn_count2 = pb%ox%dyn_count2 + 1
      endif
     endif
   else
   ! writing in chuncks
   endif
 endif

else
 if (mod(pb%it-1,pb%ot%ntout) == 0 .or. pb%it == pb%itstop) then
  if (pb%ox%i_ox_seq == 0) then

    ! <begin TP output conditional>
    ! SEISMIC: if thermal pressurisation is requested, write P and T output
    if (.not.BIN_OUTPUT) then
      if (pb%features%tp == 1) then
        write(pb%ox%unit,'(a,2i8,e14.6)')'# x t v theta dtau tau_dot slip sigma_e P T',pb%it,pb%ot%ivmax,pb%time
        ! JPA: this output should also contain y and z
        do ixout=1,pb%mesh%nn,pb%ox%nxout
          write(pb%ox%unit,'(e15.7,e24.16,8e15.7)') pb%mesh%x(ixout),pb%time,pb%v(ixout),   &
            pb%theta(ixout),pb%tau(ixout), pb%dtau_dt(ixout),pb%slip(ixout), &
            pb%sigma(ixout)-pb%tp%P(ixout), pb%tp%P(ixout), pb%tp%T(ixout)
        enddo
      else
        write(pb%ox%unit,'(a,2i8,e14.6)')'# x t v theta dtau tau_dot slip sigma',pb%it,pb%ot%ivmax,pb%time
        ! JPA: this output should also contain y and z
        do ixout=1,pb%mesh%nn,pb%ox%nxout
          write(pb%ox%unit,'(e15.7,e24.16,6e15.7)') pb%mesh%x(ixout),pb%time,pb%v(ixout),   &
            pb%theta(ixout),pb%tau(ixout), pb%dtau_dt(ixout),pb%slip(ixout), pb%sigma(ixout)
        enddo
      endif
    else
      ! SEISMIC: need to add P/T here too...
      write(pb%ox%unit) pb%time
      do ixout=1,pb%mesh%nn,pb%ox%nxout
         write(pb%ox%unit) pb%v(ixout),pb%theta(ixout),pb%tau(ixout), &
           pb%dtau_dt(ixout), pb%slip(ixout), pb%sigma(ixout)
      enddo
      if (pb%it+1 == pb%itstop) close(pb%ox%unit)

    endif
    ! <end TP output conditional>

  else
    pb%ox%unit = pb%ox%unit + 1
    write(pb%ox%unit,'(3i10,e24.14)') pb%it,pb%ot%ivmax,pb%ox%count,pb%time
    write(pb%ox%unit,'(a)') '#  x  y  z  t  v  theta  dtau  tau_dot  slip  sigma'
    do ixout=1,pb%mesh%nn,pb%ox%nxout
      write(pb%ox%unit,'(3e15.7,e24.14,6e15.7)')       &
        pb%mesh%x(ixout),pb%mesh%y(ixout),pb%mesh%z(ixout),pb%time,     &
        pb%v(ixout),pb%theta(ixout),pb%tau(ixout),   &
        pb%dtau_dt(ixout),pb%slip(ixout), pb%sigma(ixout)
    enddo
    close(pb%ox%unit)
  endif
 endif

  if (pb%ox%i_ox_dyn == 1) then

    if (pb%ox%dyn_stat2 == 0 .and. pb%v(pb%ot%ivmax) >= pb%DYN_th_on ) then
      pb%ox%dyn_stat2 = 1
      do ixout=1,pb%mesh%nn,pb%ox%nxout_dyn
        pb%ox%tau_max(ixout) = pb%tau(ixout)
        pb%ox%t_rup(ixout) = pb%time
        pb%ox%v_max(ixout) = pb%v(ixout)
        pb%ox%t_vmax(ixout) = pb%time
      enddo
      write(20001+3*pb%ox%dyn_count2,'(3i10,e24.14)')   &
            pb%it,pb%ot%ivmax,pb%ox%count,pb%time
      write(20001+3*pb%ox%dyn_count2,'(a)') '#  x  y  z  t  v  theta  dtau  tau_dot  slip sigma'
      do ixout=1,pb%mesh%nn,pb%ox%nxout_dyn
        write(20001+3*pb%ox%dyn_count2,'(3e15.7,e24.14,6e15.7)')       &
          pb%mesh%x(ixout),pb%mesh%y(ixout),pb%mesh%z(ixout),pb%time,     &
          pb%v(ixout),pb%theta(ixout),pb%tau(ixout),   &
          pb%dtau_dt(ixout),pb%slip(ixout), pb%sigma(ixout)
      enddo
      close(20001+3*pb%ox%dyn_count2)
    endif

    if (pb%ox%dyn_stat2 == 1) then
      do ixout=1,pb%mesh%nn,pb%ox%nxout_dyn
        if (pb%tau(ixout) > pb%ox%tau_max(ixout)) then
          pb%ox%tau_max(ixout) = pb%tau(ixout)
          pb%ox%t_rup(ixout) = pb%time
        endif
        if (pb%v(ixout) > pb%ox%v_max(ixout)) then
          pb%ox%v_max(ixout) = pb%v(ixout)
          pb%ox%t_vmax(ixout) = pb%time
        endif
      enddo
    endif

    if (pb%ox%dyn_stat2 == 1 .and. pb%v(pb%ot%ivmax) <= pb%DYN_th_off ) then
      pb%ox%dyn_stat2 = 0
      write(20002+3*pb%ox%dyn_count2,'(3i10,e24.14)')   &
            pb%it,pb%ot%ivmax,pb%ox%count,pb%time
      write(20002+3*pb%ox%dyn_count2,'(a)') '#  x  y  z  t  v  theta  dtau  tau_dot  slip sigma'
      do ixout=1,pb%mesh%nn,pb%ox%nxout_dyn
        write(20002+3*pb%ox%dyn_count2,'(3e15.7,e24.14,6e15.7)')       &
          pb%mesh%x(ixout),pb%mesh%y(ixout),pb%mesh%z(ixout),pb%time,     &
          pb%v(ixout),pb%theta(ixout),pb%tau(ixout),   &
          pb%dtau_dt(ixout),pb%slip(ixout),pb%sigma(ixout)
      enddo
      close(20002+3*pb%ox%dyn_count2)

      write(20003+3*pb%ox%dyn_count2,'(a)') '#  x  y  z  t_rup tau_max t_vmax vmax'
      do ixout=1,pb%mesh%nn,pb%ox%nxout_dyn
        write(20003+3*pb%ox%dyn_count2,'(3e15.7,4e28.20)')       &
          pb%mesh%x(ixout),pb%mesh%y(ixout),pb%mesh%z(ixout),       &
          pb%ox%t_rup(ixout),pb%ox%tau_max(ixout),pb%ox%t_vmax(ixout),pb%ox%v_max(ixout)
      enddo
      close(20003+3*pb%ox%dyn_count2)

      pb%ox%dyn_count2 = pb%ox%dyn_count2 + 1
    endif

  endif

  if (pb%DYN_FLAG == 1) then

    if (pb%ox%dyn_stat == 0 .and. pb%v(pb%ot%ivmax) >= pb%DYN_th_on ) then
      pb%ox%dyn_stat = 1
      OPEN (UNIT = 100, FILE='DYN_PRE.txt', STATUS='REPLACE')
      write(100,'(3i10,e24.14)') pb%it,pb%ot%ivmax,pb%ox%count,pb%time
      write(100,'(a)') '#  x  y  z  t  v  theta  dtau  tau_dot  slip sigma '
      do ixout=1,pb%mesh%nn,pb%ox%nxout_dyn
        write(pb%ox%unit,'(3e15.7,e24.14,6e15.7)')       &
          pb%mesh%x(ixout),pb%mesh%y(ixout),pb%mesh%z(ixout),pb%time,     &
          pb%v(ixout),pb%theta(ixout),pb%tau(ixout),   &
          pb%dtau_dt(ixout),pb%slip(ixout),pb%sigma(ixout)
      enddo
      pb%ox%pot_pre = pb%ot%pot
      CLOSE(100)
    endif

    if (pb%ox%dyn_stat == 1 .and. pb%v(pb%ot%ivmax) <= pb%DYN_th_off ) then
      pb%ox%dyn_stat = 0
      OPEN (UNIT = 101, FILE='DYN_POST.txt', STATUS='REPLACE')
      write(101,'(3i10,e24.14)') pb%it,pb%ot%ivmax,pb%ox%count,pb%time
      write(101,'(a)') '#  x  y  z  t  v  theta  dtau  tau_dot  slip sigma'
      do ixout=1,pb%mesh%nn,pb%ox%nxout_dyn
        write(pb%ox%unit,'(3e15.7,e24.14,6e15.7)')       &
          pb%mesh%x(ixout),pb%mesh%y(ixout),pb%mesh%z(ixout),pb%time,     &
          pb%v(ixout),pb%theta(ixout),pb%tau(ixout),   &
          pb%dtau_dt(ixout),pb%slip(ixout),pb%sigma(ixout)
      enddo
      CLOSE(101)
      if ((pb%ot%pot-pb%ox%pot_pre)*pb%smu >= pb%DYN_M) then
        pb%ox%dyn_count = pb%ox%dyn_count + 1
        if (pb%ox%dyn_count > pb%DYN_SKIP) then
          pb%itstop = pb%it
          write(222,'(3i10,e24.14)') pb%it,pb%ot%ivmax,pb%ox%count,pb%time
        endif
      endif
    endif

  endif

endif

end subroutine ox_write

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
! Collect global fault nodes to master processor for outputs
subroutine pb_global(pb)

  use mesh, only: nnLocal_perproc,nnoffset_glob_perproc
  use problem_class
  use my_mpi, only: my_mpi_rank, gather_allvdouble_root

  type(problem_type), intent(inout) :: pb
  integer :: nnLocal,nnGlobal

  nnLocal= pb%mesh%nn
  nnGlobal= pb%mesh%nnglob

  if (.not.allocated(pb%v_glob)) then

    allocate(pb%v_glob(nnGlobal),pb%theta_glob(nnGlobal),pb%tau_glob(nnGlobal),&
           pb%slip_glob(nnGlobal),pb%sigma_glob(nnGlobal),pb%dtau_dt_glob(nnGlobal))

    allocate(pb%tau_max_glob(nnGlobal),pb%t_rup_glob(nnGlobal),&
             pb%v_max_glob(nnGlobal),pb%t_vmax_glob(nnGlobal))

  endif

  pb%v_glob=0
  pb%theta_glob=0
  pb%tau_glob=0
  pb%dtau_dt_glob=0
  pb%slip_glob=0
  pb%sigma_glob=0

  call gather_allvdouble_root(pb%v,nnLocal,pb%v_glob,nnLocal_perproc, &
                           nnoffset_glob_perproc,nnGlobal)
  call gather_allvdouble_root(pb%theta,nnLocal,pb%theta_glob,nnLocal_perproc, &
                           nnoffset_glob_perproc,nnGlobal)
  call gather_allvdouble_root(pb%tau,nnLocal,pb%tau_glob,nnLocal_perproc, &
                           nnoffset_glob_perproc,nnGlobal)
  call gather_allvdouble_root(pb%slip,nnLocal,pb%slip_glob,nnLocal_perproc, &
                           nnoffset_glob_perproc,nnGlobal)
  call gather_allvdouble_root(pb%sigma,nnLocal,pb%sigma_glob,nnLocal_perproc, &
                           nnoffset_glob_perproc,nnGlobal)

  pb%tau_max_glob=0d0
  pb%t_rup_glob=0d0
  pb%v_max_glob=0d0
  pb%t_vmax_glob=0d0

end subroutine pb_global
!=====================================================================


end module output
