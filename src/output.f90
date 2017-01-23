module output

! OUTPUT: This module manages outputs
!

!  use some_module_1

  implicit none
  private
  public :: screen_init, ot_init, ox_init, &
            screen_write, ot_write, ox_write,  &
            time_write, crack_size

contains

!=====================================================================
!output initilized field to screen
subroutine screen_init(pb)

  use problem_class
  use constants, only : PI
  type (problem_type), intent(inout) :: pb

  write(6,*) '**Field Initialized.**'
  write(6,*) 'Values at selected point of the fault:'
    if (pb%mesh%dim == 0 .or.pb%mesh%dim == 1) then
      ! SEISMIC: in spring-block case, pb%kernel%k2f is not allocated
      if (pb%kernel%k2f%finite == 1 .or. pb%mesh%nn == 1) then
        write(6,*) 'K/Kc = ',(PI*pb%smu/pb%mesh%Lfault)/   &
          (pb%sigma(pb%ot%ic)*(pb%b(pb%ot%ic)-pb%a(pb%ot%ic))/pb%dc(pb%ot%ic))
        write(6,*) 'K/Kb = ',(PI*pb%smu/pb%mesh%Lfault)/   &
          (pb%sigma(pb%ot%ic)*pb%b(pb%ot%ic)/pb%dc(pb%ot%ic))
      else
        write(6,*) 'K/Kc = ',(PI*pb%smu/pb%mesh%W)/   &
          (pb%sigma(pb%ot%ic)*(pb%b(pb%ot%ic)-pb%a(pb%ot%ic))/pb%dc(pb%ot%ic))
        write(6,*) 'K/Kb = ',(PI*pb%smu/pb%mesh%W)/   &
          (pb%sigma(pb%ot%ic)*pb%b(pb%ot%ic)/pb%dc(pb%ot%ic))
      endif
    end if

! YD:  should we out put 2D mesh K?

    write(6,*)
    write(6,*) '    it,  dt (secs), time (yrs), vmax (m/s), sigma(MPa)'

end subroutine screen_init



!=====================================================================
!output one step to screen
subroutine screen_write(pb)
  
  use constants, only : YEAR
  use problem_class
  use my_mpi, only : is_MPI_parallel, is_mpi_master, max_allproc

  type (problem_type), intent(inout) :: pb
  double precision :: vtempglob
  integer :: i


  if (is_MPI_parallel()) then
    !Finding the global max
    vtempglob=0d0
    do i=1,pb%mesh%nn
     if ( pb%v(i) > vtempglob) then
       vtempglob = pb%v(i)
       pb%ot%ivmax = i
     end if
    end do
    call max_allproc(pb%v(pb%ot%ivmax),pb%vmaxglob)
    call max_allproc(pb%sigma(pb%ot%ivmax),pb%sigma_vmaxglob)
    if (is_mpi_master()) write(6,'(i7,x,4(e11.3,x),i5)') pb%it, pb%dt_did, pb%time/YEAR,&
                              pb%vmaxglob, pb%sigma_vmaxglob/1.0D6
  else
    write(6,'(i7,x,4(e11.3,x),i5)') pb%it, pb%dt_did, pb%time/YEAR,    &
                            pb%v(pb%ot%ivmax), pb%sigma(pb%ot%ivmax)/1.0D6
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
! write ot file header
subroutine ot_init(pb)

  use problem_class
  use constants, only: OUT_MASTER
  use my_mpi, only : is_MPI_parallel, is_mpi_master

  type (problem_type), intent(inout) :: pb
  integer :: i

  pb%ot%lcnew = dble(pb%mesh%nn)
  pb%ot%llocnew = dble(pb%mesh%nn)

  if (is_MPI_parallel()) then

    if (OUT_MASTER .and. is_mpi_master() ) then
      pb%ot%unit = 18
      write(pb%ot%unit,'(a)')'# macroscopic values:'
      write(pb%ot%unit,'(a)')'# 1=t'
      write(pb%ot%unit,'(a)')'# values at selected point:'
      write(pb%ot%unit,'(a)')'# 2=V, 3=theta, 4=V*theta/dc, 5=tau, 6=slip'
      close(pb%ot%unit)
    endif
! In progress
!    pb%ot%unit = 22
!    write(pb%ot%unit,'(a)')'# Seismicity record:'
!    write(pb%ot%unit,'(a)')'# 1=loc, 2=t, 3=v'
!  pb%ot%unit = 10000
!  do i=1,pb%mesh%nn
!    if (pb%iot(i) == 1) then
!      pb%ot%unit = pb%ot%unit+1
!      write(pb%ot%unit,'(a,i10)')'# nx= ', i
!      write(pb%ot%unit,'(a)')'# 1=t, 2=V, 3=theta, 4=tau, 5=slip, 6=sigma'
!    endif
!  enddo

else

  pb%ot%unit = 18
  write(pb%ot%unit,'(a)')'# macroscopic values:'
  write(pb%ot%unit,'(a)')'# 1=t,2=loc_size,3=crack_size,4=potcy,5=pot_rate'
  write(pb%ot%unit,'(a)')'# values at selected point:'
  write(pb%ot%unit,'(a)')'# 6=V, 7=theta, 8=V*theta/dc, 9=tau, 10=slip'
  write(pb%ot%unit,'(a)')'# values at max(V) location:'
  write(pb%ot%unit,'(a)')'# 11=x, 12=V, 13=theta, 14=omeg, 15=tau, 16=slip, 17=sigma'

  pb%ot%unit = 22
  write(pb%ot%unit,'(a)')'# Seismicity record:'
  write(pb%ot%unit,'(a)')'# 1=loc, 2=t, 3=v'

  pb%ot%unit = 10000
  do i=1,pb%mesh%nn
    if (pb%iot(i) == 1) then
      pb%ot%unit = pb%ot%unit+1
      write(pb%ot%unit,'(a,i10)')'# nx= ', i
      write(pb%ot%unit,'(a)')'# 1=t, 2=V, 3=theta, 4=tau, 5=slip, 6=sigma'
    endif
  enddo

endif

end subroutine ot_init



!=====================================================================
! write ox file header
subroutine ox_init(pb)

  use problem_class
  use constants, only : OUT_MASTER
  use my_mpi, only : is_MPI_parallel, is_mpi_master

  type (problem_type), intent(inout) :: pb

  integer :: i

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
    write(pb%ox%unit,'(a,i10)')'# nx= ',pb%ox%count
  endif

end subroutine ox_init

!=====================================================================
! Export timeseries
subroutine ot_write(pb)

  use problem_class
  use constants, only : OCTAVE_OUTPUT
  use my_mpi, only : is_MPI_parallel

  type (problem_type), intent(inout) :: pb
  integer :: i,ios
  character(30) :: ot_fmt

  if (is_MPI_parallel()) then

 ! Working on parallel outputs for timeseries
    ! if one station is found in this processor
    if (pb%station_found) then
      pb%ot%unit = 18
      open(pb%ot%unit,access='APPEND',status='old',iostat=ios)
      if (ios==0) then
      !Writing in File output.
      ot_fmt = '(e24.16,5e14.6)'
      write(pb%ot%unit,ot_fmt) pb%time, pb%v(pb%ot%ic), pb%theta(pb%ot%ic), &
      pb%v(pb%ot%ic)*pb%theta(pb%ot%ic)/pb%dc(pb%ot%ic), &
      pb%tau(pb%ot%ic), pb%slip(pb%ot%ic)
      else
       stop 'Error opening a fort.18 file'
      endif
    endif

  else

  pb%ot%unit = 18
  if (OCTAVE_OUTPUT) then
   ! for Octave: comma as field delimiter and no spaces
    ot_fmt = '(g0.16,16(",",g0.6))'
  else
    ot_fmt = '(e24.16,16e14.6)'
  endif
  write(pb%ot%unit,ot_fmt) pb%time, pb%ot%llocnew*pb%mesh%dx,  &
    pb%ot%lcnew*pb%mesh%dx, pb%pot, pb%pot_rate,    &
    pb%v(pb%ot%ic), pb%theta(pb%ot%ic),  &
    pb%v(pb%ot%ic)*pb%theta(pb%ot%ic)/pb%dc(pb%ot%ic), &
    pb%tau(pb%ot%ic), pb%slip(pb%ot%ic),    &
!   for ivmax
    pb%mesh%x(pb%ot%ivmax), pb%v(pb%ot%ivmax), pb%theta(pb%ot%ivmax),   &
    pb%v(pb%ot%ivmax)*pb%theta(pb%ot%ivmax)/pb%dc(pb%ot%ivmax),    &
    pb%tau(pb%ot%ivmax), pb%slip(pb%ot%ivmax), pb%sigma(pb%ot%ivmax)

  pb%ot%unit = 22
  do i=1,pb%mesh%nn
    if ((pb%iasp(i) == 1) .and. (pb%v(i) >= pb%v_th) .and.      &
        (pb%v(i) < pb%v_pre(i)) .and. (pb%v_pre(i) >= pb%v_pre2(i))) then
      write(pb%ot%unit,'(i10,2e24.16)') i, pb%time, pb%v(i)
    endif
  enddo
  pb%v_pre2=pb%v_pre
  pb%v_pre=pb%v

  pb%ot%unit = 10000

  do i=1,pb%mesh%nn
    if (pb%iot(i) == 1) then
      pb%ot%unit = pb%ot%unit+1
      write(pb%ot%unit,'(e24.16,5e14.6)') pb%time, pb%v(i),      &
      pb%theta(i), pb%tau(i), pb%slip(i), pb%sigma(i)
    endif
  enddo

endif

end subroutine ot_write

!=====================================================================
! Export snapshots
subroutine ox_write(pb)

  use problem_class
  use constants, only: OUT_MASTER
  use my_mpi, only: is_MPI_parallel, is_mpi_master, my_mpi_tag, synchronize_all

  type (problem_type), intent(inout) :: pb

  integer :: ixout,i
  double precision :: vtempglob
  character(len=256) :: fileproc

if (is_MPI_parallel()) then
! In progress
! if (mod(pb%it-1,pb%ot%ntout) == 0 .or. pb%it == pb%itstop) then
 if (mod(pb%it,pb%ot%ntout) == 0 .or. pb%it == pb%itstop) then
  if (OUT_MASTER) then
  ! Collecting global nodes
    call synchronize_all()
    call pb_global(pb)
    if (is_mpi_master()) then
      vtempglob=0d0
      do i=1,pb%mesh%nnglob
        if ( pb%v_glob(i) > vtempglob) then
          vtempglob = pb%v_glob(i)
          pb%ot%ivmaxglob = i
        end if
      end do
      ! Writing fault points in single file fort.19
      if (pb%ox%i_ox_seq == 0) then
        write(pb%ox%unit,'(2a,2i8,e14.6)') '# x y z t v theta',' V_dot/V dtau tau_dot slip sigma ',&
                                          pb%it,pb%ot%ivmaxglob,pb%time
        do ixout=1,pb%mesh%nnglob,pb%ox%nxout
          write(pb%ox%unit,'(3e15.7,e24.16,7e15.7)') pb%mesh%xglob(ixout),pb%mesh%yglob(ixout),&
          pb%mesh%zglob(ixout),pb%time,pb%v_glob(ixout),pb%theta_glob(ixout),&
          pb%dv_dt_glob(ixout)/pb%v_glob(ixout),&
          pb%tau_glob(ixout),pb%dtau_dt_glob(ixout),pb%slip_glob(ixout), pb%sigma_glob(ixout)
        enddo
      else
      !Writing in File output in snapshots. fort.XXXXXX
      pb%ox%unit = pb%ox%unit + 1
      write(pb%ox%unit,'(3i10,e24.14)') pb%it,pb%ot%ivmaxglob,pb%ox%countglob,pb%time
      write(pb%ox%unit,'(2a)') '#  x  y  z  t  v  theta','  V./V  dtau  tau_dot  slip '
      do ixout=1,pb%mesh%nnglob,pb%ox%nxout
        write(pb%ox%unit,'(3e15.7,e24.14,7e15.7)')       &
          pb%mesh%xglob(ixout),pb%mesh%yglob(ixout),pb%mesh%zglob(ixout),pb%time,     &
          pb%v_glob(ixout),pb%theta_glob(ixout),pb%dv_dt_glob(ixout)/pb%v_glob(ixout),pb%tau_glob(ixout),   &
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
      write(pb%ox%unit,'(2a)') '#  x  y  z  t  v  theta','  V./V  dtau  tau_dot  slip '
      do ixout=1,pb%mesh%nn,pb%ox%nxout
        write(pb%ox%unit,'(3e15.7,e24.14,7e15.7)')       &
          pb%mesh%x(ixout),pb%mesh%y(ixout),pb%mesh%z(ixout),pb%time,     &
          pb%v(ixout),pb%theta(ixout),pb%dv_dt(ixout)/pb%v(ixout),pb%tau(ixout),   &
          pb%dtau_dt(ixout),pb%slip(ixout), pb%sigma(ixout)
      enddo
      close(pb%ox%unit) !Closing snapshot
  endif
!  close(pb%ox%unit)
 endif

 if (pb%ox%i_ox_dyn == 1) then
   if (OUT_MASTER) then
     call synchronize_all()
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
        write(20001+3*pb%ox%dyn_count2,'(2a)') '#  x  y  z  t  v  theta',  &
              '  V./V  dtau  tau_dot  slip sigma'
        do ixout=1,pb%mesh%nnglob,pb%ox%nxout_dyn
           write(20001+3*pb%ox%dyn_count2,'(3e15.7,e24.14,7e15.7)')       &
           pb%mesh%xglob(ixout),pb%mesh%yglob(ixout),pb%mesh%zglob(ixout),pb%time,     &
           pb%v_glob(ixout),pb%theta_glob(ixout),pb%dv_dt_glob(ixout)/pb%v_glob(ixout),pb%tau_glob(ixout),   &
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
        write(20002+3*pb%ox%dyn_count2,'(2a)') '#  x  y  z  t  v  theta',  &
              '  V./V  dtau  tau_dot  slip sigma'
        do ixout=1,pb%mesh%nnglob,pb%ox%nxout_dyn
           write(20002+3*pb%ox%dyn_count2,'(3e15.7,e24.14,7e15.7)')       &
           pb%mesh%xglob(ixout),pb%mesh%yglob(ixout),pb%mesh%zglob(ixout),pb%time,     &
           pb%v_glob(ixout),pb%theta_glob(ixout),pb%dv_dt_glob(ixout)/pb%v_glob(ixout),pb%tau_glob(ixout),   &
           pb%dtau_dt_glob(ixout),pb%slip_glob(ixout),pb%sigma_glob(ixout)
        enddo
        close(20002+3*pb%ox%dyn_count2)

        write(20003+3*pb%ox%dyn_count2,'(2a)') '#  x  y  z  t_rup tau_max t_vmax vmax'
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
    ! SEISMIC: added time to header
    write(pb%ox%unit,'(2a,2i8,e14.6)')'# x t v theta',' V./V dtau tau_dot slip ',pb%it,pb%ot%ivmax,pb%time
  ! JPA: this output should also contain y and z
    do ixout=1,pb%mesh%nn,pb%ox%nxout
      write(pb%ox%unit,'(e15.7,e24.16,7e15.7)') pb%mesh%x(ixout),pb%time,pb%v(ixout),   &
        pb%theta(ixout),pb%dv_dt(ixout)/pb%v(ixout),pb%tau(ixout),   &
        pb%dtau_dt(ixout),pb%slip(ixout), pb%sigma(ixout)
    enddo
  else
    pb%ox%unit = pb%ox%unit + 1
    write(pb%ox%unit,'(3i10,e24.14)') pb%it,pb%ot%ivmax,pb%ox%count,pb%time
    write(pb%ox%unit,'(2a)') '#  x  y  z  t  v  theta','  V./V  dtau  tau_dot  slip '
    do ixout=1,pb%mesh%nn,pb%ox%nxout
      write(pb%ox%unit,'(3e15.7,e24.14,7e15.7)')       &
        pb%mesh%x(ixout),pb%mesh%y(ixout),pb%mesh%z(ixout),pb%time,     &
        pb%v(ixout),pb%theta(ixout),pb%dv_dt(ixout)/pb%v(ixout),pb%tau(ixout),   &
        pb%dtau_dt(ixout),pb%slip(ixout), pb%sigma(ixout)
    enddo
    close(pb%ox%unit)
  endif
 endif

  if (pb%ox%i_ox_dyn == 1) then

    if (pb%ox%dyn_stat2 == 0 .and. pb%v(pb%ot%ivmax) >= pb%DYN_th_on ) then
      pb%ox%dyn_stat2 = 1
      do ixout=1,pb%mesh%nn,pb%ox%nxout_dyn
        pb%tau_max(ixout) = pb%tau(ixout)
        pb%t_rup(ixout) = pb%time
        pb%v_max(ixout) = pb%v(ixout)
        pb%t_vmax(ixout) = pb%time
      enddo
      write(20001+3*pb%ox%dyn_count2,'(3i10,e24.14)')   &
            pb%it,pb%ot%ivmax,pb%ox%count,pb%time
      write(20001+3*pb%ox%dyn_count2,'(2a)') '#  x  y  z  t  v  theta',  &
            '  V./V  dtau  tau_dot  slip sigma'
      do ixout=1,pb%mesh%nn,pb%ox%nxout_dyn
        write(20001+3*pb%ox%dyn_count2,'(3e15.7,e24.14,7e15.7)')       &
          pb%mesh%x(ixout),pb%mesh%y(ixout),pb%mesh%z(ixout),pb%time,     &
          pb%v(ixout),pb%theta(ixout),pb%dv_dt(ixout)/pb%v(ixout),pb%tau(ixout),   &
          pb%dtau_dt(ixout),pb%slip(ixout), pb%sigma(ixout)
      enddo
      close(20001+3*pb%ox%dyn_count2)
    endif

    if (pb%ox%dyn_stat2 == 1) then
      do ixout=1,pb%mesh%nn,pb%ox%nxout_dyn
        if (pb%tau(ixout) > pb%tau_max(ixout)) then
          pb%tau_max(ixout) = pb%tau(ixout)
          pb%t_rup(ixout) = pb%time
        endif
        if (pb%v(ixout) > pb%v_max(ixout)) then
          pb%v_max(ixout) = pb%v(ixout)
          pb%t_vmax(ixout) = pb%time
        endif
      enddo
    endif

    if (pb%ox%dyn_stat2 == 1 .and. pb%v(pb%ot%ivmax) <= pb%DYN_th_off ) then
      pb%ox%dyn_stat2 = 0
      write(20002+3*pb%ox%dyn_count2,'(3i10,e24.14)')   &
            pb%it,pb%ot%ivmax,pb%ox%count,pb%time
      write(20002+3*pb%ox%dyn_count2,'(2a)') '#  x  y  z  t  v  theta',  &
            '  V./V  dtau  tau_dot  slip sigma'
      do ixout=1,pb%mesh%nn,pb%ox%nxout_dyn
        write(20002+3*pb%ox%dyn_count2,'(3e15.7,e24.14,7e15.7)')       &
          pb%mesh%x(ixout),pb%mesh%y(ixout),pb%mesh%z(ixout),pb%time,     &
          pb%v(ixout),pb%theta(ixout),pb%dv_dt(ixout)/pb%v(ixout),pb%tau(ixout),   &
          pb%dtau_dt(ixout),pb%slip(ixout),pb%sigma(ixout)
      enddo
      close(20002+3*pb%ox%dyn_count2)

      write(20003+3*pb%ox%dyn_count2,'(2a)') '#  x  y  z  t_rup tau_max t_vmax vmax'
      do ixout=1,pb%mesh%nn,pb%ox%nxout_dyn
        write(20003+3*pb%ox%dyn_count2,'(3e15.7,4e28.20)')       &
          pb%mesh%x(ixout),pb%mesh%y(ixout),pb%mesh%z(ixout),       &
          pb%t_rup(ixout),pb%tau_max(ixout),pb%t_vmax(ixout),pb%v_max(ixout)
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
      write(100,'(2a)') '#  x  y  z  t  v  theta','  V./V  dtau  tau_dot  slip sigma '
      do ixout=1,pb%mesh%nn,pb%ox%nxout_dyn
        write(pb%ox%unit,'(3e15.7,e24.14,6e15.7)')       &
          pb%mesh%x(ixout),pb%mesh%y(ixout),pb%mesh%z(ixout),pb%time,     &
          pb%v(ixout),pb%theta(ixout),pb%dv_dt(ixout)/pb%v(ixout),pb%tau(ixout),   &
          pb%dtau_dt(ixout),pb%slip(ixout),pb%sigma(ixout)
      enddo
      pb%pot_pre = pb%pot
      CLOSE(100)
    endif

    if (pb%ox%dyn_stat == 1 .and. pb%v(pb%ot%ivmax) <= pb%DYN_th_off ) then
      pb%ox%dyn_stat = 0
      OPEN (UNIT = 101, FILE='DYN_POST.txt', STATUS='REPLACE')
      write(101,'(3i10,e24.14)') pb%it,pb%ot%ivmax,pb%ox%count,pb%time
      write(101,'(2a)') '#  x  y  z  t  v  theta','  V./V  dtau  tau_dot  slip sigma'
      do ixout=1,pb%mesh%nn,pb%ox%nxout_dyn
        write(pb%ox%unit,'(3e15.7,e24.14,6e15.7)')       &
          pb%mesh%x(ixout),pb%mesh%y(ixout),pb%mesh%z(ixout),pb%time,     &
          pb%v(ixout),pb%theta(ixout),pb%dv_dt(ixout)/pb%v(ixout),pb%tau(ixout),   &
          pb%dtau_dt(ixout),pb%slip(ixout),pb%sigma(ixout)
      enddo
      CLOSE(101)
      if ((pb%pot-pb%pot_pre)*pb%smu >= pb%DYN_M) then
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

  if (n == 0) then
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

  use fault_stress, only: nnLocal_perproc,nnoffset_glob_perproc
  use problem_class
  use my_mpi, only: my_mpi_rank, gather_allvdouble_root

  type(problem_type), intent(inout) :: pb
  integer :: nLocal,nnGlobal

  nLocal=nnLocal_perproc(my_mpi_rank())
  nnGlobal=sum(nnLocal_perproc)

  if (.not.allocated(pb%v_glob)) then

    allocate(pb%v_glob(nnGlobal),pb%theta_glob(nnGlobal),pb%tau_glob(nnGlobal),&
           pb%slip_glob(nnGlobal),pb%sigma_glob(nnGlobal),pb%dv_dt_glob(nnGlobal),&
           pb%dtheta_dt_glob(nnGlobal),pb%dtau_dt_glob(nnGlobal))

    allocate(pb%tau_max_glob(nnGlobal),pb%t_rup_glob(nnGlobal),&
             pb%v_max_glob(nnGlobal),pb%t_vmax_glob(nnGlobal))

  endif

  pb%v_glob=0
  pb%theta_glob=0
  pb%tau_glob=0
  pb%slip_glob=0
  pb%sigma_glob=0
  pb%dv_dt_glob=0
  pb%dtheta_dt_glob=0
  pb%dtau_dt_glob=0

  pb%tau_max_glob=0
  pb%t_rup_glob=0
  pb%v_max_glob=0
  pb%t_vmax=0

  call gather_allvdouble_root(pb%v,nLocal,pb%v_glob,nnLocal_perproc, &
                           nnoffset_glob_perproc,nnGlobal)
  call gather_allvdouble_root(pb%dv_dt,nLocal,pb%dv_dt_glob,nnLocal_perproc, &
                           nnoffset_glob_perproc,nnGlobal)
  call gather_allvdouble_root(pb%theta,nLocal,pb%theta_glob,nnLocal_perproc, &
                           nnoffset_glob_perproc,nnGlobal)
  call gather_allvdouble_root(pb%dtheta_dt,nLocal,pb%dtheta_dt_glob,nnLocal_perproc, &
                           nnoffset_glob_perproc,nnGlobal)
  call gather_allvdouble_root(pb%tau,nLocal,pb%tau_glob,nnLocal_perproc, &
                           nnoffset_glob_perproc,nnGlobal)
  call gather_allvdouble_root(pb%dtau_dt,nLocal,pb%dtau_dt_glob,nnLocal_perproc, &
                           nnoffset_glob_perproc,nnGlobal)
  call gather_allvdouble_root(pb%slip,nLocal,pb%slip_glob,nnLocal_perproc, &
                           nnoffset_glob_perproc,nnGlobal)
  call gather_allvdouble_root(pb%sigma,nLocal,pb%sigma_glob,nnLocal_perproc, &
                           nnoffset_glob_perproc,nnGlobal)

  end subroutine pb_global
!=====================================================================


end module output
