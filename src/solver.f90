! Solve_master

module solver

  use problem_class, only : problem_type

  implicit none
  private
  
  integer(kind=8), save :: iktotal

  public  :: solve 

contains

!=====================================================================
! Master Solver    
!  
subroutine solve(pb)
  
  use output, only : screen_init, screen_write, ox_write, ot_write
  use my_mpi, only: is_MPI_parallel, is_mpi_master, finalize_mpi
  
  type(problem_type), intent(inout)  :: pb

  if (is_mpi_master()) call screen_init(pb)
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
!
subroutine do_bsstep(pb)

  use derivs_all
  use ode_bs

  type(problem_type), intent(inout) :: pb

  double precision, dimension(pb%neqs*pb%mesh%nn) :: yt, dydt, yt_scale
  integer :: ik

  ! Pack v, theta into yt
! yt(2::pb%neqs) = pb%v(pb%rs_nodes) ! JPA Coulomb
  yt(2::pb%neqs) = pb%v
  yt(1::pb%neqs) = pb%theta
  dydt(2::pb%neqs) = pb%dv_dt
  dydt(1::pb%neqs) = pb%dtheta_dt
  if ( pb%neqs == 3) then           ! Temp solution for normal stress coupling
    yt(3::pb%neqs) = pb%sigma
    dydt(3::pb%neqs) = pb%dsigma_dt
  endif

  ! this update of derivatives is only needed to set up the scaling (yt_scale)
  call derivs(pb%time,yt,dydt,pb)
  yt_scale=dabs(yt)+dabs(pb%dt_try*dydt)
  ! One step 
  call bsstep(yt,dydt,pb%neqs*pb%mesh%nn,pb%time,pb%dt_try,pb%acc,yt_scale,pb%dt_did,pb%dt_next,pb,ik)
!PG: Here is necessary a global min, or dt_next and dt_max is the same in all processors?.
  if (pb%dt_max >  0.d0) then
    pb%dt_try = min(pb%dt_next,pb%dt_max)
  else
    pb%dt_try = pb%dt_next
  endif
  iktotal=ik+iktotal
!  if (MY_RANK==0) write(6,*) 'iktotal=',iktotal,'pb%time=',pb%time
! Unpack yt into v, theta
!  pb%v(pb%rs_nodes) = yt(2::pb%neqs) ! JPA Coulomb
  pb%v = yt(2::pb%neqs)
  pb%theta = yt(1::pb%neqs)
  pb%dv_dt = dydt(2::pb%neqs)
  pb%dtheta_dt = dydt(1::pb%neqs) 
  if ( pb%neqs == 3) then           ! Temp solution for normal stress coupling
    pb%sigma = yt(3::pb%neqs)  
    pb%dsigma_dt = dydt(3::pb%neqs)  
  endif
  
end subroutine do_bsstep


!=====================================================================
! Update field: slip, tau, potency potency rate, crack,    
!
subroutine update_field(pb)
  
  use output, only : crack_size
  use friction, only : friction_mu
  use my_mpi, only: max_allproc, is_MPI_parallel

  type(problem_type), intent(inout) :: pb

  integer :: i,ix,iw  

  ! Update slip, stress. 
  pb%slip = pb%slip + pb%v*pb%dt_did
  pb%tau = pb%sigma * friction_mu(pb%v,pb%theta,pb) + pb%coh
  ! update potency and potency rate
  pb%pot=0d0;
  pb%pot_rate=0d0;
  if (pb%mesh%dim == 0 .or. pb%mesh%dim == 1) then
    pb%pot = sum(pb%slip) * pb%mesh%dx
    pb%pot_rate = sum(pb%v) * pb%mesh%dx
  else
    do iw=1,pb%mesh%nw
      do ix=1,pb%mesh%nx
        i=(iw-1)*pb%mesh%nx+ix
        pb%pot = pb%pot + pb%slip(i) * pb%mesh%dx * pb%mesh%dw(iw)
        pb%pot_rate = pb%pot_rate + pb%v(i) * pb%mesh%dx * pb%mesh%dw(iw)
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
! Finding global vmax
   call max_allproc(pb%v(pb%ot%ivmax),pb%vmaxglob)
!   if.not.(vtemp==vtempglob) pb%ot%ivmax=-1 !This processor does not host the maximum vel.
 endif

end subroutine update_field

!=====================================================================
! check stop: 
!
subroutine check_stop(pb)

  use output, only : time_write
  use my_mpi, only: is_MPI_parallel

  type(problem_type), intent(inout) :: pb

  double precision, save :: vmax_old = 0d0, vmax_older = 0d0

if (is_MPI_parallel()) then 
! In progress
  if (pb%itstop == -1) then
      !         STOP soon after end of slip localization 
    if (pb%NSTOP == 1) then
    !  if (pb%ot%llocnew > pb%ot%llocold) pb%itstop=pb%it+2*pb%ot%ntout

      ! STOP soon after maximum slip rate
    elseif (pb%NSTOP == 2) then

    !  if (pb%it > 2 .and. vmax_old > vmax_older .and. pb%v(pb%ot%ivmax) < vmax_old)  &
    !      pb%itstop = pb%it+10*pb%ot%ntout
    !  vmax_older = vmax_old
    !  vmax_old = pb%v(pb%ot%ivmax)

        !         STOP at a slip rate threshold
    elseif (pb%NSTOP == 3) then 
      if (pb%vmaxglob > pb%tmax) pb%itstop = pb%it    !here tmax is threshhold velocity
        !         STOP if time > tmax
    else
!      if (MY_RANK==0) call time_write(pb)
      if (pb%tmax > 0.d0 .and. pb%time > pb%tmax) pb%itstop = pb%it
    endif
  endif

else

  if (pb%itstop == -1) then
      !         STOP soon after end of slip localization 
    if (pb%NSTOP == 1) then
      if (pb%ot%llocnew > pb%ot%llocold) pb%itstop=pb%it+2*pb%ot%ntout

      ! STOP soon after maximum slip rate
    elseif (pb%NSTOP == 2) then

      if (pb%it > 2 .and. vmax_old > vmax_older .and. pb%v(pb%ot%ivmax) < vmax_old)  &
          pb%itstop = pb%it+10*pb%ot%ntout
      vmax_older = vmax_old
      vmax_old = pb%v(pb%ot%ivmax)

        !         STOP at a slip rate threshold
    elseif (pb%NSTOP == 3) then    
      if (pb%v(pb%ot%ivmax) > pb%tmax) pb%itstop = pb%it    !here tmax is threshhold velocity

        !         STOP if time > tmax
    else
      call time_write(pb)
      if (pb%tmax > 0.d0 .and. pb%time > pb%tmax) pb%itstop = pb%it
    endif
  endif

endif
    
end subroutine check_stop



end module solver
