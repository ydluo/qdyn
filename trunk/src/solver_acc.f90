! functions for solver

module solver_acc

  implicit none

  private 

  public :: stop_check

contains

!=====================================================================
! stop_check: 
!
subroutine stop_check(pb,it)

  use problem_class
  use output, only : time_write,
  type(problem_type), intent(inout) :: pb
  integer :: it
  double precision :: vmax_old = 0d0, vmax_older = 0d0
  save vmax_old, vmax_older

  if (pb%itstop == 0) then
      !         STOP soon after end of slip localization 
    if (pb%NSTOP == 1) then
      if (pb%output%llocnew > pb%output%llocold) pb%itstop=it+2*pb%output%ntout

      ! STOP soon after maximum slip rate
    elseif (pb%NSTOP == 2) then

      if (it > 2 .and. vmax_old > vmax_older .and. pb%v(pb%ot%ivmax) < vmax_old)  &
          pb%itstop = it+10*pb%output%ntout
      vmax_older = vmax_old
      vmax_old = pb%v(pb%ot%ivmax)

        !         STOP at a slip rate threshold
    elseif (pb%NSTOP == 3) then    
      if (pb%v(pb%oyt%ivmax) > pb%tmax) pb%itstop = it    !here tmax is threshhold velocity

        !         STOP if time > tmax
    else
      call time_write(pb)
      if (pb%tmax > 0.d0 .and. pb%time > pb%tmax) pb%itstop = it
    endif
  endif
    
end subroutine stop_check

!=====================================================================
! do bs_step, pack and unpack 
!
subroutine do_bsstep(pb)

  use problem_class
  type(problem_type), intent(inout) :: pb

  double precision, dimension(:), allocatable ::  yt, dydt, yt_scale

  allocate (yt(pb%neqs*pb%mesh%nn))
  allocate (dydt(pb%neqs*pb%mesh%nn))
  allocate (yt_scale(pb%neqs*pb%mesh%nn))


  !-------Pack v, theta into yt---------------------------------    
  yt(1:pb%neqs:) = pb%v
  yt(2:pb%neqs:) = pb%theta
  dydt(1:pb%neqs:) = pb%dv_dt
  dydt(2:pb%neqs:) = pb%dtheta_dt
  !-------Pack v, theta into yt--------------------------------- 

  ! One step
 
  !--------Call EXT routine bsstep [Bulirsch-Stoer Method] --------------
  !-------- 
  yt_scale=dabs(yt)+dabs(pb%dt_try*dydt)
  call bsstep(yt,dydt,pb%neqs*pb%mesh%nn,pb%time,pb%acc,yt_scale,derivs,pb)
  if (pb%dt_max >  0.d0) then
    pb%dt_try = min(pb%dt_next,pb%dt_max)
  else
    pb%dt_try = pb%dt_next
  endif

  !-------Unpack yt into v, theta--------------------------------- 
  pb%v = yt(1:pb%neqs:)
  pb%theta = yt(2:pb%neqs:)
  pb%dv_dt = dydt(1:pb%neqs:)
  pb%dtheta_dt = dydt(2:pb%neqs:) 
  !-------Unpack yt into v, theta--------------------------------- 


  deallocate(yt,dydt,yt_scale)
end subroutine do_bsstep


end module solver_acc
