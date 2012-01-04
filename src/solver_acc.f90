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



end module solver_acc
