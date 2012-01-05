! Solve_master

module solve_master

  implicit none
 
  private

  public   solve 
  

contains



!=====================================================================
! Master Solver    
!  
subroutine solve(pb)
  
  use problem_class
  use derivs_all
  use solver_acc
  use output
  
  type(problem_type), intent(inout)  :: pb
  integer :: it
  !=======================Time loop. START===========================
  ! Time loop

  do while (pb%it /= pb%itstop)

    pb%it = pb%it + 1
    call derivs(pb)
    call do_bsstep(pb)
    call update_field(pb)
    call ot_write(pb)
    call stop_check(pb)   ! here itstop will change
    !--------Output onestep to screen and ox file(snap_shot)
    if(mod(pb%it-1,pb%output%ntout) == 0 .or. pb%it == pb%itstop) then
      call screen_write(pb)
      call ox_write(pb)
    endif

  enddo

end subroutine solve



end module solve_master

