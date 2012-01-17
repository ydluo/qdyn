! Solve_master

module solve_master

  implicit none
 
  private

  public  :: solve 

contains

!=====================================================================
! Master Solver    
!  
subroutine solve(pb)
  
  use problem_class
  use solver_acc
  use output, only : screen_init, screen_write, ox_write, ot_write
  
  type(problem_type), intent(inout)  :: pb

  call screen_init(pb)

  ! Time loop
  do while (pb%it /= pb%itstop)

    pb%it = pb%it + 1
    call do_bsstep(pb)
    if (pb%it == 1) then
      write(6,*) pb%dtau_dt
    end if
    call update_field(pb)
    call ot_write(pb)
    call check_stop(pb)   ! here itstop will change
    !--------Output onestep to screen and ox file(snap_shot)
    if(mod(pb%it-1,pb%ot%ntout) == 0 .or. pb%it == pb%itstop) then
      call screen_write(pb)
      call ox_write(pb)
    endif

  enddo

end subroutine solve



end module solve_master

