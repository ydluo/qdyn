! Solve_master

module solve_master

  implicit none
 
  private

  public   solve 
  

contains
  
  subroutine solve(pb)
  

  use solver_1d
  
  use problem_class
  
  type(problem_type) :: pb
  
  if (pb%mesh%kind == 0 .and. pb%kernel%kind == 0 .and. pb%kernel%k2f%kind == 0) then     ! solve 1D
      call solve_1D_fft(pb)
  end if

  end subroutine solve




end module solve_master
