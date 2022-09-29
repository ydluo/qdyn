module logger

    implicit none
    private

    public :: log_screen

contains

!=====================================================================
! Writes a line to screen
subroutine log_screen(msg)

    use constants, only : FID_SCREEN
    use my_mpi, only : is_mpi_master
  
    character(*) :: msg
  
    if (is_MPI_master()) then
      write(FID_SCREEN, *) msg
    endif
  
end subroutine log_screen
!=====================================================================

end module logger