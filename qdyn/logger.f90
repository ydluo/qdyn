module logger

  implicit none
  private

  public :: init_log, log_msg

! Custom operator definition
interface operator( .exists. )
    module procedure file_exists
end interface

contains

!=====================================================================
! Function to check if a given file exists
function file_exists(filename) result(res)
  implicit none
  character(len=*),intent(in) :: filename
  logical :: res
  inquire(file=trim(filename), exist=res)
end function
!=====================================================================

!=====================================================================
! Initiate the log
subroutine init_log()

  use constants, only : FID_LOG, FILE_LOG, RESTART, VERBOSE
  use my_mpi, only : is_mpi_master, my_mpi_NPROCS

  logical :: file_exists = .false.
  character(100) :: msg

  ! If this is not the master node: do nothing
  ! If user requests logging to screen, nothing needs to be done
  if (.not. is_MPI_master() .or. VERBOSE .eqv. .true.) return

  file_exists = .exists. FILE_LOG

  ! First we check if the file exists. 
  ! If not: create new file and our job is done
  if (.not. file_exists) then
    open(FID_LOG, file=FILE_LOG, status="new")
  ! If the file exists but we are not restarting the
  ! simulation: wipe existing file
  elseif(file_exists .and. .not. RESTART) then
    open(FID_LOG, file=FILE_LOG, status="replace")
  ! If the simulation is being restarted, append to existing log
  elseif(file_exists .and. RESTART) then
    open(FID_LOG, file=FILE_LOG, status="old", access="append")
  endif

  close(FID_LOG)

  write(msg, *) "Number of processors = ", my_mpi_NPROCS()
  call log_msg(msg)

end subroutine init_log
!=====================================================================

!=====================================================================
! Add a log entry
subroutine log_msg(msg)

  use constants, only : FID_LOG, FILE_LOG, VERBOSE
  use my_mpi, only : is_mpi_master

  character(*) :: msg

  ! If this is not the master node: do nothing
  if (.not. is_mpi_master()) return

  ! If verbose: write to screen
  if (VERBOSE) then
    write(FID_LOG, *) msg
  ! Else append to log file
  else
    open(FID_LOG, file=FILE_LOG, status="old", access="append")
    write(FID_LOG, *) trim(msg)
    close(FID_LOG)
  endif
  
end subroutine log_msg
!=====================================================================

end module logger