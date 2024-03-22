module logger

  implicit none
  private

  public :: init_log, log_msg, log_debug

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

  use constants, only : FID_LOG, FILE_LOG, DEBUG, RESTART, VERBOSE
  use my_mpi, only : is_mpi_master, my_mpi_NPROCS, my_mpi_rank, my_mpi_tag

  logical :: file_exists = .false.
  character(100) :: msg


  ! If this is not the master node: do nothing
  ! If user requests logging to screen, nothing needs to be done
  ! If we're in debugging mode, keep going
  if ((.not. is_MPI_master() .or. VERBOSE .eqv. .true.) .and. .not. DEBUG) return

  ! Redefine FILE_LOG and FID_LOG based on MPI tag
  write(msg, "(A,A)") trim(FILE_LOG), trim(my_mpi_tag())
  FILE_LOG = trim(msg)
  FID_LOG = FID_LOG + my_mpi_rank()

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

  use constants, only : FID_LOG, FILE_LOG, VERBOSE, DEBUG
  use my_mpi, only : is_mpi_master

  character(*) :: msg

  ! If this is not the master node: do nothing
  if (.not. is_mpi_master() .and. .not. DEBUG) return

  ! If verbose: write to screen
  if (VERBOSE .and. .not. DEBUG) then
    write(FID_LOG, *) msg
  ! Else append to log file
  else
    open(FID_LOG, file=FILE_LOG, status="old", access="append")
    write(FID_LOG, *) trim(msg)
    close(FID_LOG)
  endif
  
end subroutine log_msg
!=====================================================================

!=====================================================================
! Add a log entry
subroutine log_debug(msg, step)

  use constants, only : FID_LOG, FILE_LOG
  use my_mpi, only : is_mpi_master

  character(255) :: msg
  character(255) :: msg2
  integer :: step

  write(msg2, "(A,I0,A,A)") "[", step, "]", trim(msg)

  open(FID_LOG, file=FILE_LOG, status="old", access="append")
  write(FID_LOG, *) trim(msg2)
  close(FID_LOG)
  
end subroutine log_debug
!=====================================================================

end module logger