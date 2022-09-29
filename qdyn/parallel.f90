!----------------------------------------------------------------------------
! Wrapper to MPI parallel routines for QDYN
! Modified from SPECFEM3D's parallel.f90

module my_mpi

  use mpi

  implicit none
  private

  integer, save :: MY_RANK=0, NPROCS=1

  public :: init_mpi, finalize_mpi, &
            my_mpi_tag, my_mpi_rank, my_mpi_NPROCS, is_mpi_parallel, &
            is_mpi_master, gather_allv, gather_allvdouble, &
            gather_allvdouble_root, gather_allvi_root, gather_alli, &
            synchronize_all, sum_allreduce, max_allproc, min_allproc

contains

!----------------------------------------------------------------------------
! Subroutine to initialize MPI
! If only one processor is found, MPI is finalized
subroutine init_mpi()

  integer :: ier

  call MPI_INIT(ier)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, NPROCS, ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD, MY_RANK, ier)

  if (ier /= 0 ) stop 'Error initializing MPI'

  if (NPROCS<2) call MPI_FINALIZE(ier)
  if (MY_RANK==0) write(6,*) 'Number of processors = ',NPROCS

end subroutine init_mpi

!-------------------------------------------------------------------------------------------------
! Finalize mpi
subroutine finalize_mpi()

  integer :: ier

! do NOT remove the barrier here, it is critical in order for the failsafe mechanism to work fine when it is activated
  call MPI_BARRIER(MPI_COMM_WORLD,ier)

! stop all the MPI processes, and exit
  call MPI_FINALIZE(ier)
  if (ier /= 0) stop 'Error finalizing MPI'

end subroutine finalize_mpi


!----------------------------------------------------------------------------
! Make a 6-character text tag based on processor rank
! The tag is blank if this is not an MPI parallel run.
! Use trim(my_mpi_tag) to include rank in file names if one file per processor is needed
function my_mpi_tag() result(iprocnum)

  character(6) :: iprocnum

  if (NPROCS>1) then
    write(iprocnum,'(i6.6)') MY_RANK
  else
    iprocnum = ""
  endif

end function my_mpi_tag

!----------------------------------------------------------------------------
integer function my_mpi_rank()
  my_mpi_rank = MY_RANK
end function my_mpi_rank

!----------------------------------------------------------------------------
integer function my_mpi_NPROCS()
  my_mpi_NPROCS = NPROCS
end function my_mpi_NPROCS

!----------------------------------------------------------------------------
logical function is_mpi_parallel()
  is_mpi_parallel = (NPROCS>1)
end function is_mpi_parallel

!----------------------------------------------------------------------------
logical function is_mpi_master()
  is_mpi_master = (MY_RANK==0)
end function is_mpi_master

!----------------------------------------------------------------------------
!Gather all MPI.
subroutine gather_allv(sendbuf, scounts, recvbufall, recvcountsall, recvoffsetall,recvcountstotal)

  use constants, only: CUSTOM_REAL

  implicit none

  integer :: scounts,recvcountstotal
  real(kind=CUSTOM_REAL), dimension(scounts) :: sendbuf
  real(kind=CUSTOM_REAL), dimension(recvcountstotal) :: recvbufall
  integer, dimension(0:NPROCS-1) :: recvcountsall,recvoffsetall

  integer ier

  !PG: sending the each processor data to the corresponding index in the recvbufall array.

  call MPI_ALLGATHERV(sendbuf,scounts,MPI_DOUBLE_PRECISION,recvbufall,recvcountsall,&
                      recvoffsetall,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ier)

end subroutine gather_allv

!-------------------------------------------------------------------------------------------------
!Gather all MPI.
subroutine gather_allvdouble(sendbuf, scounts, recvbufall, recvcountsall, recvoffsetall,recvcountstotal )

  use constants

  implicit none

  integer :: scounts,recvcountstotal
  double precision, dimension(scounts) :: sendbuf
  double precision, dimension(recvcountstotal) :: recvbufall
  integer, dimension(0:NPROCS-1) :: recvcountsall,recvoffsetall

  integer ier

  !PG: sending the each processor data to the corresponding index in the recvbufall array.

  call MPI_ALLGATHERV(sendbuf,scounts,MPI_DOUBLE_PRECISION,recvbufall,recvcountsall,&
                      recvoffsetall,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ier)

end subroutine gather_allvdouble

!-------------------------------------------------------------------------------------------------
!Gather all MPI to root processor=0
subroutine gather_allvdouble_root(sendbuf, scounts, recvbufall, recvcountsall, recvoffsetall,recvcountstotal)

  integer :: scounts,recvcountstotal
  double precision, dimension(scounts) :: sendbuf
  double precision, dimension(recvcountstotal) :: recvbufall
  integer, dimension(0:NPROCS-1) :: recvcountsall,recvoffsetall

  integer ier

  !PG: sending the each processor data to the corresponding index in the recvbufall array.

  call MPI_GATHERV(sendbuf,scounts,MPI_DOUBLE_PRECISION,recvbufall,recvcountsall,&
                      recvoffsetall,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)

end subroutine gather_allvdouble_root

!-------------------------------------------------------------------------------------------------
!Gather all MPI to root processor=0
subroutine gather_allvi_root(sendbuf, scounts, recvbufall, recvcountsall, recvoffsetall,recvcountstotal)

  integer :: scounts, recvcountstotal
  integer, dimension(scounts) :: sendbuf
  integer, dimension(recvcountstotal) :: recvbufall
  integer, dimension(0:NPROCS-1) :: recvcountsall, recvoffsetall

  integer ier

  !PG: sending the each processor data to the corresponding index in the recvbufall array.

  call MPI_GATHERV( sendbuf, scounts, MPI_INTEGER, recvbufall, recvcountsall,&
                    recvoffsetall, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)

end subroutine gather_allvi_root

!-------------------------------------------------------------------------------------------------
subroutine gather_alli(sendbuf, recvbuf)

  integer :: sendbuf
  integer, dimension(0:NPROCS-1) :: recvbuf

  integer :: ier

  call MPI_ALLGATHER(sendbuf,1,MPI_INTEGER, &
                  recvbuf,1,MPI_INTEGER, &
                  MPI_COMM_WORLD,ier)

end subroutine gather_alli

!-------------------------------------------------------------------------------------------------
subroutine synchronize_all()

  integer :: ier

  ! synchronizes MPI processes
  call MPI_BARRIER(MPI_COMM_WORLD,ier)
  if (ier /= 0 ) stop 'Error synchronize MPI processes'

end subroutine synchronize_all

!-------------------------------------------------------------------------------------------------
!
  subroutine sum_allreduce(buffer,countval)

  integer :: countval
  double precision, dimension(countval),intent(inout) :: buffer

  ! local parameters
  integer :: ier
  double precision,dimension(countval) :: send

  ! Comment in memory of Dimitri Komatitsch...
  !
  ! seems not to be supported on all kind of MPI implementations...
  !! DK DK: yes, I confirm, using MPI_IN_PLACE is tricky
  !! DK DK (see the answer at http://stackoverflow.com/questions/17741574/in-place-mpi-reduce-crashes-with-openmpi
  !! DK DK      for how to use it right)
  !call MPI_ALLREDUCE(MPI_IN_PLACE, buffer, countval, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ier)
  !

  send(:) = buffer(:)

  call MPI_ALLREDUCE(send, buffer, countval, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ier)
  if (ier /= 0) stop 'Allreduce to get summed values failed.'

  end subroutine sum_allreduce

!-------------------------------------------------------------------------------------------------
!
  subroutine max_allproc(sendbuf, recvbuf)

  double precision :: sendbuf, recvbuf
  integer ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ier)

  end subroutine max_allproc

!-------------------------------------------------------------------------------------------------
!
  subroutine min_allproc(sendbuf, recvbuf)

  double precision :: sendbuf, recvbuf
  integer ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ier)

  end subroutine min_allproc

!-------------------------------------------------------------------------------------------------

end module my_mpi
