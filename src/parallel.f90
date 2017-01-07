!----------------------------------------------------------------------------
! MPI parallel routines for QDYN. 
! This file containes all the MPI routines used in QDYN to run in parallel.
! Modified from Parallel.f90 SPECFEM3D.

module my_mpi

  use mpi

  implicit none

  integer, public :: MY_RANK=0, NPROCS=1

contains

!----------------------------------------------------------------------------
! Subroutine to initialize MPI 
subroutine init_mpi()

  integer :: ier 

  call MPI_INIT(ier)
  if (ier /= 0 ) stop 'Error initializing MPI'
  call world_rank(MY_RANK)
  call world_size(NPROCS)

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

!-------------------------------------------------------------------------------------------------
subroutine world_size(sizeval)

  integer,intent(out) :: sizeval

  ! local parameters
  integer :: ier

  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeval,ier)
  if (ier /= 0 ) stop 'Error getting MPI world size'

end subroutine world_size

!----------------------------------------------------------------------------
! Retrieve processor number
subroutine world_rank(rank)

  integer,intent(out) :: rank

  integer :: ier

  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ier)
  if (ier /= 0 ) stop 'Error getting MPI rank'

end subroutine world_rank

!-------------------------------------------------------------------------------------------------
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

  ! seems not to be supported on all kind of MPI implementations...
  !! DK DK: yes, I confirm, using MPI_IN_PLACE is tricky
  !! DK DK (see the answer at http://stackoverflow.com/questions/17741574/in-place-mpi-reduce-crashes-with-openmpi
  !! DK DK      for how to use it right)
  !call MPI_ALLREDUCE(MPI_IN_PLACE, buffer, countval, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ier)

  send(:) = buffer(:)

  call MPI_ALLREDUCE(send, buffer, countval, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ier)
  if (ier /= 0) stop 'Allreduce to get max values failed.'

  end subroutine sum_allreduce

!-------------------------------------------------------------------------------------------------
!
  subroutine max_allproc(sendbuf, recvbuf)

  double precision :: sendbuf, recvbuf
  integer ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ier)

  end subroutine max_allproc

!-------------------------------------------------------------------------------------------------
end module my_mpi
