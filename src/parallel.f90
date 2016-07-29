!----------------------------------------------------------------------------
! MPI parallel routines for QDYN. 
! This file containes all the MPI routines used in QDYN to run in parallel.
! Modified from Parallel.f90 SPECFEM3D.

module my_mpi

  use mpi

  implicit none

  integer :: MY_RANK=0, NPROCS=1

  integer :: my_local_mpi_comm_world, my_local_mpi_comm_for_bcast

  public MY_RANK, NPROCS

end module my_mpi

!----------------------------------------------------------------------------
! Subroutine to initialize MPI 
subroutine init_mpi()

  use my_mpi
  
  implicit none
  
  integer :: ier 

  call MPI_INIT(ier)
  if (ier /= 0 ) stop 'Error initializing MPI'

end subroutine init_mpi 

!-------------------------------------------------------------------------------------------------
subroutine world_size(sizeval)

  use my_mpi

  implicit none

  integer,intent(out) :: sizeval

  ! local parameters
  integer :: ier

  call MPI_COMM_SIZE(my_local_mpi_comm_world,sizeval,ier)
  if (ier /= 0 ) stop 'Error getting MPI world size'

end subroutine world_size
!----------------------------------------------------------------------------
! Retrieve processor number
subroutine world_rank(rank)

  use my_mpi

  implicit none

  integer,intent(out) :: rank

  integer :: ier

  call MPI_COMM_RANK(my_local_mpi_comm_world,rank,ier)
  if (ier /= 0 ) stop 'Error getting MPI rank'

end subroutine world_rank

!-------------------------------------------------------------------------------------------------
!Gather all MPI.
subroutine gather_allv(sendbuf, scounts, recvbufall, recvcountsall, recvoffsetall,recvcountstotal, NPROC)

  use my_mpi
  use constants, only: CUSTOM_REAL

  implicit none
  
  integer :: NPROC
  integer :: scounts,recvcountstotal
  real(kind=CUSTOM_REAL), dimension(scounts) :: sendbuf
  real(kind=CUSTOM_REAL), dimension(recvcountstotal) :: recvbufall
  integer, dimension(0:NPROC-1) :: recvcountsall,recvoffsetall

  integer ier

  !PG: sending the each processor data to the corresponding index in the recvbufall array.
  
  call MPI_ALLGATHERV(sendbuf,scounts,MPI_DOUBLE_PRECISION,recvbufall,recvcountsall,&
                      recvoffsetall,MPI_DOUBLE_PRECISION,my_local_mpi_comm_world,ier)

end subroutine gather_allv

!-------------------------------------------------------------------------------------------------
!Gather all MPI.
subroutine gather_allvdouble(sendbuf, scounts, recvbufall, recvcountsall, recvoffsetall,recvcountstotal, NPROC)

  use my_mpi
  use constants

  implicit none
  
  integer :: NPROC
  integer :: scounts,recvcountstotal
  double precision, dimension(scounts) :: sendbuf
  double precision, dimension(recvcountstotal) :: recvbufall
  integer, dimension(0:NPROC-1) :: recvcountsall,recvoffsetall

  integer ier

  !PG: sending the each processor data to the corresponding index in the recvbufall array.
  
  call MPI_ALLGATHERV(sendbuf,scounts,MPI_DOUBLE_PRECISION,recvbufall,recvcountsall,&
                      recvoffsetall,MPI_DOUBLE_PRECISION,my_local_mpi_comm_world,ier)

end subroutine gather_allvdouble

!-------------------------------------------------------------------------------------------------
!Gather all MPI to root processor=0
subroutine gather_allvdouble_root(sendbuf, scounts, recvbufall, recvcountsall, recvoffsetall,recvcountstotal, NPROC)

  use my_mpi

  implicit none
  
  integer :: NPROC
  integer :: scounts,recvcountstotal
  double precision, dimension(scounts) :: sendbuf
  double precision, dimension(recvcountstotal) :: recvbufall
  integer, dimension(0:NPROC-1) :: recvcountsall,recvoffsetall

  integer ier

  !PG: sending the each processor data to the corresponding index in the recvbufall array.
  
  call MPI_GATHERV(sendbuf,scounts,MPI_DOUBLE_PRECISION,recvbufall,recvcountsall,&
                      recvoffsetall,MPI_DOUBLE_PRECISION,0,my_local_mpi_comm_world,ier)

end subroutine gather_allvdouble_root

!-------------------------------------------------------------------------------------------------
subroutine gather_alli(sendbuf, recvbuf, NPROC)

  use my_mpi

  implicit none

  integer :: NPROC
  integer :: sendbuf
  integer, dimension(0:NPROC-1) :: recvbuf

  integer :: ier

  call MPI_ALLGATHER(sendbuf,1,MPI_INTEGER, &
                  recvbuf,1,MPI_INTEGER, &
                  my_local_mpi_comm_world,ier)

end subroutine gather_alli

!-------------------------------------------------------------------------------------------------
subroutine synchronize_all()

  use my_mpi

  implicit none

  integer :: ier

  ! synchronizes MPI processes
  call MPI_BARRIER(my_local_mpi_comm_world,ier)
  if (ier /= 0 ) stop 'Error synchronize MPI processes'

end subroutine synchronize_all

!-------------------------------------------------------------------------------------------------
!
  subroutine sum_allreduce(buffer,countval)

  use my_mpi

  implicit none

  integer :: countval
  double precision, dimension(countval),intent(inout) :: buffer

  ! local parameters
  integer :: ier
  double precision,dimension(countval) :: send

  ! seems not to be supported on all kind of MPI implementations...
  !! DK DK: yes, I confirm, using MPI_IN_PLACE is tricky
  !! DK DK (see the answer at http://stackoverflow.com/questions/17741574/in-place-mpi-reduce-crashes-with-openmpi
  !! DK DK      for how to use it right)
  !call MPI_ALLREDUCE(MPI_IN_PLACE, buffer, countval, MPI_INTEGER, MPI_MAX, my_local_mpi_comm_world, ier)

  send(:) = buffer(:)

  call MPI_ALLREDUCE(send, buffer, countval, MPI_DOUBLE_PRECISION, MPI_SUM, my_local_mpi_comm_world, ier)
  if (ier /= 0) stop 'Allreduce to get max values failed.'

  end subroutine sum_allreduce

!-------------------------------------------------------------------------------------------------
!
  subroutine max_allproc(sendbuf, recvbuf)

  use my_mpi

  implicit none

  double precision :: sendbuf, recvbuf
  integer ier

  call MPI_ALLREDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION,MPI_MAX,my_local_mpi_comm_world,ier)

  end subroutine max_allproc

!-------------------------------------------------------------------------------------------------
