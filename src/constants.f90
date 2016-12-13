module constants

! By convention let's write the name of constants in capitals,
! but beware that Fortran does not make a difference between capitalized and non-capitalized names
! (i.e. pi and PI are the same)

!--- User Settings ---

 ! set OCTAVE_OUTPUT = .true. if you are using Octave instead of Matlab
  logical, parameter :: OCTAVE_OUTPUT = .false.

 ! set the name of the kernel*.tab file, including its full path
  character(*), parameter :: KERNEL_FILE = "~/3D_RUPTURE/qdyn/trunk/src/kernel_I_32768.tab"

! set the type of faulting:
!   1 : strike-slip (right-lateral)
!  -1 : strike-slip (left-lateral)
!   2 : thrust
!  -2 : Normal
  integer, parameter :: FAULT_TYPE = 2

! set usage of FFT in 3D
!   0 : no FFT
!   1 : FFT along-strike
!   2 : FFT along-strike and along-dip, only works for vertical faults
  integer, parameter :: FFT_TYPE = 1

!  integer :: MY_RANK=0, NPROCS=1 !Default for serial 
                                 !but automatically changes if MPI_parallel=.true.

! MPI run in parallel
!   true  : run MPI parallel
!   false : run serial or openMP
   logical, parameter :: MPI_parallel = .false.
   logical, parameter :: OUT_MASTER = .true. !To write ouput with the master MY_RANK=0

! Adding real precision and type for MPI runs.
!   CUSTOM_REAL = 4  (single precision)
!   CUSTOM_REAL = 8  (double precision)
   integer, parameter :: CUSTOM_REAL = 8

!--- END of User Settings ---

  double precision, parameter :: PI = 3.141592653589793d0
  double precision, parameter :: DAY = 3600.d0*24.d0, &
                                 WEEK = 7d0*DAY, &
                                 MONTH = 30*DAY, &
                                 YEAR = 365*DAY

end module constants
