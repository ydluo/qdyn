module constants

! By convention let's write the name of constants in capitals,
! but beware that Fortran does not make a difference between capitalized and non-capitalized names
! (i.e. pi and PI are the same)

!--- User Settings ---

 ! set OCTAVE_OUTPUT = .true. if you are using Octave instead of Matlab
  logical, parameter :: OCTAVE_OUTPUT = .false.

! set the type of faulting:
!   1 : strike-slip (right-lateral)
!  -1 : strike-slip (left-lateral)
!   2 : thrust
!  -2 : Normal
  integer, parameter :: FAULT_TYPE = 1 
!JPA this should be an input in qdyn.in, not a parameter here

! set usage of FFT in 3D
!   0 : no FFT
!   1 : FFT along-strike
!   2 : FFT along-strike and along-dip, only works for vertical faults
  integer, parameter :: FFT_TYPE = 1

! For parallel MPI runs: to write global outputs only by the master processor
   logical, parameter :: OUT_MASTER = .true. 

! Adding real precision and type for MPI runs.
!   CUSTOM_REAL = 4  (single precision)
!   CUSTOM_REAL = 8  (double precision)
   integer, parameter :: CUSTOM_REAL = 8

!--- END of User Settings ---

 ! Path to directory containing the source files and kernel file.
 ! Needed to read the kernel file when FINITE=1.
 ! The variable _FPP_SRC_PATH_ must be defined at compile time
 ! using the preprocessor directive -D_FPP_SRC_PATH_="'...your_src_path...'"
 ! The Makefile sets this automatically as  -D_FPP_SRC_PATH_="'$(CURDIR)'"
  character(*), parameter :: SRC_PATH = _FPP_SRC_PATH_

  double precision, parameter :: PI = 3.141592653589793d0
  double precision, parameter :: DAY = 3600.d0*24.d0, &
                                 WEEK = 7d0*DAY, &
                                 MONTH = 30*DAY, &
                                 YEAR = 365*DAY

end module constants
