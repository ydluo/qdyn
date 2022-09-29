module constants

! By convention let's write the name of constants in capitals,
! but beware that Fortran does not make a difference between capitalized and non-capitalized names
! (i.e. pi and PI are the same)

!--- User Settings ---

! set BIN_OUTPUT = .true. for snapshot and time serie binary outputs.
! Only available for MESHDIM = 1, OX_SEQ = 0, OCTAVE_OUTPUT = .false., DYN_FLAG = 0
logical, parameter :: BIN_OUTPUT = .false.

! set usage of FFT in 3D
!   0 : no FFT
!   1 : FFT along-strike
!   2 : FFT along-strike and along-dip, only works for vertical faults
integer, parameter :: FFT_TYPE = 1

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

! The below "constants" are set via qdyn.in (see input.f90). Default values
! are set to 0 so that an error is thrown in case of an erroneous input file

! Type of faulting:
!   1 : strike-slip (right-lateral)
!  -1 : strike-slip (left-lateral)
!   2 : thrust
!  -2 : Normal
integer :: FAULT_TYPE = 0

! Which ODE solver to use:
!   1 : Bulirsch-Stoer
!   2 : Runge-Kutta-Fehlberg
integer :: SOLVER_TYPE = 0

! Input unit
integer, parameter :: FID_IN = 15

! Output units
integer, parameter :: FID_SCREEN = 6
integer, parameter :: FID_OT = 18
integer, parameter :: FID_OX = 19
integer, parameter :: FID_VMAX = 22
integer, parameter :: FID_IASP = 23
integer, parameter :: FID_QSB_PRE = 100
integer, parameter :: FID_QSB_POST = 101
integer, parameter :: FID_TIME = 121
integer, parameter :: FID_STATIONS = 200
integer, parameter :: FID_MW = 222
integer, parameter :: FID_OX_DYN = 20000

! Output names
character(*), parameter :: FILE_OX = "output_ox"
character(*), parameter :: FILE_OX_DYN_PRE = "output_dyn_pre_"
character(*), parameter :: FILE_OX_DYN_POST = "output_dyn_post_"
character(*), parameter :: FILE_OX_DYN_MAX = "output_dyn_max_"
character(*), parameter :: FILE_OT = "output_ot_"
character(*), parameter :: FILE_IASP = "output_iasp"
character(*), parameter :: FILE_VMAX = "output_vmax"

end module constants
