! initialize all include parameters, kernel and fields

module initialize
  
  implicit none
  private

  public :: init_field

contains

!=============================================================
subroutine init_field(pb)
  
  use problem_class
  use constants, only : PI
  use output, only : ot_init, ox_init

  type(problem_type), intent(inout) :: pb

!YD This part we may want to modify it later to be able to
!impose more complicated loading/pertubation
!functions involved: problem_class/problem_type; input/read_main 
!                    initialize/init_field;  derivs_all/derivs
 
    !---------------------- dt_max & perturbation------------------
  if (pb%Aper /= 0.d0 .and. pb%Tper > 0.d0) then
    if (pb%dt_max > 0.d0) then
      pb%dt_max = min(pb%dt_max,0.2d0*pb%Tper)
    else
      pb%dt_max = 0.2d0*pb%Tper
    endif
  endif
  if (pb%Tper > 0.d0) then
    pb%Omper = 2.d0*PI/pb%Tper
  else
    pb%Omper = 0.d0
  endif
  !---------------------- dt_max & perturbation------------------

  !---------------------- impedance ------------------
  if (pb%beta > 0d0) then 
    pb%zimpedance = 0.5d0*pb%smu/pb%beta
  else
    pb%zimpedance = 0.d0
  endif
  write(6,*) 'impedance = ', pb%zimpedance
  !---------------------- impedance ------------------
     
  !---------------------- ref_value ------------------        
  pb%theta_star = pb%dc/pb%v2
  pb%tau_init = pb%sigma *   &
    (pb%mu_star- pb%a*log(pb%v1/pb%v+1d0)+ pb%b*log(pb%theta/pb%theta_star+1d0))
  pb%tau = pb%tau_init
  pb%slip = 0d0 
  !---------------------- ref_value ----------------- 

  !---------------------- init_value for solver ----------------- 
  pb%time = 0.d0
  pb%itstop = -1
  pb%it = 0
  !---------------------- init_value for solver ----------------- 
   
  call ot_init(pb)
  call ox_init(pb)
  write(6,*) 'Field Initialized'


end subroutine init_field


end module initialize
