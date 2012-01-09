!qdyn input all parameters

module input

  implicit none
  private 

  public :: read_main

contains
!=====================================================================
! read in all parameters
! 
subroutine read_main(pb)
  
  use problem_class
  
  type(problem_type), intent(inout)  :: pb

  integer :: i
  
  write(6,*) 'Start reading input: ...'

  open(unit=15,FILE= 'qdyn.in') 
  read(15,*)pb%mesh%kind

  if (pb%mesh%kind==0) then
    read(15,*)pb%mesh%nn
    read(15,*)pb%mesh%Lfault, pb%mesh%W 
   
    read(15,*)pb%kernel%kind
    if (pb%kernel%kind==0) then
      read(15,*) pb%kernel%k2f%finite
    end if
   
    read(15,*)pb%itheta_law
    read(15,*)pb%neqs 

!JPA neqs should not be setup explicitly by the user
!    It should be inferred from the type of problem:
!    neq=2 if problem in homogeneous medium without free surface
!    neq=3 if bimaterial problem, or with free surface (for which normal stress changes
!          are coupled to slip, the 3rd variable is the normal stress)
!YD we may want to determine the neqs and other value in matlab script, then read-in
! as far as it will not make conflict with other variables/parameters
! because matlab is using more human-like language 
!However, it will be safer to deal with variable/parameters here
!--?? Leave AS IS till we complete benchmark this 2D version ??---

    read(15,*)pb%ot%ntout, pb%ot%ic, pb%ox%nxout
    read(15,*)pb%beta, pb%smu

!YD This part we may want to modify it later to be able to
!impose more complicated loading/pertubation
!functions involved: problem_class/problem_type; input/read_main 
!                    initialize/init_field;  derivs_all/derivs 


    read(15,*)pb%Tper, pb%Aper
    read(15,*)pb%dt_try, pb%dt_max,pb%tmax, pb%acc
    read(15,*)pb%NSTOP

    allocate (pb%tau(pb%mesh%nn),     &
             pb%tau_init(pb%mesh%nn), pb%sigma(pb%mesh%nn), &
             pb%slip(pb%mesh%nn), pb%v(pb%mesh%nn), pb%dv_dt(pb%mesh%nn), &
             pb%theta(pb%mesh%nn),  pb%dtheta_dt(pb%mesh%nn),  &
             pb%a(pb%mesh%nn), pb%b(pb%mesh%nn), pb%dc(pb%mesh%nn),   &
             pb%mesh%x(pb%mesh%nn),  &
             pb%v1(pb%mesh%nn), pb%v2(pb%mesh%nn), pb%mu_star(pb%mesh%nn),& 
             pb%v_star(pb%mesh%nn), pb%theta_star(pb%mesh%nn))
 
    do i=1,pb%mesh%nn
      read(15,*)pb%sigma(i), pb%v(i), pb%theta(i),  &
                pb%a(i), pb%b(i), pb%dc(i), pb%v1(i), &
                pb%v2(i), pb%mu_star(i), pb%v_star(i)                 
    end do


  end if

  close(15)
  write(6,*) 'Input complete'

end subroutine read_main


end module input
