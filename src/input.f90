!qdyn input all parameters

module input

  implicit none
  private 

  public   read_main




  integer :: iinput
  
contains
!=====================================================================
! read in all parameters
! 
  subroutine read_main(pb)
  
  use problem_class
  
  type(problem_type), intent(inout)  :: pb
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
    read(15,*)pb%ot%ntout, pb%ot%ic, pb%ox%nxout
    read(15,*)pb%beta, pb%smu
    read(15,*)pb%Tper, pb%Aper
    read(15,*)pb%dt_try, pb%dt_max, &
              pb%tmax, pb%acc
    read(15,*)pb%NSTOP

    allocate (pb%tau(pb%mesh%nn), pb%dtau_dt(pb%mesh%nn),    &
             pb%tau_init(pb%mesh%nn), pb%sigma(pb%mesh%nn), &
             pb%slip(pb%mesh%nn), pb%v(pb%mesh%nn), pb%dv_dt(pb%mesh%nn), &
             pb%theta(pb%mesh%nn),  pb%dtheta_dt(pb%mesh%nn),  &
             pb%a(pb%mesh%nn), pb%b(pb%mesh%nn), pb%dc(pb%mesh%nn),   &
             pb%mesh%x(pb%mesh%nn),  &
             pb%v1(pb%mesh%nn), pb%v2(pb%mesh%nn), pb%mu_star(pb%mesh%nn),& 
             pb%v_star(pb%mesh%nn), pb%theta_star(pb%mesh%nn))
 
    do iinput=1,pb%mesh%nn
      read(15,*)pb%sigma(iinput), pb%v(iinput), pb%theta(iinput),  &
                pb%a(iinput), pb%b(iinput), pb%dc(iinput), pb%v1(iinput), &
                pb%v2(iinput), pb%mu_star(iinput), pb%v_star(iinput)                 
    end do


  end if

  close(15)
  write(6,*) 'Input complete'


  end subroutine read_main


end module input
