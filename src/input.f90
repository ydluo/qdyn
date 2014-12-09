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
  use mesh, only : read_mesh, mesh_get_size
  
  type(problem_type), intent(inout)  :: pb

  integer :: i,n
  
  write(6,*) 'Start reading input: ...'

  open(unit=15,FILE= 'qdyn.in') 

  call read_mesh(15,pb%mesh)
  n = mesh_get_size(pb%mesh)

  write(6,*) '   Mesh input complete'
  pb%kernel%kind = pb%mesh%dim+1 
  if (pb%mesh%nx < 4 .and. pb%mesh%dim==2) then
    write(6,*) 'nx < 4, FFT disabled'
    pb%kernel%kind = pb%kernel%kind+1
  end if      
     
  if (pb%kernel%kind==2) then 
    allocate(pb%kernel%k2f)
    read(15,*) pb%kernel%k2f%finite
  elseif (pb%kernel%kind==3) then 
    allocate(pb%kernel%k3f)
  elseif (pb%kernel%kind==4) then     ! 
    allocate(pb%kernel%k3)
  elseif (pb%kernel%kind==5) then
    allocate(pb%kernel%k3f2)
  end if
   
  read(15,*) pb%itheta_law
  read(15,*) pb%i_rns_law
  read(15,*) pb%kernel%i_sigma_cpl
  read(15,*) pb%neqs 

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

  read(15,*)pb%ot%ntout, pb%ot%ic, pb%ox%nxout, pb%ox%nxout_dyn,    &
            pb%ox%i_ox_seq, pb%ox%i_ox_dyn
  read(15,*)pb%beta, pb%smu, pb%lam, pb%v_th

!YD This part we may want to modify it later to be able to
!impose more complicated loading/pertubation
!functions involved: problem_class/problem_type; input/read_main 
!                    initialize/init_field;  derivs_all/derivs 

  read(15,*)pb%Tper, pb%Aper
  read(15,*)pb%dt_try, pb%dt_max,pb%tmax, pb%acc
  read(15,*)pb%NSTOP
  read(15,*)pb%DYN_FLAG,pb%DYN_SKIP
  read(15,*)pb%DYN_M,pb%DYN_th_on,pb%DYN_th_off
  write(6,*) '  Flags input complete'

!JPA some of these arrays should be allocated elsewhere, unless it's better to do it
!    here to optimize memory access
  allocate ( pb%tau(n),     &
             pb%dtau_dt(n), pb%dsigma_dt(n), &
             pb%tau_init(n), pb%sigma(n), &
             pb%slip(n), pb%v(n), pb%dv_dt(n), &
             pb%theta(n),  pb%dtheta_dt(n),  &
             pb%a(n), pb%b(n), pb%dc(n),   &
             pb%v_pre(n), pb%v_pre2(n), &
             pb%t_rup(n), pb%tau_max(n),   &
             pb%v_max(n), pb%t_vmax(n),   &
             pb%v1(n), pb%v2(n), pb%mu_star(n),& 
             pb%v_star(n), pb%theta_star(n),   &
             pb%iot(n),pb%iasp(n))
 
  do i=1,n
    read(15,*)pb%sigma(i), pb%v(i), pb%theta(i),  &
              pb%a(i), pb%b(i), pb%dc(i), pb%v1(i), &
              pb%v2(i), pb%mu_star(i), pb%v_star(i), pb%iot(i), pb%iasp(i)                 
  end do

  close(15)
  write(6,*) 'Input complete'

end subroutine read_main


end module input
