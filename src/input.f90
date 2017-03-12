!qdyn input all parameters

module input

  implicit none
  private 
 
  public :: read_main

contains
!=====================================================================
! Read input file qdyn.in
! If MPI, there is one qdyn*.in file per processor
! 
subroutine read_main(pb)
  
  use problem_class
  use mesh, only : read_mesh, mesh_get_size
  use my_mpi, only: my_mpi_tag, is_MPI_parallel
 
  type(problem_type), intent(inout)  :: pb

  integer :: i,n,nsta,ista,ik,i_sigma_cpl
  double precision :: xsta, ysta, zsta, dmin2, d2
  
  write(6,*) 'Start reading input ...'

  open(unit=15, file='qdyn'//trim(my_mpi_tag())//'.in', action='read')
  call read_mesh(15,pb%mesh)
  write(6,*) '   Mesh input complete'

  if (pb%mesh%dim==1) read(15,*) pb%finite
  read(15,*) pb%itheta_law
  read(15,*) pb%i_rns_law
  read(15,*) pb%i_sigma_cpl
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
!--?? Leave AS IS till we complete benchmark this 2D version ??--

  read(15,*)pb%ot%ntout, pb%ot%ic, pb%ox%nxout, pb%ox%nxout_dyn,    &
            pb%ox%i_ox_seq, pb%ox%i_ox_dyn
  read(15,*)pb%beta, pb%smu, pb%lam, pb%D, pb%H, pb%v_th

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

!JPA some of these arrays should be allocated in initialize.f90, 
!    unless it's better to do it here to optimize memory access

  n = mesh_get_size(pb%mesh) ! number of nodes in this processor
  allocate ( pb%tau(n), pb%dtau_dt(n), &
             pb%sigma(n), &
             pb%slip(n), pb%v(n), pb%theta(n),  &
             pb%a(n), pb%b(n), pb%dc(n),   &
             pb%v_pre(n), pb%v_pre2(n), &
             pb%t_rup(n), pb%tau_max(n),   &
             pb%v_max(n), pb%t_vmax(n),   &
             pb%v1(n), pb%v2(n), pb%mu_star(n),& 
             pb%v_star(n), pb%theta_star(n),   &
             pb%iot(n),pb%iasp(n),pb%coh(n))
 
  do i=1,n
    read(15,*)pb%sigma(i), pb%v(i), pb%theta(i),  &
              pb%a(i), pb%b(i), pb%dc(i), pb%v1(i), &
              pb%v2(i), pb%mu_star(i), pb%v_star(i), &
              pb%iot(i), pb%iasp(i), pb%coh(i)                 
  end do

  if (is_MPI_parallel()) then 
   !JPA Move this to a separate input file, one per processor
   !JPA and read it in subroutine read_mesh
    allocate(pb%mesh%x(pb%mesh%nn), pb%mesh%y(pb%mesh%nn),& 
             pb%mesh%z(pb%mesh%nn), pb%mesh%dip(pb%mesh%nn))     
    pb%mesh%dx = pb%mesh%Lfault/pb%mesh%nx !Check if dx is needed 
    do i=1,n
      read(15,*) pb%mesh%x(i),pb%mesh%y(i),pb%mesh%z(i),pb%mesh%dip(i)
    enddo

   ! Find station in this processor
   ! Keep only the last station found in this processor (ot%ic is not an array)
    if (pb%ot%ic == 0) then 
      dmin2 = 0.01d0*pb%mesh%dx ! distance tolerance = 1% grid size
      dmin2 = dmin2*dmin2
      open(unit=200,file='stations.dat',action='read')
      read(200,*) nsta
      do ista=1,nsta 
        read(200,*) xsta, ysta, zsta
        do ik=1,pb%mesh%nn
          d2=(pb%mesh%x(ik)-xsta)**2+(pb%mesh%y(ik)-ysta)**2+(pb%mesh%z(ik)-zsta)**2
          if (d2 < dmin2) then
            pb%ot%ic=ik
            write(6,*) 'Processor: ',my_mpi_tag(),', station ',ista, &
                       ' found, distance mismatch = ',d2
            exit
          endif
        enddo
      enddo
      close(200)
    endif

  endif

  close(15)
  write(6,*) 'Input complete'

end subroutine read_main


end module input
