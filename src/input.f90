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
  use mesh, only : read_mesh_parameters, read_mesh_nodes, mesh_get_size
  use output, only : ot_read_stations
  use my_mpi, only: my_mpi_tag, is_MPI_parallel
 
  type(problem_type), intent(inout)  :: pb

  integer :: i,n
  
  write(6,*) 'Start reading input ...'

  open(unit=15, file='qdyn'//trim(my_mpi_tag())//'.in', action='read')

  call read_mesh_parameters(15,pb%mesh)
  write(6,*) '   Mesh input complete'

  if (pb%mesh%dim==1) read(15,*) pb%finite
  read(15,*) pb%itheta_law
  read(15,*) pb%i_rns_law
  read(15,*) pb%i_sigma_cpl
  read(15,*) pb%ot%ntout, pb%ot%ic, pb%ox%nxout, pb%ox%nxout_dyn,    &
             pb%ox%i_ox_seq, pb%ox%i_ox_dyn
  read(15,*) pb%beta, pb%smu, pb%lam, pb%D, pb%H, pb%v_th

!YD This part we may want to modify it later to be able to
!impose more complicated loading/pertubation
!functions involved: problem_class/problem_type; input/read_main 
!                    initialize/init_field;  derivs_all/derivs 

  read(15,*) pb%Tper, pb%Aper
  read(15,*) pb%dt_try, pb%dt_max,pb%tmax, pb%acc
  read(15,*) pb%NSTOP
  read(15,*) pb%DYN_FLAG,pb%DYN_SKIP
  read(15,*) pb%DYN_M,pb%DYN_th_on,pb%DYN_th_off
  write(6,*) '  Flags input complete'

  n = mesh_get_size(pb%mesh) ! number of nodes in this processor
  allocate ( pb%sigma(n), pb%v(n), pb%theta(n),  &
             pb%a(n), pb%b(n), pb%dc(n), pb%v1(n), &
             pb%v2(n), pb%mu_star(n), pb%v_star(n), &
             pb%iot(n),pb%iasp(n),pb%coh(n))
  do i=1,n
    read(15,*)pb%sigma(i), pb%v(i), pb%theta(i),  &
              pb%a(i), pb%b(i), pb%dc(i), pb%v1(i), &
              pb%v2(i), pb%mu_star(i), pb%v_star(i), &
              pb%iot(i), pb%iasp(i), pb%coh(i)                 
  end do

  if (is_MPI_parallel()) then 
    call read_mesh_nodes(15,pb%mesh)
    call ot_read_stations(pb%ot)
  endif

  close(15)
  write(6,*) 'Input complete'

end subroutine read_main


end module input
