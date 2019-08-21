! Collection of auxiliary functions

module utils

  use problem_class, only : problem_type
  implicit none
  public
 !   double precision, allocatable, save :: Vec
 !   public :: Vec

contains

!===============================================================================
! Helper routine to pack variables for solver
! Storage conventions:
!
! theta = yt(1::pb%neqs)
! v or tau = yt(2::pb%neqs) (depending on friction law)
! + feature specific variables
!
! dtheta/dt = dydt(1::pb%neqs)
! dv/dt or dtau/dt = dydt(2::pb%neqs)
! + feature specific variables
!===============================================================================
subroutine pack(yt, theta, main_var, sigma, theta2, pb)

  type(problem_type), intent(inout) :: pb
  double precision, dimension(pb%neqs*pb%mesh%nn), intent(out) :: yt
  double precision, dimension(pb%mesh%nn), intent(in) :: theta, main_var, sigma
  double precision, dimension(pb%mesh%nn), intent(in) :: theta2
  integer :: nmax, ind_stress_coupling, ind_localisation, ind_tp

  nmax = pb%neqs*pb%mesh%nn
  ! Define the indices of yt and dydt based on which
  ! features are requested (defined in input file)
  ind_stress_coupling = 2 + pb%features%stress_coupling
  ind_localisation = ind_stress_coupling + pb%features%localisation
  ind_tp = ind_localisation + pb%features%tp

  yt(1:nmax:pb%neqs) = theta
  yt(2:nmax:pb%neqs) = main_var

  if (pb%features%stress_coupling == 1) then
    yt(ind_stress_coupling:nmax:pb%neqs) = sigma
  endif

  if (pb%features%localisation == 1) then
    yt(ind_localisation:nmax:pb%neqs) = theta2
  endif

end subroutine pack

!===============================================================================
! Helper routine to unpack variables from solver
!===============================================================================
subroutine unpack(yt, theta, main_var, sigma, theta2, pb)

  type(problem_type), intent(inout) :: pb
  double precision, dimension(pb%neqs*pb%mesh%nn), intent(in) :: yt
  double precision, dimension(pb%mesh%nn) :: theta, main_var, sigma, theta2
  integer :: nmax, ind_stress_coupling, ind_localisation, ind_tp

  nmax = pb%neqs*pb%mesh%nn
  ind_stress_coupling = 2 + pb%features%stress_coupling
  ind_localisation = ind_stress_coupling + pb%features%localisation
  ind_tp = ind_localisation + pb%features%tp

  theta = yt(1:nmax:pb%neqs)
  main_var = yt(2:nmax:pb%neqs)

  if (pb%features%stress_coupling == 1) then
    sigma = yt(ind_stress_coupling:nmax:pb%neqs)
  else
    sigma = pb%sigma
  endif

  if (pb%features%localisation == 1) then
    theta2 = yt(ind_localisation:nmax:pb%neqs)
  else
    theta2 = 0d0
  endif

end subroutine unpack

!
!=====================================================================
!
subroutine save_array(x,y,z,V,iproc,typ,nw,nx)

  double precision, dimension(nw*nx), intent(in) :: V
  double precision, dimension(nw*nx), intent(in) :: x,y,z
  character(len=256) :: fileproc
  character(len=16)   :: typ !PG, 'loc' or 'glo'

  integer :: i,iproc,nx,nw

  write(fileproc,'(a,i6.6,a)') 'snap_v_',iproc,typ

  open(101,file=fileproc(1:len_trim(fileproc)),status='replace',form='formatted',action='write')

  do i=1,nw*nx
      write(101,'(4(D15.7))') x(i),y(i),z(i),V(i)
  enddo

  close(101)

end subroutine save_array

!------------------------------------------------


subroutine save_vectorV(x,y,z,V,iproc,typ,nw,nx)

  double precision, dimension(nw,nx), intent(in) :: V
  double precision, dimension(nw,nx), intent(in) :: x,y,z
  character(len=256) :: fileproc
  character(len=16)   :: typ !PG, 'loc' or 'glo'

  integer :: i,j,iproc,nx,nw

  write(fileproc,'(a,i6.6,a)') 'snap_v_',iproc,typ

  open(101,file=fileproc(1:len_trim(fileproc)),status='replace',form='formatted',action='write')

  do i=1,nw
    do j=1,nx
      write(101,'(4(D15.7))') x(i,j),y(i,j),z(i,j),V(i,j)
    enddo
  enddo

  close(101)

end subroutine save_vectorV

! -----------------------------------------------
subroutine save_vector(V,iproc,typ,nw,nx)

  double precision, dimension(nw,nx), intent(in) :: V
  character(len=256) :: fileproc
  character(len=16)   :: typ !PG, 'loc' or 'glo'

  integer :: i,j,iproc,nx,nw

  write(fileproc,'(a,i6.6,a)') 'snap_v_',iproc,typ

  open(101,file=fileproc(1:len_trim(fileproc)),status='replace',form='formatted',action='write')

  do i=1,nw
    do j=1,nx
      write(101,'(D15.7)') V(i,j)
    enddo
  enddo

  close(101)

end subroutine save_vector

! -----------------------------------------------
subroutine save_vector3(V,iproc,typ,nwloc,nwglob,nx)

  double precision, dimension(nwloc,nwglob,nx), intent(in) :: V
  character(len=256) :: fileproc
  character(len=16)   :: typ !PG, 'loc' or 'glo'

  integer :: i,j,k,iproc,nx,nwloc,nwglob

  write(fileproc,'(a,i6.6,a)') 'snap_v_',iproc,typ

  open(101,file=fileproc(1:len_trim(fileproc)),status='replace',form='formatted',action='write')

 do k=1,nwloc
  do i=1,nwglob
    do j=1,nx
      write(101,'(D15.7)') V(k,i,j)
    enddo
  enddo
 enddo

  close(101)

end subroutine save_vector3

end module utils
