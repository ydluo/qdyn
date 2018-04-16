module utils 

  implicit none
 !   double precision, allocatable, save :: Vec
 !   public :: Vec

contains
!
!=====================================================================
!
subroutine save_array(x,y,z,V,iproc,typ,nw,nx)

  double precision, dimension(nw*nx), intent(in) :: V
  double precision, dimension(nw*nx), intent(in) :: x,y,z
  character(len=256) :: fileproc
  character(len=16)   :: typ !PG, 'loc' or 'glo'

  integer :: i,j,iproc,nx,nw

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
