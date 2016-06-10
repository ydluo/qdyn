module utils 

  implicit none
 !   double precision, allocatable, save :: Vec
 !   public :: Vec

contains
!
!=====================================================================
!
subroutine save_vectorV(x,y,z,V,iproc,typ,nw,nx)

  double precision, dimension(nw,nx), intent(in) :: V
  double precision, dimension(nw,nx), intent(in) :: x,y,z
  character(len=256) :: filenameproc
  character(len=3)   :: typ !PG, 'loc' or 'glo'

  integer :: i,j,n,iproc,nx,nw

 ! nw=size(V,1)
 ! nx=size(V,2)

  write(filenameproc,'(a,a,i1,a)') 'snap_v_',typ,iproc,'.dat'
 
  open(101,file=trim(filenameproc),status='replace',form='formatted',action='write')
  
  do i=1,nw
    do j=1,nx
!      write(101,*) "x, y, z, data"
      write(101,'(4(D15.7))') x(i,j),y(i,j),z(i,j),V(i,j)  
    enddo
  enddo

  close(101)

end subroutine save_vectorV

end module utils
