module mesh

  implicit none
  private

  type mesh_type
    integer :: dim = 0  ! dim = 1, 2 ,3 ~xD
    integer :: nx, nw, nn ! along-strike, along-dip, total grid number
    double precision :: dx !along-strike grid size(constant)  
    double precision :: Lfault, W, Z_CORNER ! fault length, width, lower-left corner z (follow Okada's convention)
    double precision, allocatable :: dw(:), DIP_W(:) !along-dip grid size and dip (adjustable), nw count
    double precision, allocatable :: x(:), y(:), z(:), dip(:) !coordinates and dip of every grid (nx*nw count)
    double precision, allocatable :: temp_xx(:), temp_ww(:)  !temp for Okada kernel output (nx*nw count)
    double precision :: temp_mu, temp_lam, temp_tau, temp_sigma_n !temp for Okada kernel output  
    integer :: temp_iret !temp for Okada kernel output 

  end type mesh_type

  public :: mesh_type, read_mesh, init_mesh, mesh_get_size

contains

!=============================================================
subroutine read_mesh(iin,m)

  use okada, only : compute_kernel

  type(mesh_type), intent(inout) :: m
  integer, intent(in) :: iin

  integer :: i,j

  ! problem dimension (1D, 2D or 3D), mesh type
  read(iin,*) m%dim

  !spring-block or 2d problem
  select case (m%dim)

  case(0,1)
    read(iin,*) m%nn
    read(iin,*) m%Lfault, m%W 

  case(2) !3d problem
    read(iin,*) m%nx,m%nw
    m%nn = m%nx * m%nw
    read(iin,*) m%Lfault, m%W , m%Z_CORNER   ! JPA m%W is not used in this case, remove it
    allocate(m%dw(m%nw), m%DIP_W(m%nw))
    do i=1,m%nw
      read(iin,*) m%dw(i), m%DIP_W(i)
    end do

 ! JPA The two cases below are not for earthquake cycle computations, 
 !     they should be in a separate code in utils/ linked to modules in src/
  case(67) 
    write(6,*) 'Calculate Okada Kernel'
    read(iin,*) m%nn
    read(iin,*) m%temp_lam, m%temp_mu
    allocate(m%x(m%nn),m%y(m%nn),m%z(m%nn),m%dip(m%nn))
    allocate(m%temp_xx(m%nn),m%temp_ww(m%nn))
    do i =1,m%nn
       read(iin,*) m%x(i),m%y(i),m%z(i),m%dip(i),m%temp_xx(i),m%temp_ww(i)
    end do
    do i=1,m%nn
      do j=1,m%nn
        call compute_kernel(m%temp_lam,m%temp_mu,m%x(i),m%y(i),m%z(i),  &
               m%dip(i),m%temp_xx(i),m%temp_ww(i),   &
               m%x(j),m%y(j),m%z(j),m%dip(j),m%temp_iret,m%temp_tau,m%temp_sigma_n)
        if (m%temp_iret == 0) then
          write(67,*) m%temp_tau
        else 
          write(6,*) '!!WARNING!! : Kernel Singular, set value to 0,(i,j)',i,j
          write(67,*) 0d0
        endif
      end do
    end do
    stop 'Kernel calculation completed and stored in fort.67'
        
  case(68)
    write(6,*) 'Calculate Okada Kernel: Const DX'
    read(iin,*) m%nn,m%nw,m%nx
    read(iin,*) m%temp_lam, m%temp_mu
    allocate(m%x(m%nn),m%y(m%nn),m%z(m%nn),m%dip(m%nn))
    allocate(m%temp_xx(m%nn),m%temp_ww(m%nn))
    do i =1,m%nn
       read(iin,*) m%x(i),m%y(i),m%z(i),m%dip(i),m%temp_xx(i),m%temp_ww(i)
    end do
    do i=1,m%nw
      do j=1,m%nn
        call compute_kernel(m%temp_lam,m%temp_mu,m%x(j),m%y(j),m%z(j),  &
               m%dip(j),m%temp_xx(j),m%temp_ww(j),   &
               m%x(1+(i-1)*m%nx),m%y(1+(i-1)*m%nx),m%z(1+(i-1)*m%nx),   &
               m%dip(1+(i-1)*m%nx),m%temp_iret,m%temp_tau,m%temp_sigma_n)
        if (m%temp_iret == 0) then
!          write(68,'(4e24.16)') m%temp_tau,m%z(j),m%z(1+(i-1)*m%nx),(m%x(j)-m%x(1+(i-1)*m%nx))
          write(68,'(e16.8)') m%temp_tau
        else
          write(6,*) '!!WARNING!! : Kernel Singular, set value to 0,(i,j)',i,j
!          write(68,'(4e24.16)') 0d0,m%z(j),m%z(1+(i-1)*m%nx),(m%x(j)-m%x(1+(i-1)*m%nx))
          write(68,'(e16.8)') 0d0
        endif
      end do
    end do
    stop 'Kernel calculation completed and stored in fort.68'    


  case default
    write(6,*) 'mesh dimension should be 0, 1, 2 or 67/68'

  end select

  if (m%dim==0) m%nn = 1

end subroutine read_mesh

!=============================================================
function mesh_get_size(m) result(n)
  type(mesh_type), intent(inout) :: m
  integer :: n
  n = m%nn
end function mesh_get_size

!=============================================================

subroutine init_mesh(m)

  type(mesh_type), intent(inout) :: m

  write(6,*) 'Initializing mesh ...'

  select case (m%dim)
    case(0); call init_mesh_0D(m)
    case(1); call init_mesh_1D(m)
    case(2); call init_mesh_2D(m)
    case(4); call init_mesh_2D(m)
  end select

end subroutine init_mesh

!--------------------------------------------------------
! Spring-block System
subroutine init_mesh_0D(m)

  type(mesh_type), intent(inout) :: m

  write(6,*) 'Spring-block System' 
  allocate(m%x(m%nn), m%y(m%nn), m%z(m%nn))
  m%dx = m%Lfault
  m%x = 0d0
  m%y = 0d0
  m%z = 0d0

end subroutine init_mesh_0D

!--------------------------------------------------------
! 1D fault, uniform grid
subroutine init_mesh_1D(m)

  type(mesh_type), intent(inout) :: m

  integer :: i

  write(6,*) '1D fault, uniform grid' 
  m%dx = m%Lfault/m%nn
  allocate(m%x(m%nn), m%y(m%nn), m%z(m%nn))
  do i=1,m%nn
    m%x(i) = (i-m%nn*0.5d0-0.5d0)*m%dx
    ! Assuming nn is even (usually a power of 2), 
    ! the center of the two middle elements (i=nn/2 and nn/2+1) 
    ! are located at x=-dx/2 and x=dx/2, respectively 
    m%y(i) = 0d0
    m%z(i) = 0d0
  enddo

end subroutine init_mesh_1D

!--------------------------------------------------------
  ! 2D fault, uniform grid along-strike
  ! Assumptions: 
  !   + the fault trace is parallel to x 
  !   + the lower "left" corner of the fault is at (0,0,Z_CORNER)
  !   + z is negative downwards (-depth)
  ! Storage scheme: faster index runs along-strike (x)
subroutine init_mesh_2D(m)

  use constants, only : PI

  type(mesh_type), intent(inout) :: m

  double precision :: cd, sd, cd0, sd0
  integer :: i, j, j0

  write(6,*) '2D fault, uniform grid along-strike'
  m%dx = m%Lfault/m%nx
  allocate(m%x(m%nn), m%y(m%nn), m%z(m%nn), m%dip(m%nn)) 

  ! set x, y, z, dip of first row
  cd = cos(m%DIP_W(1)/180d0*PI)
  sd = sin(m%DIP_W(1)/180d0*PI)
  do j = 1,m%nx
    m%x(j) = 0d0+(0.5d0+dble(j-1))*m%dx
  end do
  m%y(1:m%nx) = 0d0+0.5d0*m%dw(1)*cd
  m%z(1:m%nx) = m%Z_CORNER+0.5d0*m%dw(1)*sd
  m%dip(1:m%nx) = m%DIP_W(1)

  ! set x, y, z, dip of row 2 to nw
  do i = 2,m%nw
    cd0 = cd
    sd0 = sd
    cd = cos(m%DIP_W(i)/180d0*PI)
    sd = sin(m%DIP_W(i)/180d0*PI)
    j0 = (i-1)*m%nx
    m%x(j0+1:j0+m%nx) = m%x(1:m%nx)
    m%y(j0+1:j0+m%nx) = m%y(j0) + 0.5d0*m%dw(i-1)*cd0 + 0.5d0*m%dw(i)*cd
    m%z(j0+1:j0+m%nx) = m%z(j0) + 0.5d0*m%dw(i-1)*sd0 + 0.5d0*m%dw(i)*sd
    write(66,*) m%z(j0+1:j0+m%nx) 
    m%dip(j0+1:j0+m%nx) = m%DIP_W(i)
  end do


end subroutine init_mesh_2D

end module mesh

