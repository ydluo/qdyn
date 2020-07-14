module mesh

  implicit none
  private

  type mesh_type
    integer :: dim = 0  ! dim = 1, 2 ,3 ~xD
    integer :: nx, nw, nn, nnglob, nwglob ! along-strike, along-dip, total grid number
    double precision :: Lfault, W, Z_CORNER ! fault length, width, lower-left corner z (follow Okada's convention)
    double precision, allocatable :: DIP_W(:) ! along-dip grid size and dip (adjustable), nw count
    double precision, pointer :: time(:) => null() ! time (same for every element)
    ! Local mesh coordinates
    double precision, pointer :: x(:) => null(), y(:) => null(), z(:) => null()
    ! Local dip, and along-strike and along-dip size of grid cells (nx*nw count)
    ! TODO: grid sizes should be constant in order for FFT to work!
    double precision, pointer :: dip(:) => null(), dx(:) => null(), dw(:) => null()
    ! Global mesh coordinates (nx*nwglobal count)
    double precision, pointer :: xglob(:) => null(), yglob(:) => null(), zglob(:) => null()
    double precision, pointer :: dipglob(:) => null(), dwglob(:) => null()
  end type mesh_type

  ! SEISMIC: discretisation of spectral domain for diffusion solver
  type spectral_mesh_type
    double precision, dimension(:), allocatable :: lw, F_inv ! dimensionless wavenumbers, inv Fourier kernel
    double precision :: Dlogl=0.8, lw_max=10.0 ! logarithmic grid spacing, max dimensionless wavenumber
    integer :: Nl=20 ! number of mesh elements
  end type spectral_mesh_type

  integer, allocatable, save :: nnLocal_perproc(:),nnoffset_glob_perproc(:)

  public :: mesh_type, read_mesh_parameters, read_mesh_nodes, init_mesh, mesh_get_size
  public :: nnLocal_perproc, nnoffset_glob_perproc
  public :: spectral_mesh_type

contains

!=============================================================
subroutine read_mesh_parameters(iin,m)

  type(mesh_type), intent(inout) :: m
  integer, intent(in) :: iin

  integer :: i

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

  case default
    write(6,*) 'mesh dimension should be 0, 1 or 2'

  end select

  if (m%dim==0) m%nn = 1

end subroutine read_mesh_parameters

!=============================================================
subroutine read_mesh_nodes(iin,m)

  type(mesh_type), intent(inout) :: m
  integer, intent(in) :: iin

  integer :: i

  allocate(m%x(m%nn), m%y(m%nn), m%z(m%nn), m%dip(m%nn))
  do i=1,m%nn
    read(iin,*) m%x(i),m%y(i),m%z(i),m%dip(i)
  enddo

end subroutine read_mesh_nodes

!=============================================================
integer function mesh_get_size(m) result(n)
  type(mesh_type), intent(in) :: m
  n = m%nn
end function mesh_get_size

!=============================================================

subroutine init_mesh(m)

  type(mesh_type), intent(inout) :: m

  write(6,*) 'Initializing mesh ...'

  ! Allocate time array (same for all m%dim)
  allocate(m%time(m%nn))
  m%time = 0d0
  m%nnglob = m%nn

  select case (m%dim)
    case(0); call init_mesh_0D(m)
    case(1); call init_mesh_1D(m)
    case(2); call init_mesh_2D(m)
  end select

end subroutine init_mesh

!--------------------------------------------------------
! Spring-block System
subroutine init_mesh_0D(m)

  type(mesh_type), intent(inout) :: m

  write(6,*) 'Spring-block System'
  allocate(m%x(m%nn), m%y(m%nn), m%z(m%nn))
  allocate(m%dx(1))
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
  allocate(m%dx(1))
  m%dx = m%Lfault/m%nn
  allocate(m%x(m%nn), m%y(m%nn), m%z(m%nn))
  do i=1,m%nn
    m%x(i) = (i-m%nn*0.5d0-0.5d0)*m%dx(1)
    ! Assuming nn is even (usually a power of 2),
    ! the center of the two middle elements (i=nn/2 and nn/2+1)
    ! are located at x=-dx/2 and x=dx/2, respectively
    m%y(i) = 0d0
    m%z(i) = 0d0
  enddo

  m%xglob => m%x
  m%yglob => m%y
  m%zglob => m%z

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
  use my_mpi, only: is_MPI_parallel, my_mpi_NPROCS, gather_alli, gather_allvdouble, my_mpi_tag

  type(mesh_type), intent(inout) :: m

  double precision :: cd, sd, cd0, sd0
  integer :: i, j, j0
  integer :: iproc, nwLocal, nnLocal, nwGlobal, nnGlobal, NPROCS
  integer, allocatable :: nwLocal_perproc(:), nwoffset_glob_perproc(:)

  write(6,*) '2D fault, uniform grid along-strike'

  allocate(m%dx(1))

if (.not.is_MPI_parallel()) then

  m%dx = m%Lfault/m%nx
  allocate(m%x(m%nn), m%y(m%nn), m%z(m%nn), m%dip(m%nn))

  ! set x, y, z, dip of first row
  cd = cos(m%DIP_W(1)/180d0*PI)
  sd = sin(m%DIP_W(1)/180d0*PI)
  do j = 1,m%nx
    m%x(j) = 0d0+(0.5d0+dble(j-1))*m%dx(1)
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
!    write(66,*) m%z(j0+1:j0+m%nx) !JPA Who is using this output? Shall we remove it?
    m%dip(j0+1:j0+m%nx) = m%DIP_W(i)
  end do

    m%xglob => m%x
    m%yglob => m%y
    m%zglob => m%z
    m%dipglob => m%dip
    m%dwglob => m%dw
    m%nwglob = m%nw

else
! If MPI parallel, the mesh for each processor has been already read
! from qdynxxx.in and initialized in input.f90
! Here we gather the whole fault to compute the kernel
    m%dx = m%Lfault/m%nx
    NPROCS = my_mpi_NPROCS()
!PG, for debugging:
    nwLocal = m%nw
    write(6,*) 'iproc,nwLocal:',my_mpi_tag(),nwLocal
    allocate(nwLocal_perproc(0:NPROCS-1))
    call gather_alli(nwLocal,nwLocal_perproc)
    allocate(nwoffset_glob_perproc(0:NPROCS-1))
    do iproc=0,NPROCS-1
      nwoffset_glob_perproc(iproc)=sum(nwLocal_perproc(0:iproc))-nwLocal_perproc(iproc)
    enddo
    nwGlobal=sum(nwLocal_perproc)
    m%nwglob = nwGlobal
    write(6,*) 'iproc,nwGlobal:',my_mpi_tag(),nwGlobal

    nnLocal= m%nx*nwLocal
    allocate(nnLocal_perproc(0:NPROCS-1))
    nnLocal_perproc = nwLocal_perproc * m%nx
    allocate(nnoffset_glob_perproc(0:NPROCS-1))
    nnoffset_glob_perproc = nwoffset_glob_perproc * m%nx
    nnGlobal=nwGlobal*m%nx
    m%nnglob = nnGlobal

   ! global mesh for computing the kernel and for outputs
    allocate(m%xglob(nnGlobal),m%yglob(nnGlobal),&
             m%zglob(nnGlobal),m%dwglob(nwGlobal),m%dipglob(nnGlobal))
    call gather_allvdouble(m%x, nnLocal,m%xglob,nnLocal_perproc,nnoffset_glob_perproc,nnGlobal)
    call gather_allvdouble(m%y, nnLocal,m%yglob,nnLocal_perproc,nnoffset_glob_perproc,nnGlobal)
    call gather_allvdouble(m%z, nnLocal,m%zglob,nnLocal_perproc,nnoffset_glob_perproc,nnGlobal)
    call gather_allvdouble(m%dip, nnLocal,m%dipglob,nnLocal_perproc,nnoffset_glob_perproc,nnGlobal)
    call gather_allvdouble(m%dw, nwLocal,m%dwglob,nwLocal_perproc,nwoffset_glob_perproc,nwGlobal)

endif

end subroutine init_mesh_2D

end module mesh
