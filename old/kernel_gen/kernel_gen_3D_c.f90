! Generate 3D kernel for qdyn

program main

  use calc_kernel

  implicit none
  type mesh_type
    integer :: w_kind     ! w_kind = 0 uniform grid along-dip 
    integer :: nx, nw     ! along-strike, along-dip grid number
    double precision :: L, W, Z_CORNER ! fault length, width, lower-left corner z (follow Okada's convension)
    double precision :: dx   !along-strike grid size(constant)
    double precision, allocatable :: dw(:), DIP_W(:) !along-dip grid size and dip (adjustable)
    double precision, allocatable :: x(:), y(:), z(:), dip(:) !coordinates and dip of every grid (nx*nw count)
  end type mesh_type
  
!kernel(i,j): response at i of source at j
  double precision, allocatable :: kernel(:,:), temp(:)  ! 2d compact kernel (nx*nw*nw) count

  double precision :: tau = 0d0   ! return kernel from compute_kernel
  double precision :: PI = 3.1415926535897932384626d0
  double precision :: dz0, dip0, cd, sd, cd0, sd0
  double precision :: MU, LAM
  integer :: i, j, IRET, l, mo
  
  type(mesh_type) :: mesh
   
  
  dip0 = 30d0
  MU = 30d9
  LAM = 30d9
  

  mesh%w_kind = 0
  mesh%L = 100d3
  mesh%W = 1d3
  mesh%nx = 1
  mesh%nw = 100
  mesh%Z_CORNER = -50d3
  
  mesh%dx = mesh%L/mesh%nx
  

  allocate (mesh%dw(mesh%nw), mesh%DIP_W(mesh%nw),   &
    mesh%x(mesh%nw*mesh%nx),mesh%y(mesh%nw*mesh%nx),  &
    mesh%z(mesh%nw*mesh%nx), mesh%dip(mesh%nw*mesh%nx),  &
    kernel(mesh%nw*mesh%nx,mesh%nw*mesh%nx), temp(mesh%nw*mesh%nx))
  
!------ give value to dw, DIP_W -------------------
!!!!!!!!! should read in if w_kind != 0 
  if (mesh%w_kind == 0) then
    dz0 = mesh%W/mesh%nw
    do i = 1,mesh%nw
      mesh%dw(i) = dz0
      mesh%DIP_W(i) = dip0
    end do
  end if
!------ give value to dw, DIP_W -------------------
!!!!!!!!! should read in if w_kind != 0 

!------ give value to x, y, z , dip of first row -------------------
  cd = dcos(mesh%DIP_W(1)/180d0*PI)
  sd = dsin(mesh%DIP_W(1)/180d0*PI)
  do j = 1,mesh%nx
    mesh%x(j) = 0d0+(0.5d0+dble(j-1))*mesh%dx
    mesh%y(j) = 0d0+0.5d0*mesh%dw(1)*cd
    mesh%z(j) = mesh%Z_CORNER+0.5d0*mesh%dw(1)*sd
    mesh%dip(j) = mesh%DIP_W(1)
  end do
!------ give value to x, y, z , dip of first row -------------------


!------ give value to x, y, z , dip of row 2 to nw------------------- 
  do i = 2,mesh%nw
    cd0 = dcos(mesh%DIP_W(i-1)/180d0*PI)
    sd0 = dsin(mesh%DIP_W(i-1)/180d0*PI)
    cd = dcos(mesh%DIP_W(i)/180d0*PI)
    sd = dsin(mesh%DIP_W(i)/180d0*PI)
    do j = 1,mesh%nx
      mesh%x((i-1)*mesh%nx+j) = 0d0+(0.5d0+dble(j-1))*mesh%dx
      mesh%y((i-1)*mesh%nx+j) = mesh%y((i-2)*mesh%nx+j)+0.5d0*mesh%dw(i-1)*cd0+0.5d0*mesh%dw(i)*cd
      mesh%z((i-1)*mesh%nx+j) = mesh%z((i-2)*mesh%nx+j)+0.5d0*mesh%dw(i-1)*sd0+0.5d0*mesh%dw(i)*sd
      mesh%dip((i-1)*mesh%nx+j) = mesh%DIP_W(i)
    end do
  end do
!------ give value to x, y, z , dip of row 2 to nw------------------- 

!------ calculate kernel -----------------
!kernel(i,j): response at i of source at j
!because dx = constant, only need to calculate i at first column
  do i = 1,mesh%nw
    do j = 1,mesh%nw*mesh%nx
      call compute_kernel(LAM,MU,mesh%x(j),mesh%y(j),mesh%z(j),  &
                      mesh%dip(j),mesh%dx,mesh%dw((j-1)/mesh%nx+1),   &
                      mesh%x(1+(i-1)*mesh%nx),mesh%y(1+(i-1)*mesh%nx),   &
                      mesh%z(1+(i-1)*mesh%nx),mesh%dip(1+(i-1)*mesh%nx),IRET,tau)
      if (IRET == 0) then
        kernel(i,j) = tau
      else
        write(6,*) 'Kernel Singular, set value to 0'
        kernel(i,j) = 0
      end if
    end do
  end do
 
  write(6,*) 'x'
  write(6,*) mesh%x 
  write(6,*) 'y'
  write(6,*) mesh%y 
  write(6,*) 'z'
  write(6,*) mesh%z 
  write(6,*) 'dip'
  write(6,*) mesh%dip 
  write(6,*) 'kernel(i,j): response at i of source at j: store first column only'
  do i = 1,mesh%nw
    temp = kernel(i,:)
    write(6,*) 'row=',i,' column=',1  
    write(6,*) temp
    do l = 2, mesh%nx
      write(6,*) l
      do j =1, mesh%nw*mesh%nx
        mo = mod(j,mesh%nx)
        if (mo >= l .or. mo == 0) then
          temp(j) = kernel(i,j-l+1)
        else
          temp(j) = kernel(i,j+1+l-2*mo)
        end if
      end do
      write(6,*) 'row=',i,' column=',l  
      write(6,*) temp
    end do     
  end do
   write(6,*) mesh%nx

end program main
