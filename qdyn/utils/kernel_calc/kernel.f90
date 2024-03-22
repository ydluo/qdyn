! Computes Okada shear stress kernel 
!
! To compile: 
! 1. In qdyn/src/, run make or 
!    ifort -O3 -c okada.f90
! 2. In utils/kernel_calc/, run: 
!    ifort -O3 -o kernel.exe  -I../../src/ kernel.f90 ../../src/okada.o
!
! Input file kernel.in. Format line by line:
! 1. fault type (integer)
!     1 = right-lateral strike-slip 
!    -1 = left-lateral strike-slip 
!     2 = thrust dip-slip
!    -2 = normal dip-slip
! 2. lambda, mu (two reals) = elastic moduli
! 3. mesh type (integer)
!     1 = general mesh
!     2 = mesh has constant strike and dx, 
! 4. number of elements (integer) if mesh type = 1
!    number of elements along strike, along dip (two integers) if mesh type = 2
! 5-END. Information for each element (6 reals, space delimited)
!    x, y, z, dip, dx (size along strike), dw (size along dip)
!    If mesh type = 1, no particular order is assumed
!    If mesh type = 2, elements must be ordered as (along-strike,along-dip)
!
! Output binary files kernel.s.out and kernel.n.out
! contain the kernel matrices for shear and normal stress, respectively.
! If mesh type = 1, the whole matrix K(i,j) is computed
! If mesh type = 2, only K(x,z,0,z') is computed

program kernel

  use okada, only : compute_kernel
  use logger, only : log_msg

  integer, parameter :: iin =15, iout_s=16, iout_n=17

  double precision, dimension(:), allocatable :: x,y,z,dip,xx,ww 
  double precision :: mu, lam, tau, sigma
  integer :: FAULT_TYPE
  integer :: iret, i,j,di, mode, nn, nw, nx
  character(100) :: msg

! read inputs
  open(unit=iin,FILE= 'kernel.in') 

  read(iin,*) FAULT_TYPE
  read(iin,*) lam, mu
  read(iin,*) mode

  select case (mode)
  case(1) 
    call log_msg("Calculate Okada kernel for general mesh")
    read(iin,*) nn
    di = 1
  case(2)
    call log_msg("Calculate Okada kernel for constant strike and DX")
    read(iin,*) nx,nw
    nn = nw*nx
    di = nx
  case default
    call log_msg("Input mode should be 1 (general mesh) or 2 (fixed strike and dx)")
  end select

  allocate(x(nn),y(nn),z(nn),dip(nn),xx(nn),ww(nn))
  do i =1,nn
    read(iin,*) x(i),y(i),z(i),dip(i),xx(i),ww(i)
  end do

  close(iin)

! compute and output kernels
  open(unit=iout_s,FILE= 'kernel.s.out',access='stream',status='replace') 
  open(unit=iout_n,FILE= 'kernel.n.out',access='stream',status='replace') 
        
  do i=1,nn,di  ! receivers 
  do j=1,nn     ! sources 
    call compute_kernel(lam,mu,x(j),y(j),z(j), dip(j),xx(j),ww(j), &
                        x(i),y(i),z(i), dip(i), iret,tau,sigma,FAULT_TYPE)
    if (iret /= 0) then
      write(msg, *) "WARNING : Kernel singular, set value to 0, (i,j)=", i, j
      call log_msg(msg)
      tau = 0d0
      sigma = 0d0
    endif
    write(iout_s) tau
    write(iout_n) sigma
  end do
  end do

  close(iout_s)
  close(iout_n)

  call log_msg("Kernel calculation completed and stored in kernel.*.out")
  
end program kernel
