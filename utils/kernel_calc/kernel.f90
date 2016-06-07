! Computes Okada shear stress kernel 
! To compile: ifort -o kernel.exe  -I../../src/ kernel.f90 ../../src/okada.o

program kernel

  use okada, only : compute_kernel

  integer, parameter :: iin =15, iout=16

  double precision, dimension(:), allocatable :: x,y,z,dip,xx,ww 
  double precision :: mu, lam, tau, sigma
  integer :: iret, i,j,ii, mode, nn, nw, nx

  open(unit=iin,FILE= 'kernel.in') 
  open(unit=iout,FILE= 'kernel.out') 

  read(iin,*) mode

  select case (mode)

  ! mode = 1 : general mesh, compute all (i,j) pairs
  case(1) 
    write(6,*) 'Calculate Okada kernel for general mesh'
    read(iin,*) nn
    read(iin,*) lam, mu
    allocate(x(nn),y(nn),z(nn),dip(nn),xx(nn),ww(nn))
    do i =1,nn
      read(iin,*) x(i),y(i),z(i),dip(i),xx(i),ww(i)
    end do
    do i=1,nn
      do j=1,nn
        call compute_kernel(lam,mu,x(i),y(i),z(i),dip(i),xx(i),ww(i), &
               x(j),y(j),z(j),dip(j),iret,tau,sigma)
        if (iret /= 0) then
          write(6,*) 'WARNING : Kernel singular, set value to 0, (i,j)=',i,j
          tau = 0d0
        endif
        write(iout,*) tau
      end do
    end do
        
  ! mode = 2 : mesh has constant strike and dx, compute only non-redundant (i,j) pairs
  case(2)
    write(6,*) 'Calculate Okada kernel for constant strike and DX'
    read(iin,*) nn,nw,nx
    read(iin,*) lam, mu
    allocate(x(nn),y(nn),z(nn),dip(nn),xx(nn),ww(nn))
    do i =1,nn
      read(iin,*) x(i),y(i),z(i),dip(i),xx(i),ww(i)
    end do
    do i=1,nw
      ii = 1+(i-1)*nx
      do j=1,nn
        call compute_kernel(lam,mu,x(j),y(j),z(j), dip(j),xx(j),ww(j), &
               x(ii),y(ii),z(ii), dip(ii), iret,tau,sigma)
        if (iret /= 0) then
          write(6,*) 'WARNING : Kernel singular, set value to 0, (i,j)=',i,j
          tau = 0d0
        endif
        write(iout,'(e16.8)') tau
      end do
    end do

  case default
    write(6,*) 'mode should be 1 or 2'

  end select

  write(6,*) 'Kernel calculation completed and stored in kernel.out'    
  close(iin)
  close(iout)
  
end program kernel
