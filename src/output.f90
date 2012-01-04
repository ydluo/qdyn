module output

! OUTPUT: This module manages outputs
!

!  use some_module_1

  implicit none
  public ot_type, ox_type
  private

 ! timeseries outputs: at every time step, but only macroscopic quantities
  type ot_type
    private
    double precision, pointer ::
    integer, pointer ::
    logical, pointer ::
    double precision :: lcold,lcnew,llocnew,llocold
    integer :: unit,ic
    logical ::
  end type ot_type

 ! snapshot outputs: at every fault point, but only at few selected times
  type ox_type
    private
    double precision, pointer ::
    integer, pointer ::
    logical, pointer ::
    double precision ::
    integer :: count,unit,nx
    logical ::
  end type ox_type

  public :: ot_type, ox_type, ot_init, ox_init

contains

!=====================================================================
! write ot file header
subroutine ot_init(ot,nn)

  type (ot_type), intent(inout) :: ot
 
  ot%ic = nn/2
  if (ot%ic == 0) ot%ic = 1

  ot%lcnew = dble(nn)
  ot%llocnew = dble(nn)

  ot%unit = 18
 
  write(ot%unit,'(a)')'# macroscopic values:'
  write(ot%unit,'(a)')'# 1=t,2=loc_size,3=crack_size,4=potcy,5=pot_rate'
  write(ot%unit,'(a)')'# values at center:'
  write(ot%unit,'(a)')'# 6=V, 7=theta, 8=V*theta/dc, 9=tau, 10=slip'
  write(ot%unit,'(a)')'# values at max(V) location:'
  write(ot%unit,'(a)')'# 11=x, 12=V, 13=theta, 14=omeg, 15=tau, 16=slip'

end subroutine ot_init


!=====================================================================
! write ox file header
subroutine ox_init(ox,nn)
 
  type (ox_type), intent(inout) :: ox
  integer, intent(in) :: nn

  ox%unit = 19

  ox%count=0
  do i=1,nn,ox%nx
    ox%count = ox%count+1
  enddo
  write(ox%unit,'(a,2i5)')'# nx= ',ox%count

subroutine ox_init

!=====================================================================
! Export timeseries
subroutine ot_write(ot)
 
  type (ot_type), intent(inout) :: ot

  write(ot%unit,'(e24.16,15e14.6)') pb%time, pb%output%llocnew*pb%mesh%dx,  &
    pb%output%lcnew*pb%mesh%dx, pot, pot_rate,    &
    pb%v(ot%ic), pb%theta(ot%ic), pb%v(ot%ic)*pb%theta(ot%ic)/pb%dc(ot%ic), &
    pb%tau(ot%ic), pb%slip(ot%ic),    &
    xvmax, vmax, pb%theta(ivmax), vmax*pb%theta(ivmax)/pb%dc(ivmax),    &
    pb%tau(ivmax), pb%slip(ivmax)

end subroutine ot_write

!=====================================================================
! Export snapshots
subroutine ox_write(ox,it,ivmax,pb)
 
  type (ox_type), intent(inout) :: ox

  integer :: i

  write(ox%unit,'(2a,2i5,e14.6)')'# x v theta V./V dtau tau_dot slip ',it,ivmax,pb%time
  do i=1,pb%mesh%nn,ox%nx
    write(19,'(8e15.7)') pb%mesh%x(i),pb%time,pb%v(i), &
      pb%theta(i),pb%dvdt(i)/pb%v(i),pb%tau(i),pb%dtau_dt(i),pb%slip(i)
  enddo

end subroutine ox_write

!=====================================================================
! distance between largest peak on the left half
! and largest peak on the right half
function crack_size(s,n)
       
  integer :: n
  double precision ::  s(n)
  double precision ::  smin,smax,s1,s2,xL,xR
  integer :: i,iL,iR,imin
  double precision :: crack_size

  if (n == 0) then
    crack_size = 0d0
    return
  endif

  ! assuming minimum stressing rate at nucleation point
  !      smin=s(1)
  !      do i=1,n
  !        if (s(i).le.smin) then
  !          imin=i
  !          smin=s(i)
  !        endif
  !      enddo

  ! WARNING: assumes one peak on x<0 and one peak on x>0
  imin=n/2

  iL = maxloc(s(1:imin))
  smax=s(iL)
  xL = dble(iL);
  if (iL > 1) then
    s1 = 0.5d0*( s(iL+1)-s(iL-1) )
    s2 = s(iL-1)-2d0*s(iL)+s(iL+1)
    if (s2 /= 0d0) xL = xL-s1/s2
  endif
      
  iR = maxloc(s(imin:n)
  smax=s(iR)
  xR = dble(iR);
  if (iR < n) then
    s1 = 0.5d0*( s(iR+1)-s(iR-1) )
    s2 = s(iR-1)-2d0*s(iR)+s(iR+1)
    if (s2 /= 0d0) xR = xR-s1/s2
  endif

  crack_size = xR-xL

end function crack_size


end module mod_name
