module output

! OUTPUT: This module manages outputs
!

!  use some_module_1

  implicit none
  private :: ixout
  public :: ot_init, ox_init, ot_write, ox_write

contains

!=====================================================================
!output initilized field to screen
subroutine screen_init(pb)

  use problem_class
  type (problem_type), intent(inout) :: pb
  
  write(6,*) '**Field Initialized.**'
  write(6,*) 'Values at selected point of the fault:'
    
    if (pb%kernel%k2f%finite == 1 .or. pb%mesh%nn == 1) then
      write(6,*) 'K/Kc = ',(PI*pb%smu/pb%mesh%Lfault)/   &
        (pb%sigma(pb%ot%ic)*(pb%b(pb%ot%ic)-pb%a(pb%ot%ic))/pb%dc(pb%ot%ic))
      write(6,*) 'K/Kb = ',(PI*pb%smu/pb%mesh%Lfault)/   &
        (pb%sigma(pb%ot%ic)*pb%b(pb%ot%ic)/pb%dc(pb%ot%ic))
    else
      write(6,*) 'K/Kc = ',(PI*pb%smu/pb%mesh%W)/   &
        (pb%sigma(pb%ot%ic)*(pb%b(pb%ot%ic)-pb%a(pb%ot%ic))/pb%dc(pb%ot%ic))
      write(6,*) 'K/Kb = ',(PI*pb%smu/pb%mesh%W)/   &
        (pb%sigma(pb%ot%ic)*pb%b(pb%ot%ic)/pb%dc(pb%ot%ic))
    endif
    write(6,*)
    write(6,*) '    it,  dt (secs), time (yrs), vmax (m/s)'

end subroutine screen_init



!=====================================================================
!output one step to screen
subroutine screen_write(pb,it,dt_did)
  
  use constant, only : YEAR
  use problem_class
  type (problem_type), intent(inout) :: pb
  integer :: it
  double precision :: dt_did

  write(6,'(i7,x,3(e11.3,x),i5)') it, dt_did, pb%time/YEAR, pb%v(pb%it%ivmax)

end subroutine screen_write


!=====================================================================
! write time of every step
subroutine time_write(pb)
 
  use problem_class
  type (problem_type), intent(inout) :: pb

  write(121,*) pb%time, pb%tmax

end subroutine time_write


!=====================================================================
! write ot file header
subroutine ot_init(pb)

  use problem_class
  type (problem_type), intent(inout) :: pb

  pb%ot%lcnew = dble(pb%mesh%nn)
  pb%ot%llocnew = dble(pb%mesh%nn)

  pb%ot%unit = 18
 
  write(pb%ot%unit,'(a)')'# macroscopic values:'
  write(pb%ot%unit,'(a)')'# 1=t,2=loc_size,3=crack_size,4=potcy,5=pot_rate'
  write(pb%ot%unit,'(a)')'# values at center:'
  write(pb%ot%unit,'(a)')'# 6=V, 7=theta, 8=V*theta/dc, 9=tau, 10=slip'
  write(pb%ot%unit,'(a)')'# values at max(V) location:'
  write(pb%ot%unit,'(a)')'# 11=x, 12=V, 13=theta, 14=omeg, 15=tau, 16=slip'

end subroutine ot_init



!=====================================================================
! write ox file header
subroutine ox_init(pb)
 
  use problem_class
  type (problem_type), intent(inout) :: pb

  pb%ox%unit = 19

  pb%ox%count=0
  do i=1,pb%mesh%nn,pb%ox%nxout
    pb%ox%count = pb%ox%count+1
  enddo
  write(pb%ox%unit,'(a,2i5)')'# nx= ',pb%ox%count

subroutine ox_init



!=====================================================================
! Export timeseries
subroutine ot_write(pb)
 
  use problem_class
  type (problem_type), intent(inout) :: pb

  write(ot%unit,'(e24.16,15e14.6)') pb%time, pb%output%llocnew*pb%mesh%dx,  &
    pb%output%lcnew*pb%mesh%dx, pot, pot_rate,    &
    pb%v(pb%ot%ic), pb%theta(pb%ot%ic),  &
    pb%v(pb%ot%ic)*pb%theta(pb%ot%ic)/pb%dc(pb%ot%ic), &
    pb%tau(pb%ot%ic), pb%slip(pb%ot%ic),    &
    pb%mesh%x(pb%ot%ivmax), pb%v(pb%ot%ivmax), pb%theta(pb%ot%ivmax),   &
    pb%v(pb%ot%ivmax)*pb%theta(pb%ot%ivmax)/pb%dc(pb%ot%ivmax),    &
    pb%tau(pb%ot%ivmax), pb%slip(pb%ot%ivmax)

end subroutine ot_write




!=====================================================================
! Export snapshots
subroutine ox_write(it,pb)
 
  use problem_class
  type (problem_type), intent(inout) :: pb

  integer :: it,ixout

  write(ox%unit,'(2a,2i5,e14.6)')'# x v theta V./V dtau tau_dot slip ',it,pb%ot%ivmax,pb%time
  do ixout=1,pb%mesh%nn,pb%ox%nxout
    write(19,'(8e15.7)') pb%mesh%x(ixout),pb%time,pb%v(ixout),   &
      pb%theta(ixout),pb%dvdt(ixout)/pb%v(ixout),pb%tau(ixout),   &
      pb%dtau_dt(ixout),pb%slip(ixout)
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
