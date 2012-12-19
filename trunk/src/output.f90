module output

! OUTPUT: This module manages outputs
!

!  use some_module_1

  implicit none
  private  
  public :: screen_init, ot_init, ox_init, &
            screen_write, ot_write, ox_write,  &
            time_write, crack_size

contains

!=====================================================================
!output initilized field to screen
subroutine screen_init(pb)

  use problem_class
  use constants, only : PI
  type (problem_type), intent(inout) :: pb
  
  write(6,*) '**Field Initialized.**'
  write(6,*) 'Values at selected point of the fault:'
    if (pb%mesh%dim == 0 .or.pb%mesh%dim == 1) then   
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
    end if

! YD:  should we out put 2D mesh K?

    write(6,*)
    write(6,*) '    it,  dt (secs), time (yrs), vmax (m/s)'

end subroutine screen_init



!=====================================================================
!output one step to screen
subroutine screen_write(pb)
  
  use constants, only : YEAR
  use problem_class
  type (problem_type), intent(inout) :: pb

  write(6,'(i7,x,3(e11.3,x),i5)') pb%it, pb%dt_did, pb%time/YEAR, pb%v(pb%ot%ivmax)

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
  integer :: i
  pb%ot%lcnew = dble(pb%mesh%nn)
  pb%ot%llocnew = dble(pb%mesh%nn)

  pb%ot%unit = 18 
  write(pb%ot%unit,'(a)')'# macroscopic values:'
  write(pb%ot%unit,'(a)')'# 1=t,2=loc_size,3=crack_size,4=potcy,5=pot_rate'
  write(pb%ot%unit,'(a)')'# values at center:'
  write(pb%ot%unit,'(a)')'# 6=V, 7=theta, 8=V*theta/dc, 9=tau, 10=slip'
  write(pb%ot%unit,'(a)')'# values at max(V) location:'
  write(pb%ot%unit,'(a)')'# 11=x, 12=V, 13=theta, 14=omeg, 15=tau, 16=slip'

  pb%ot%unit = 22
  write(pb%ot%unit,'(a)')'# Seismicity record:' 
  write(pb%ot%unit,'(a)')'# 1=loc, 2=t, 3=v'

  pb%ot%unit = 10000
  do i=1,pb%mesh%nn
    if (pb%iot(i) == 1) then
      pb%ot%unit = pb%ot%unit+1
      write(pb%ot%unit,'(a,i10)')'# nx= ', i
      write(pb%ot%unit,'(a)')'# 1=t, 2=V, 3=theta, 4=tau, 5=slip'
    endif
  enddo


end subroutine ot_init



!=====================================================================
! write ox file header
subroutine ox_init(pb)
 
  use problem_class
  type (problem_type), intent(inout) :: pb
  integer :: i
  if (pb%ox%i_ox_seq == 0) then
    pb%ox%unit = 19
  else
    pb%ox%unit = 1000
  endif

  pb%ox%count=0
  do i=1,pb%mesh%nn,pb%ox%nxout
    pb%ox%count = pb%ox%count+1
  enddo
  write(pb%ox%unit,'(a,i10)')'# nx= ',pb%ox%count

end subroutine ox_init



!=====================================================================
! Export timeseries
subroutine ot_write(pb)
 
  use problem_class
  type (problem_type), intent(inout) :: pb
  integer :: i

  pb%ot%unit = 18
  write(pb%ot%unit,'(e24.16,15e14.6)') pb%time, pb%ot%llocnew*pb%mesh%dx,  &
    pb%ot%lcnew*pb%mesh%dx, pb%pot, pb%pot_rate,    &
    pb%v(pb%ot%ic), pb%theta(pb%ot%ic),  &
    pb%v(pb%ot%ic)*pb%theta(pb%ot%ic)/pb%dc(pb%ot%ic), &
    pb%tau(pb%ot%ic), pb%slip(pb%ot%ic),    &
    pb%mesh%x(pb%ot%ivmax), pb%v(pb%ot%ivmax), pb%theta(pb%ot%ivmax),   &
    pb%v(pb%ot%ivmax)*pb%theta(pb%ot%ivmax)/pb%dc(pb%ot%ivmax),    &
    pb%tau(pb%ot%ivmax), pb%slip(pb%ot%ivmax)


  pb%ot%unit = 22
  do i=1,pb%mesh%nn
    if ((pb%iasp(i) == 1) .and. (pb%v(i) >= pb%v_th) .and.      &
        (pb%v(i) < pb%v_pre(i)) .and. (pb%v_pre(i) >= pb%v_pre2(i))) then
      write(pb%ot%unit,'(i10,2e24.16)') i, pb%time, pb%v(i)
    endif
  enddo
  pb%v_pre2=pb%v_pre
  pb%v_pre=pb%v



  pb%ot%unit = 10000

  do i=1,pb%mesh%nn
    if (pb%iot(i) == 1) then
      pb%ot%unit = pb%ot%unit+1
      write(pb%ot%unit,'(e24.16,4e14.6)') pb%time, pb%v(i),      &
      pb%theta(i), pb%tau(i), pb%slip(i)
    endif
  enddo




end subroutine ot_write




!=====================================================================
! Export snapshots
subroutine ox_write(pb)
 
  use problem_class
  type (problem_type), intent(inout) :: pb

  integer :: ixout

  if (pb%ox%i_ox_seq == 0) then
    write(pb%ox%unit,'(2a,2i8,e14.6)')'# x v theta',' V./V dtau tau_dot slip ',pb%it,pb%ot%ivmax,pb%time
  ! JPA: this output should also contain y and z
    do ixout=1,pb%mesh%nn,pb%ox%nxout
      write(pb%ox%unit,'(e15.7,e24.16,6e15.7)') pb%mesh%x(ixout),pb%time,pb%v(ixout),   &
        pb%theta(ixout),pb%dv_dt(ixout)/pb%v(ixout),pb%tau(ixout),   &
        pb%dtau_dt(ixout),pb%slip(ixout)
    enddo
  else
    pb%ox%unit = pb%ox%unit + 1
    write(pb%ox%unit,'(3i10,e24.14)') pb%it,pb%ot%ivmax,pb%ox%count,pb%time
    write(pb%ox%unit,'(2a)') '#  x  y  z  t  v  theta','  V./V  dtau  tau_dot  slip '
    do ixout=1,pb%mesh%nn,pb%ox%nxout
      write(pb%ox%unit,'(3e15.7,e24.14,6e15.7)')       &
        pb%mesh%x(ixout),pb%mesh%y(ixout),pb%mesh%z(ixout),pb%time,     &
        pb%v(ixout),pb%theta(ixout),pb%dv_dt(ixout)/pb%v(ixout),pb%tau(ixout),   &
        pb%dtau_dt(ixout),pb%slip(ixout)
    enddo
    close(pb%ox%unit)
  endif
  


  if (pb%DYN_FLAG == 1)

    if (pb%ox%dyn_stat == 0) .and. (pb%v(pb%ot%ivmax) >= pb%DYN_th_on )
      pb%ox%dyn_stat = 1
      OPEN (UNIT = 100, FILE='DYN_PRE.txt', STATUS='REPLACE')
      write(100,'(3i10,e24.14)') pb%it,pb%ot%ivmax,pb%ox%count,pb%time
      write(100,'(2a)') '#  x  y  z  t  v  theta','  V./V  dtau  tau_dot  slip '
      do ixout=1,pb%mesh%nn,pb%ox%nxout
        write(pb%ox%unit,'(3e15.7,e24.14,6e15.7)')       &
          pb%mesh%x(ixout),pb%mesh%y(ixout),pb%mesh%z(ixout),pb%time,     &
          pb%v(ixout),pb%theta(ixout),pb%dv_dt(ixout)/pb%v(ixout),pb%tau(ixout),   &
          pb%dtau_dt(ixout),pb%slip(ixout)
      enddo
      pb%pot_pre = pb%pot
      CLOSE(100)
    endif

    if (pb%ox%dyn_stat == 1) .and. (pb%v(pb%ot%ivmax) <= pb%DYN_th_off )
      pb%ox%dyn_stat = 0
      OPEN (UNIT = 101, FILE='DYN_POST.txt', STATUS='REPLACE')
      write(101,'(3i10,e24.14)') pb%it,pb%ot%ivmax,pb%ox%count,pb%time
      write(101,'(2a)') '#  x  y  z  t  v  theta','  V./V  dtau  tau_dot  slip '
      do ixout=1,pb%mesh%nn,pb%ox%nxout
        write(pb%ox%unit,'(3e15.7,e24.14,6e15.7)')       &
          pb%mesh%x(ixout),pb%mesh%y(ixout),pb%mesh%z(ixout),pb%time,     &
          pb%v(ixout),pb%theta(ixout),pb%dv_dt(ixout)/pb%v(ixout),pb%tau(ixout),   &
          pb%dtau_dt(ixout),pb%slip(ixout)
      enddo
      CLOSE(101)
      if ((pb%pot-pb%pot_pre)*pb%smu >= pb%DYN_M)
        pb%ox%dyn_count = pb%ox%dyn_count + 1
        if (pb%ox%dyn_count > pb%DYN_SKIP)
          pb%itstop = pb%it
          write(222,'(3i10,e24.14)') pb%it,pb%ot%ivmax,pb%ox%count,pb%time
        endif
      endif
    endif 
   
  endif  
  

end subroutine ox_write

!=====================================================================
! distance between largest peak on the left half
! and largest peak on the right half
function crack_size(s,n)
       
  integer :: n
  double precision ::  s(n)
  double precision ::  stemp,smax,s1,s2,xL,xR
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

  iL=1
  stemp=0d0;
  do i=1,imin
    if (s(i) > stemp) then
      iL = i
      stemp = s(i)
    end if
  end do

  iR=imin
  stemp=0d0
  do i=imin,n
    if (s(i) >= stemp) then
      iR = i
      stemp = s(i)
    end if
  end do

   
  !iL = maxloc(s(1:imin))
  smax=s(iL)
  xL = dble(iL);
  if (iL > 1) then
    s1 = 0.5d0*( s(iL+1)-s(iL-1) )
    s2 = s(iL-1)-2d0*s(iL)+s(iL+1)
    if (s2 /= 0d0) xL = xL-s1/s2
  endif
      
  !iR = maxloc(s(imin:n))
  smax=s(iR)
  xR = dble(iR);
  if (iR < n) then
    s1 = 0.5d0*( s(iR+1)-s(iR-1) )
    s2 = s(iR-1)-2d0*s(iR)+s(iR+1)
    if (s2 /= 0d0) xR = xR-s1/s2
  endif

  crack_size = xR-xL

end function crack_size


end module output
