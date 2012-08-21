module friction

! This is the only module that needs modifications to implement a new friction law 
! Follow the instructions in the comment blocks below that start with "! new friction law:"
!
! Assumptions:
!   The friction coefficient (mu) and the state variable rate (dtheta/dt)
!   depend on slip velocity (v) and state variable (theta)
!     mu = f(v,theta)
!     dtheta/dt = g(v,theta)
!   All friction properties can be spatially non-uniform

  use problem_class, only : problem_type

  implicit none
  private

  public  :: set_theta_star, friction_mu, dmu_dv_dtheta, dtheta_dt

contains

!--------------------------------------------------------------------------------------
subroutine set_theta_star(pb)

  type(problem_type), intent(inout) :: pb

  select case (pb%i_rns_law)

  case (0)
    pb%theta_star = pb%dc/pb%v_star

  case (1)
    pb%theta_star = pb%dc/pb%v2

! new friction law:
!  case(xxx)
!    implement here your definition of theta_star (could be none)
!    pb%theta_star = ...

  case default 
    stop 'set_theta_star: unknown friction law type'
  end select

end subroutine set_theta_star

!--------------------------------------------------------------------------------------
function friction_mu(v,theta,pb) result(mu)

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: v, theta
  double precision, dimension(pb%mesh%nn) :: mu

  select case (pb%i_rns_law)

  case (0)
    mu = pb%mu_star - pb%a*log(pb%v_star/pb%v) + pb%b*log(pb%theta/pb%theta_star)

  case (1)
    mu = pb%mu_star - pb%a*log(pb%v1/pb%v+1d0) + pb%b*log(pb%theta/pb%theta_star+1d0)

! new friction law:
!  case(xxx)
!    implement here your friction coefficient: mu = f(v,theta)
!    mu = ...

  case default 
    stop 'friction_mu: unknown friction law type'
  end select

end function friction_mu
  
!--------------------------------------------------------------------------------------
function dtheta_dt(v,theta,pb) result(dth_dt)

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: v, theta

  double precision, dimension(pb%mesh%nn) :: dth_dt, omega
  
  omega = v * theta / pb%dc
  select case (pb%itheta_law)

  case(0) ! "aging" in the no-healing approximation
    dth_dt = -omega

  case(1) ! "aging" law
    dth_dt = 1.d0-omega

  case(2) ! "slip" law
    dth_dt = -omega*log(omega)

! new friction law:
!  case(xxx)
!    implement here your state evolution law: dtheta/dt = g(v,theta)
!    dth_dt = ...

  case default
    stop 'dtheta_dt: unknown state evolution law type'
  end select

end function dtheta_dt

!--------------------------------------------------------------------------------------
subroutine dmu_dv_dtheta(dmu_dv,dmu_dtheta,v,theta,pb)

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: v, theta
  double precision, dimension(pb%mesh%nn), intent(out) :: dmu_dv,dmu_dtheta

  double precision, dimension(pb%mesh%nn) :: dth_dt, omega

  select case (pb%i_rns_law)

  case(0)
    dmu_dtheta = pb%b / theta
    dmu_dv = pb%a / v

  case(1)
    dmu_dtheta = pb%b * pb%v2 / ( pb%v2*theta + pb%dc )
    dmu_dv = pb%a * pb%v1 / v / ( pb%v1 + v )

! new friction law:
!  case(xxx)
!    implement here the partial derivatives of the friction coefficient
!    dmu_dtheta = ...
!    dmu_dv = ...

  case default
    stop 'dmu_dv_dtheta: unknown friction law type'
  end select
  
end subroutine dmu_dv_dtheta

end module friction
