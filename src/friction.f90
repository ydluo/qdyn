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

  case (2) ! 2018 SCEC Benchmark
    pb%theta_star = pb%dc/pb%v_star

  case (3) ! SEISMIC: the CNS friction law does not use theta_star
    pb%theta_star = 1

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
    mu = pb%mu_star - pb%a*log(pb%v_star/v) + pb%b*log(theta/pb%theta_star)

  case (1)
    mu = pb%mu_star - pb%a*log(pb%v1/v+1d0) + pb%b*log(theta/pb%theta_star+1d0)

  case (2) ! SCEC 2018 benchmark
    mu = pb%a*asinh( v/(2*pb%v_star)*exp( (pb%mu_star + pb%b*log(theta/pb%theta_star))/pb%a ) )  

  case (3) ! SEISMIC: CNS model
    write (6,*) "friction.f90::friction_mu is deprecated for the CNS model"
    stop

! new friction law:
!  case(xxx)
!    implement here your friction coefficient: mu = f(v,theta)
!    mu = ...

  case default
    stop 'friction_mu: unknown friction law type'
  end select

end function friction_mu

!--------------------------------------------------------------------------------------
subroutine dtheta_dt(v,tau,sigma,theta,theta2,dth_dt,dth2_dt,pb)

  use friction_cns, only : dphi_dt

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: v, tau, sigma
  double precision, dimension(pb%mesh%nn), intent(in) :: theta, theta2
  double precision, dimension(pb%mesh%nn) :: dth_dt, dth2_dt, omega

  ! SEISMIC: If the CNS model is selected
  if (pb%i_rns_law == 3) then
    call dphi_dt(v,tau,sigma,theta,theta2,dth_dt,dth2_dt,pb)
  ! SEISMIC: Else, the RSF model is selected (with various theta laws)
  else

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

  endif

end subroutine dtheta_dt

!--------------------------------------------------------------------------------------
subroutine dmu_dv_dtheta(dmu_dv,dmu_dtheta,v,tau,sigma,theta,theta2,pb)

  type(problem_type), intent(in) :: pb
  double precision, dimension(pb%mesh%nn), intent(in) :: v, tau, sigma, theta, theta2
  double precision, dimension(pb%mesh%nn), intent(out) :: dmu_dv, dmu_dtheta
  double precision :: z(pb%mesh%nn)

  select case (pb%i_rns_law)

  case(0)
    dmu_dtheta = pb%b / theta
    dmu_dv = pb%a / v

  case(1)
    dmu_dtheta = pb%b * pb%v2 / ( pb%v2*theta + pb%dc )
    dmu_dv = pb%a * pb%v1 / v / ( pb%v1 + v )
  
  case(2) ! 2018 SCEC Benchmark
    z = v/(2*pb%v_star)*exp( (pb%mu_star + pb%b*log(theta/pb%theta_star))/pb%a )
    dmu_dtheta = pb%b*v/(2*theta*pb%v_star)/sqrt(1 + z**2)*exp( (pb%mu_star + pb%b*log(theta/pb%theta_star))/pb%a )
    dmu_dv = pb%a/sqrt(1 + z**2)*exp( (pb%mu_star + pb%b*log(theta/pb%theta_star))/pb%a )/(2*pb%v_star)
 
  case(3) ! SEISMIC: CNS model
    write (6,*) "friction.f90::dmu_dv_dtheta is deprecated for the CNS model"
    stop

  case default
    write (6,*) "dmu_dv_dtheta: unkown friction law type"
    stop
  end select

end subroutine dmu_dv_dtheta

end module friction
