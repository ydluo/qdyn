! Auxiliary functions

module aux

  use problem_class, only : problem_type

  implicit none
  public

contains

!===============================================================================
! Helper routine to pack variables for solver
! Storage conventions:
!
! theta = yt(1::pb%neqs)
! v or tau = yt(2::pb%neqs) (depending on friction law)
! + feature specific variables
!
! dtheta/dt = dydt(1::pb%neqs)
! dv/dt or dtau/dt = dydt(2::pb%neqs)
! + feature specific variables
!===============================================================================
subroutine pack(yt, theta, main_var, sigma, theta2, P, pb)

  type(problem_type), intent(inout) :: pb
  double precision, dimension(pb%neqs*pb%mesh%nn), intent(out) :: yt
  double precision, dimension(pb%mesh%nn), intent(in) :: theta, main_var, sigma
  double precision, dimension(pb%mesh%nn), intent(in) :: theta2, P
  integer :: nmax, ind_stress_coupling, ind_localisation, ind_tp

  nmax = pb%neqs*pb%mesh%nn
  ! Define the indices of yt and dydt based on which
  ! features are requested (defined in input file)
  ind_stress_coupling = 2 + pb%features%stress_coupling
  ind_localisation = ind_stress_coupling + pb%features%localisation
  ind_tp = ind_localisation + pb%features%tp

  yt(1:nmax:pb%neqs) = theta
  yt(2:nmax:pb%neqs) = main_var

  if (pb%features%stress_coupling == 1) then
    yt(ind_stress_coupling:nmax:pb%neqs) = sigma
  endif

  if (pb%features%localisation == 1) then
    yt(ind_localisation:nmax:pb%neqs) = theta2
  endif

  if (pb%features%tp == 1) then
    yt(ind_tp:nmax:pb%neqs) = P
  endif

end subroutine pack

!===============================================================================
! Helper routine to unpack variables from solver
!===============================================================================
subroutine unpack(yt, theta, main_var, sigma, theta2, P, pb)

  type(problem_type), intent(inout) :: pb
  double precision, dimension(pb%neqs*pb%mesh%nn), intent(in) :: yt
  double precision, dimension(pb%mesh%nn) :: theta, main_var, sigma, theta2, P
  integer :: nmax, ind_stress_coupling, ind_localisation, ind_tp

  nmax = pb%neqs*pb%mesh%nn
  ind_stress_coupling = 2 + pb%features%stress_coupling
  ind_localisation = ind_stress_coupling + pb%features%localisation
  ind_tp = ind_localisation + pb%features%tp

  theta = yt(1:nmax:pb%neqs)
  main_var = yt(2:nmax:pb%neqs)

  if (pb%features%stress_coupling == 1) then
    sigma = yt(ind_stress_coupling:nmax:pb%neqs)
  else
    sigma = pb%sigma
  endif

  if (pb%features%localisation == 1) then
    theta2 = yt(ind_localisation:nmax:pb%neqs)
  else
    theta2 = 0d0
  endif

  if (pb%features%tp == 1) then
    P = yt(ind_tp:nmax:pb%neqs)
  else
    P = 0d0
  endif

end subroutine unpack

end module aux
