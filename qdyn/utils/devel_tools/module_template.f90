module mod_name

! MOD_NAME: 
!

  use some_module_1
  use some_module_2

  implicit none
  private

  type type_name
    private
    double precision, pointer ::
    integer, pointer ::
    logical, pointer ::
    double precision ::
    integer ::
    logical ::
  end type type_name

  public :: type_name,sub_name

contains

!=====================================================================
! SUB_NAME:
!

subroutine sub_name()

  use some_module, only : some_variable

  type(), intent() ::
  double precision, intent() ::
  integer, intent() ::
  logical, intent() ::

  type() ::
  double precision ::
  integer ::
  logical ::

end subroutine

end module mod_name
