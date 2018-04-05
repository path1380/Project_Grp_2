module problemsetup
  use type_defs
  implicit none
  integer, parameter :: nvar = 1
  integer, parameter :: q = 1
  integer, parameter :: nint = 10
  !number of intervals in theta and r
  integer, parameter :: nt = 2, nr = 2
  real(kind = dp), parameter :: CFL = 0.1d0
  real(kind = dp), parameter :: tend = 1.d0
  real(kind = dp) :: bc(10:99,nvar)

  integer, parameter :: nplot = 4
  logical, parameter :: upwind = .true.
  logical, parameter :: plot = .true.

contains

  subroutine set_bc
    implicit none
    !
    ! This routine is used to set boundary conditions
    ! on boundary curve xx
    !
    ! Default is Dirichlet for all boundaries.
    bc = 1.d0
    ! DEAA: FILL IN

  end subroutine set_bc

  real(kind = dp) function init_u(x,y)
    use type_defs
    implicit none
    real(kind = dp) :: x,y
    real(kind = dp), parameter :: pi = acos(-1.d0)
    !init_u = sin(2.d0*pi*x)*sin(2.d0*pi*y)
    init_u = x*y + 2.0_dp
    return
    ! DEAA: Change to fit with your problem
  end function init_u

  subroutine pis(xy,s,xy_start,xy_end,curve_type)
    use type_defs
    implicit none
    real(kind=dp) :: xy(2),xy_start(2),xy_end(2),s
    real(kind=dp) :: theta_s, theta_e
    integer :: curve_type
    if (curve_type .eq. 10 ) then
     ! Straight lines.
     xy = xy_start + 0.5d0*(1.d0 + s)*(xy_end - xy_start)
    elseif (curve_type .eq. 12) then
     ! Circle of radius 0.5 and center in (0,0).
     theta_s = atan2(xy_start(2),xy_start(1))
     theta_e = atan2(xy_end(2),xy_end(1))
     xy(1) = 0.5d0*cos(theta_s + 0.5d0*(1.d0 + s)*(theta_e - theta_s))
     xy(2) = 0.5d0*sin(theta_s + 0.5d0*(1.d0 + s)*(theta_e - theta_s))
    elseif (curve_type .eq. 14) then
     ! Circle of radius 1 and center in (0,0).
     theta_s = atan2(xy_start(2),xy_start(1))
     theta_e = atan2(xy_end(2),xy_end(1))
     xy(1) = cos(theta_s + 0.5d0*(1.d0 + s)*(theta_e - theta_s))
     xy(2) = sin(theta_s + 0.5d0*(1.d0 + s)*(theta_e - theta_s))
    elseif (curve_type .eq. 100 ) then
     ! Straight line for internal boundaries, used for the mapping of
     ! curved elements.
     xy = xy_start + 0.5d0*(1.d0 + s)*(xy_end - xy_start)
    end if
  end subroutine pis

end module problemsetup
