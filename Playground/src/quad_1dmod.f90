!============================================================================================
!
!quad1dmod : This module defines a datatype quad_1d which stores values of different scalar/
!            vector fields in the given element. The values are evaluated at the Gauss-
!            Legendre quadrature nodes of the element. The name is not changed to quad_2d
!            in order to avoid conflict with other files.
!
!Variables : nvars - Number of dimensions of the vector field to be stored
!            qx - Maximum degree of Legendre polynomials in X-direction
!            qy - Maximum degree of Legendre polynomials in Y-direction
!            a_2d - We can store the values of the field as follows :
!                   Let C be a scalar field to be stored. Then we will have a variable C of 
!                   type quad_1d and the array a_2d will have values which follow:
!                   a(x_coord,y_coord,1) = value
!
!============================================================================================

module quad_1dmod
  use type_defs
  implicit none

  type quad_1d
     integer :: nvars, q, qx, qy
     real(kind=dp) :: lt_endpt, rt_endpt, lt_trace, rt_trace
     real(kind=dp), allocatable, dimension(:,:) :: a
     real(kind=dp), allocatable, dimension(:,:,:) :: a_2d
  end type quad_1d

contains
  
  subroutine allocate_quad1d(el)
    type(quad_1d)::el
    allocate(el%a(0:el%q,el%nvars))
    allocate(el%a_2d(0:el%qx,0:el%qy,0:el%nvars))
  end subroutine allocate_quad1d
  
  subroutine deallocate_quad1d(el)
    type(quad_1d)::el
    deallocate(el%a)
    deallocate(el%a_2d)
  end subroutine deallocate_quad1d

end module quad_1dmod
