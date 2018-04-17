module solve_quad
  use type_defs
  use quad_element
  use problemsetup
  implicit none
contains
  subroutine solve_quad(qd)
  type(quad), intent(inout) :: qd
  
  integer :: elt_num,elt_num_x,elt_num_y
  elt_num_x = ((elt_num-1)/ny) + 1
  elt_num_y = mod(elt_num-1,nx) + 1
  
end module solve_quad
