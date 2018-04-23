module solve_quad_flux
  use type_defs
  use quad_element
  use problemsetup
  implicit none
contains
  subroutine solve_flux(qd,beta)
     type(quad), intent(inout) :: qd
     real(kind=dp) :: beta
  
     integer :: elt_num,elt_num_x,elt_num_y

     integer :: i,j,k

     elt_num = qd%my_ind
     elt_num_x = ((elt_num-1)/ny) + 1
     elt_num_y = ny - mod(elt_num-1,nx)

     do i=1,4
        do j=1,2
           qd%nbr(i,j) = 0.0_dp
        end do
        do k=1,qd%n_gll
           qd%u_out(k,i,qd%nvar) = 0.0_dp
        end do
     end do

     do k=1,qd%n_gll
        qd%u_in(k,1,qd%nvar) = qd%u(k-1,qd%q,qd%nvar)
        qd%u_in(k,2,qd%nvar) = qd%u(0,k-1,qd%nvar)
        qd%u_in(k,3,qd%nvar) = qd%u(k-1,0,qd%nvar)
        qd%u_in(k,4,qd%nvar) = qd%u(qd%q,k-1,qd%nvar)
     end do

     if (elt_num_y /= ny) then
        qd%nbr(1,1) = elt_num - 1
        qd%nbr(1,2) = 3
        !write neighbouring element solutions on qd%u_out. Then, calculate the qd%dl_face using beta. Do the same in other conditions.
     end if

     if (elt_num_x /= 1) then
        qd%nbr(2,1) = elt_num - ny
        qd%nbr(2,2) = 4
     end if

     if (elt_num_y /= 1) then
        qd%nbr(3,1) = elt_num + 1
        qd%nbr(3,2) = 1
     end if

     if (elt_num_x /= nx) then
        qd%nbr(4,1) = elt_num + ny
        qd%nbr(4,2) = 2
     end if

     write(*,*) maxval(qd%u_in)
  end subroutine solve_flux
  
end module solve_quad_flux
