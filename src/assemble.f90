subroutine assemble(qd,nint,P,DERP,weights)
  use type_defs
  use quad_element
  use problemsetup, only : q
  implicit none
  integer :: nint
  type(quad) :: qd
  real(kind=dp) :: P(0:nint,0:q),DERP(0:nint,0:q) ! Legendre and derivative of L. at quadrature nodes
  real(kind=dp) :: fint(0:nint),weights(0:nint)
  integer :: i,j,k,l,iy,row,col
  !
  ! This routine assembles the matrices M and S in the system
  !
  ! M*u_t + Su = "flux-stuff",
  !
  ! We assume the modes are ordered in column major order with the "1" coordinate first:
  ! u_00,1, u_10,1, u_20,1,..., u_01,1,..., u_qq,1, u_00,2, u_10,2, u_20,2,..., u_01,2,..., u_qq,2....
  !
  ! Assemble Mass and Stiffness matrices
  ! i,k is index in r. j,l is index in s.
  ! i,j for phi
  ! k,l for u
  !
  ! M = [ \int \phi_{i,j} v_{k,l}
  qd%M(:,:) = 0.d0
  ! First diagonal block
  ! Row index (eq. number)
  do j = 0,q
   do i = 0,q
    row = i+1 + j*(q+1)
    do l = 0,q
     do k = 0,q
      col = k + 1 + l*(q+1)
      ! Integrate in r for each s
      do iy = 0,nint
       fint(iy) = sum(weights*qd%jac(:,iy)&
         *P(:,i)*P(iy,j)&   ! \phi-part, test
         *P(:,k)*P(iy,l))   !    v-part, trial
      end do
      ! Then integrate in s
      qd%M(row,col) = sum(weights*fint)
     end do
    end do
   end do
  end do

  qd%S(:,:) = 0.d0
  ! S. DEAA this is a \phi_x u_x term...
  do j = 0,q ! 2-dir phi
   do i = 0,q ! 1-dir
    row = i+1 + j*(q+1)
    do l = 0,q ! 2-dir u
     do k = 0,q ! 1-dir
      col = k + 1 + l*(q+1)
      ! Recall w_1 = r_1 w_r + s_1 w_s
      do iy = 1,nint
       fint(iy) = sum(weights*&
         qd%jac(:,iy)*&
         (qd%rx(:,iy)*DERP(:,i)*P(iy,j)&
         +qd%sx(:,iy)*P(:,i)*DERP(iy,j))*&
         (qd%rx(:,iy)*DERP(:,k)*P(iy,l)&
         +qd%sx(:,iy)*P(:,k)*DERP(iy,l)))
      end do
      qd%S(row,col) = sum(weights*fint)
      ! Also, w_y = r_y w_r +s_y w_s
      ! DEEA, not added yet...
     end do
    end do
   end do
  end do

end subroutine assemble