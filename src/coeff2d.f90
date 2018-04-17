module coeff2d
  use type_defs
  use quad_element
  use problemsetup
  use legendre_module
  implicit none
contains
  subroutine set_initial_data(qd)
      ! ========================================================
      ! Inputs: - qd         : quad element containing the  
      !                        information of a given element
      !         - weights    : array containing the weights 
      !                        for gaussian quadrature
      !         - leg_mat    : matrix containing the evaluation
      !                        of each Leg poly at the quadrature
      !                        nodes
      !         - nint       : number of intervals given by 
      !                        xnodes (i.e. # of nodes - 1)   
      !
      ! Output:   Store the coefficients for the projection
      !           of initial data into the space of Legendre
      !           polynomials of degree q specified in 
      !           problemsetup.f90.
      ! ========================================================
      use type_defs
      use quad_element
      use problemsetup
      implicit none
      type(quad), intent(inout) :: qd
      integer :: i, j, n_gll, row, i1, i2
      integer :: INFO
      real(kind=dp) :: u_loc(0:nint,0:nint)
      real(kind=dp) :: weights(0:nint), xnodes(0:nint),diffmat(0:nint,0:nint),&
                   leg_mat(0:nint,0:q),leg_der_mat(0:nint,0:q),BFWeights(0:nint,2)
      real(kind=dp) :: b(0:(q+1)**2 - 1)

      n_gll = nint + 1
      b(:) = 0.0_dp
      !evaluate initial condition at nodes
      do j = 1,n_gll
        do i =1,n_gll
          u_loc(i-1,j-1) = init_u(qd%x(i,j),qd%y(i,j))
        end do 
      end do

      do j = 0,q
        do i = 0,nint
          leg_mat(i,j) = legendre(xnodes(i),j)
          leg_der_mat(i,j) = derlegendre(xnodes(i),j)
        end do
      end do

      !build RHS vector b
      do j = 0, q
        do i =0,q 
          row = i + j*(q+1)
          do i2 = 0,nint
            do i1 =0, nint
              b(row) = b(row) + leg_mat(i1,i)*leg_mat(i2,j)*weights(i1)&
                       *weights(i2)*u_loc(i1,i2)*qd%jac(i1+1,i2+1)
            end do 
          end do
        end do 
      end do

      !build matrices 
      CALL assemble(qd,nint,leg_mat,weights) 
      !here we'll need to backsolve the matrix to find the coefficients

      !build the LU decomposition of the mass matrix and 
      !backsolve for the coefficients
      call DGETRF((q+1)**2,(q+1)**2,qd%M,(q+1)**2,qd%IPIV,INFO)
      call DGETRS('N',(q+1)**2,1,qd%M,(q+1)**2,qd%IPIV,b,(q+1)**2,INFO)
      !Reshape and overwrite our quad with the coefficients
      do j = 0,q
        do i =0,q 
          qd%u(i,j,nvar) = b(i + j*(q+1))
        end do
      end do


    end subroutine set_initial_data
end module
