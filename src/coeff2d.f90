program coeff2d
!This is a test program to calculate the coefficients 
!of a 2D projection of initial data. Here we will choose
!the function 
!           u(x,y) = 0.25*(x^2 + y^2)
!so that 
!           u_xx + u_yy = 1
!Then we will attempt to compute the Laplacian and check our 
!work. First, we changed quad_element to be 0 indexed, let us 
!fix those bugs first.

  use type_defs
  use quad_element
  use problemsetup
  use legendre_module
  implicit none

  ! list of all elements
  type(quad), DIMENSION(:), ALLOCATABLE :: qds
  real(kind=dp), DIMENSION(:), ALLOCATABLE :: err_vec, err_vec_lap
  real(kind = dp), parameter :: pi = acos(-1.d0)

  integer :: i,j,ind, n_gll, k, l, INFO, step, n, neighbor_ind
  integer :: num_quads, sz2
  real(kind=dp) :: weights(0:nint), xnodes(0:nint),diffmat(0:nint,0:nint),&
                   leg_mat(0:nint,0:q),leg_der_mat(0:nint,0:q),BFWeights(0:nint,2)
  real(kind=dp) :: true_sol(0:nint,0:nint), approx_sol(0:nint,0:nint)
  real(kind=dp) :: u_x(0:(q+1)**2-1), gradx(0:(q+1)**2-1), grady(0:(q+1)**2-1)
  real(kind=dp) :: temp(0:(q+1)**2-1), hx, hy, lt_x, rt_x, lt_y, rt_y
  num_quads = nelemx*nelemy

  ! Weights for quadrature and differentiation on the elements.
  call lglnodes(xnodes,weights,nint)
  n_gll = nint + 1

  !build matrices with values of Legendre polys and 
  !their derivatives
  do j = 0,q
   do i = 0,nint
    leg_mat(i,j) = legendre(xnodes(i),j)
    leg_der_mat(i,j) = derlegendre(xnodes(i),j)
   end do
  end do

  ! Differentiation matrix for the metric.
  do i = 0,nint
   call weights1(xnodes(i),xnodes,nint,nint,1,BFWEIGHTS)
   DiffMat(i,:) = BFWEIGHTS(:,2)
  end do

  !allocate quad array
  ALLOCATE(qds(1:num_quads))
  ALLOCATE(err_vec(1:num_quads))
  ALLOCATE(err_vec_lap(1:num_quads))

  hx = 2.d0/DBLE(nelemx)
  hy = 2.d0/DBLE(nelemy)

  !loop and initialize our quad array
  do j = 1,nelemy
    
    !build corner of quad in y direction
    lt_y = -1 + DBLE(j-1)*hy
    rt_y = -1 + DBLE(j)*hy

    do i =1,nelemx

      !assign each quad an ID
      ind = i + (j-1)*nelemx
      call allocate_quad(qds(ind),q,n_gll,1)

      !insist on straight line boundaries
      qds(ind)%my_ind = ind 
      qds(ind)%bc_type(:) = 10

      !build corners of current quad in x direction
      lt_x = -1 + DBLE(i-1)*hx
      rt_x = -1 + DBLE(i)*hx

      !define corners of quad
      qds(ind)%xy(1,:) = (/rt_x, rt_y/)
      qds(ind)%xy(2,:) = (/lt_x, rt_y/)
      qds(ind)%xy(3,:) = (/lt_x, lt_y/)
      qds(ind)%xy(4,:) = (/rt_x, lt_y/)

      !compute and store the metric
      call set_metric(qds(ind),xnodes,diffmat,nint)
      call set_initial_data(qds(ind))

      !set information about neighbors and boundary 
      !conditions. NOTE: If instead of a neighbor 
      !the element has a physical boundary, we let
      !the value of the neighbor be -2 (later we'll)
      !replace with MPI_PROC_NULL

      !top neighbor
      neighbor_ind = ind + nelemx 
      qds(ind)%nbr(1,:) = (/neighbor_ind, 3/)
      qds(ind)%bc_type(1) = 1

      !bottom neighbor
      neighbor_ind = ind - nelemx 
      qds(ind)%nbr(3,:) = (/neighbor_ind, 1/)
      qds(ind)%bc_type(3) = 1

      !left neighbor
      neighbor_ind = ind - 1 
      qds(ind)%nbr(2,:) = (/neighbor_ind, 4/)
      qds(ind)%bc_type(2) = 1

      !right neighbor
      neighbor_ind = ind + 1 
      qds(ind)%nbr(4,:) = (/neighbor_ind, 2/)
      qds(ind)%bc_type(4) = 1

      !Default assume no boundary
      qds(ind)%has_physical_bc = .FALSE.

      ! ====== Account for physical boundaries ====== !
      
      !If we are on the bottom row
      IF( j .eq. 1) THEN
        !bottom "neighbor" is a boundary
        qds(ind)%nbr(3,:) = (/-2, 1/)
        qds(ind)%bc_type(3) = -1
        qds(ind)%has_physical_bc = .TRUE.
      END IF

      !If we are on the top row
      IF( j .eq. nelemy) THEN
        !top "neighbor" is a boundary
        qds(ind)%nbr(1,:) = (/-2, 3/)
        qds(ind)%bc_type(1) = -1
        qds(ind)%has_physical_bc = .TRUE.
      END IF

      !If we are on the left column
      IF( i .eq. 1) THEN
        !left "neighbor" is a boundary
        qds(ind)%nbr(2,:) = (/-2, 4/)
        qds(ind)%bc_type(2) = -1
        qds(ind)%has_physical_bc = .TRUE.
      END IF

      !If we are on the right column
      IF( i .eq. nelemx) THEN
        !right "neighbor" is a boundary
        qds(ind)%nbr(4,:) = (/-2, 2/)
        qds(ind)%bc_type(4) = -1
        qds(ind)%has_physical_bc = .TRUE.
      END IF
    end do
  end do 


  !Now let's double check the projection is 
  !working correctly for multiple elements
  do n = 1,num_quads

    !build coefficient vector from matrix version
    do j = 0,q
      do i =0,q 
        u_x(i + j*(q+1)) = qds(n)%u(i,j,nvar)
      end do
    end do

    !build true solution on the given grid
    do j = 0,n_gll-1
      do i =0,n_gll-1
        true_sol(i,j) = init_u(qds(n)%x(i,j),qds(n)%y(i,j))

        !build approximation
        approx_sol(i,j) = 0.0_dp
        do l=0,q
          do k=0,q 
            approx_sol(i,j) = approx_sol(i,j) + qds(n)%u(k,l,nvar)*leg_mat(i,k)*leg_mat(j,l)
          end do 
        end do 
      end do
    end do
    err_vec(n) = NORM2(true_sol - approx_sol) 
  end do
  write(*,*) 'Here is the error for projecting down to Legendre polys: \n'
  write(*,*) MAXVAL(err_vec)

  !Needed for DGEMV
  sz2 = (q+1)**2
  step = 1
 
  !Now let's compute the Laplacian in each quad
  do n = 1,num_quads
    
    !build coefficient vector from matrix version
    do j = 0,q
      do i =0,q 
        u_x(i + j*(q+1)) = qds(n)%u(i,j,nvar)
      end do
    end do

    !take a single derivative in x
    call DGEMV('N',sz2,sz2,1.0_dp,qds(n)%Diff_x,sz2,&
         u_x,step,0.0_dp,gradx,step)
    call DGETRS('N',(q+1)**2,1,qds(n)%M,(q+1)**2,qds(n)%IPIV,gradx,(q+1)**2,INFO)
    
    temp = 0.0_dp
    !take a second derivative in x
    call DGEMV('N',sz2,sz2,1.0_dp,qds(n)%Diff_x,sz2,&
         gradx,step,0.0_dp,temp,step)
    call DGETRS('N',(q+1)**2,1,qds(n)%M,(q+1)**2,qds(n)%IPIV,temp,(q+1)**2,INFO)
    gradx = temp

    !take a single derivative in y
    call DGEMV('N',sz2,sz2,1.0_dp,qds(n)%Diff_y,sz2,&
         u_x,step,0.0_dp,grady,step)
    call DGETRS('N',(q+1)**2,1,qds(n)%M,(q+1)**2,qds(n)%IPIV,grady,(q+1)**2,INFO)
    temp = 0.0_dp
    !take a second derivative in x
    call DGEMV('N',sz2,sz2,1.0_dp,qds(n)%Diff_y,sz2,&
         grady,step,0.0_dp,temp,step)
    call DGETRS('N',(q+1)**2,1,qds(n)%M,(q+1)**2,qds(n)%IPIV,temp,(q+1)**2,INFO)
    grady = temp

    !Let us now check the derivative
    !build true solution on the given grid
    ! true_sol = 1.d0
    do j = 0,n_gll-1
      do i =0,n_gll-1

        true_sol(i,j) = -2.0_dp*(pi**2.0_dp)*init_u(qds(n)%x(i,j),qds(n)%y(i,j))
        !build approximation
        approx_sol(i,j) = 0.0_dp
        do l=0,q
          do k=0,q 
            approx_sol(i,j) = approx_sol(i,j) + (gradx(k + l*(q+1)) + grady(k + l*(q+1)))*leg_mat(i,k)*leg_mat(j,l)
          end do 
        end do 
      end do 
    end do
    err_vec_lap(n) = NORM2(true_sol - approx_sol)
    gradx = 0.0_dp
    grady = 0.0_dp
  end do

  write(*,*) 'Here is the error in computing the Laplacian: \n'
  write(*,*) MAXVAL(err_vec_lap)


  !Here we should now have a subroutine that computes the Laplacian 
  !and here we'll have a time step loop.

  !We have an initial approximation of the Laplacian, now 
  !we will do a single RK4 time step and see if we can get
  !that with decent accuracy first. Then we'll time step 
  !more generally.

  do i = 1,num_quads
    CALL deallocate_quad(qds(i))
  end do
  DEALLOCATE(qds)
  DEALLOCATE(err_vec)


  ! !build coefficient vector from matrix version
  ! do j = 0,q
  !   do i =0,q 
  !     u_x(i + j*(q+1)) = qds(ind)%u(i,j,nvar)
  !   end do
  ! end do

  ! sz2 = (q+1)**2
  ! step = 1

  ! !take a single derivative in x
  ! call DGEMV('N',sz2,sz2,1.d0,qds(ind)%Diff_x,sz2,&
  !      u_x,step,0.d0,gradx,step)
  ! call DGETRS('N',(q+1)**2,1,qds(ind)%M,(q+1)**2,qds(ind)%IPIV,gradx,(q+1)**2,INFO)
  ! temp = 0.d0
  ! !take a second derivative in x
  ! call DGEMV('N',sz2,sz2,1.d0,qds(ind)%Diff_x,sz2,&
  !      gradx,step,0.d0,temp,step)
  ! call DGETRS('N',(q+1)**2,1,qds(ind)%M,(q+1)**2,qds(ind)%IPIV,temp,(q+1)**2,INFO)
  ! gradx = temp

  ! !take a single derivative in y
  ! call DGEMV('N',sz2,sz2,1.d0,qds(ind)%Diff_y,sz2,&
  !      u_x,step,0.d0,grady,step)
  ! call DGETRS('N',(q+1)**2,1,qds(ind)%M,(q+1)**2,qds(ind)%IPIV,grady,(q+1)**2,INFO)
  ! temp = 0.d0
  ! !take a second derivative in x
  ! call DGEMV('N',sz2,sz2,1.d0,qds(ind)%Diff_y,sz2,&
  !      grady,step,0.d0,temp,step)
  ! call DGETRS('N',(q+1)**2,1,qds(ind)%M,(q+1)**2,qds(ind)%IPIV,temp,(q+1)**2,INFO)
  ! grady = temp

  ! !Let us now check the derivative
  ! !build true solution on the given grid
  ! true_sol = 1.d0
  ! do j = 0,n_gll-1
  !   do i =0,n_gll-1
  !     ! true_sol(i,j) = 12.d0*(qds%x(i,j)**2.d0)*(qds%y(i,j)**3.d0)

  !     !build approximation
  !     approx_sol(i,j) = 0.0_dp
  !     do l=0,q
  !       do k=0,q 
  !         approx_sol(i,j) = approx_sol(i,j) + (gradx(k + l*(q+1)) + grady(k + l*(q+1)))*leg_mat(i,k)*leg_mat(j,l)
  !       end do 
  !     end do 
  !   end do 
  ! end do

  ! write(*,*) MAXVAL(ABS(true_sol - approx_sol))
  ! call deallocate_quad(qds(ind))

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
      real(kind=dp) :: b(0:(q+1)**2 - 1)

      n_gll = nint + 1
      b(:) = 0.0_dp

      !evaluate initial condition at nodes
      do j = 0,n_gll-1
        do i =0,n_gll-1
          u_loc(i,j) = init_u(qd%x(i,j),qd%y(i,j))
        end do 
      end do

      !build RHS vector b
      do j = 0, q !along s
        do i =0,q !along r
          row = i + j*(q+1)
          do i2 = 0,nint
            do i1 =0, nint
              b(row) = b(row) + leg_mat(i1,i)*leg_mat(i2,j)*weights(i1)&
                      *weights(i2)*u_loc(i1,i2)*qd%jac(i1,i2)
            end do 
          end do
        end do 
      end do

      !build matrices 
      CALL assemble(qd,nint,leg_mat,leg_der_mat,weights) 

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

    subroutine set_metric(qd,xnodes,diffmat,nint)
      ! ========================================================
      ! Inputs: - qd         : quad element containing the  
      !                        geometrical information of a 
      !                        given element
      !         - xnodes     : array containing the Legendre
      !                        quadrature nodes
      !         - diffmat    : matrix containing derivative 
      !                        approximations (from weights.f)
      !                        to compute the metric
      !         - nint       : number of intervals given by 
      !                        xnodes (i.e. # of nodes - 1)   
      !
      ! Output:   Store the local coordinates, unit normals,
      !           and metric for a given quad element.
      ! ========================================================
      use type_defs
      use quad_element
      use problemsetup, only: pis
      implicit none
      type(quad) :: qd
      integer :: nint
      integer :: ix,iy
      real(kind=dp) :: xnodes(0:nint)
      real(kind=dp) :: x_coord_elem(0:nint,0:nint) ,y_coord_elem(0:nint,0:nint),diffmat(0:nint,0:nint)
      real(kind=dp) :: pi1(2),pi2(2),pi3(2),pi4(2),pi2_m(2),pi2_p(2),pi4_m(2),pi4_p(2)
      real(kind=dp) :: xy_loc(2),xy_s(2),xy_e(2),eta,xi

      ! Compute metric
      ! We use a Gordon-Hall mapping
      ! The soubroutine pis must contain the approproate information
      ! for the parametrization of the curvilinear elements
      !
      xy_s = qd%xy(3,1:2)
      xy_e = qd%xy(2,1:2)
      eta = 1.d0
      call pis(pi2_p,eta,xy_s,xy_e,qd%bc_type(2))
      eta = -1.d0
      call pis(pi2_m,eta,xy_s,xy_e,qd%bc_type(2))
      !
      xy_s = qd%xy(4,1:2)
      xy_e = qd%xy(1,1:2)
      eta = 1.d0
      call pis(pi4_p,eta,xy_s,xy_e,qd%bc_type(4))
      eta = -1.d0
      call pis(pi4_m,eta,xy_s,xy_e,qd%bc_type(4))
      !
      do iy = 0,nint
         eta = xnodes(iy)
         !
         xy_s = qd%xy(3,1:2)
         xy_e = qd%xy(2,1:2)
         call pis(pi2,eta,xy_s,xy_e,qd%bc_type(2))
         !
         xy_s = qd%xy(4,1:2)
         xy_e = qd%xy(1,1:2)
         call pis(pi4,eta,xy_s,xy_e,qd%bc_type(4))
         do ix = 0,nint
            xi  = xnodes(ix)
            !
            xy_s = qd%xy(2,1:2)
            xy_e = qd%xy(1,1:2)
            call pis(pi1,xi,xy_s,xy_e,qd%bc_type(1))
            !
            xy_s = qd%xy(3,1:2)
            xy_e = qd%xy(4,1:2)
            call pis(pi3,xi,xy_s,xy_e,qd%bc_type(3))
            xy_loc = (1.d0-eta)/2.d0*pi3+(1.d0+eta)/2.d0*pi1&
                 +(1.d0-xi)/2.d0*(pi2-(1.d0+eta)/2.d0*pi2_p-(1.d0-eta)/2.d0*pi2_m)&
                 +(1.d0+xi)/2.d0*(pi4-(1.d0+eta)/2.d0*pi4_p-(1.d0-eta)/2.d0*pi4_m)
            x_coord_elem(ix,iy) = xy_loc(1)
            y_coord_elem(ix,iy) = xy_loc(2)
         end do
      end do

      !Store local coordinates in quad
      qd%x = x_coord_elem
      qd%y = y_coord_elem

      !Compute the derivatives and Jacobian of the metric
      call compute_curve_metric(qd%rx,qd%sx,qd%ry,qd%sy,qd%jac,&
           x_coord_elem,y_coord_elem,Diffmat,nint)
      
      ! ======= Compute normals and line elements on all sides ====== !

      ! Face 1. corresponds to s = 1 and r \in [-1,1].
      ! Thus the normal is (s_x,s_y) / \sqrt(s_x^2+s_y^2).
      ! The line integral element is dl = \sqrt(x_r^2+y_r^2)| = J * \sqrt(s_x^2+s_y^2).
      qd%dl_face(:,1) = sqrt(qd%sx(:,n_gll-1)**2+qd%sy(:,n_gll-1)**2)
      ! Compute outward pointing unit normal.
      qd%nx_in(:,1)   = qd%sx(:,n_gll-1)/qd%dl_face(:,1)
      qd%ny_in(:,1)   = qd%sy(:,n_gll-1)/qd%dl_face(:,1)
      ! Scale by Jacobian to get the metric.
      qd%dl_face(:,1) = qd%dl_face(:,1)*qd%jac(:,n_gll-1)

      ! Face 2. corresponds to r = -1 and s \in [-1,1].
      ! Thus the normal is (-r_x,-r_y) / \sqrt(r_x^2+r_y^2).
      ! The line integral element is dl = \sqrt(x_s^2+y_s^2)| = J * \sqrt(r_x^2+r_y^2).
      qd%dl_face(:,2) = sqrt(qd%rx(0,:)**2+qd%ry(0,:)**2)
      qd%nx_in(:,2)   = -1.0_dp*qd%rx(0,:)/qd%dl_face(:,2)
      qd%ny_in(:,2)   = -1.0_dp*qd%ry(0,:)/qd%dl_face(:,2)
      qd%dl_face(:,2) = qd%dl_face(:,2)*qd%jac(0,:)

      ! Face 3. corresponds to s = -1 and r \in [-1,1].
      ! Thus the normal is (-s_x,-s_y) / \sqrt(s_x^2+s_y^2).
      ! The line integral element is dl = \sqrt(x_r^2+y_r^2)| = J * \sqrt(s_x^2+s_y^2).
      qd%dl_face(:,3) = sqrt(qd%sx(:,0)**2+qd%sy(:,0)**2)
      qd%nx_in(:,3)   = -1.0_dp*qd%sx(:,0)/qd%dl_face(:,3)
      qd%ny_in(:,3)   = -1.0_dp*qd%sy(:,0)/qd%dl_face(:,3)
      qd%dl_face(:,3) = qd%dl_face(:,3)*qd%jac(:,0)

      ! Face 4. corresponds to r = 1 and s \in [-1,1].
      ! Thus the normal is (r_x,r_y) / \sqrt(r_x^2+r_y^2).
      ! The line integral element is dl = \sqrt(x_s^2+y_s^2)| = J * \sqrt(r_x^2+r_y^2).
      qd%dl_face(:,4) = sqrt(qd%rx(n_gll-1,:)**2+qd%ry(n_gll-1,:)**2)
      qd%nx_in(:,4)   = qd%rx(n_gll-1,:)/qd%dl_face(:,4)
      qd%ny_in(:,4)   = qd%ry(n_gll-1,:)/qd%dl_face(:,4)
      qd%dl_face(:,4) = qd%dl_face(:,4)*qd%jac(n_gll-1,:)
    end subroutine set_metric

    subroutine compute_curve_metric(r_x,s_x,r_y,s_y,jac,X,Y,D,n)
      ! ========================================================
      ! Inputs: - X         : local X coordinates in quadrature
      !                       grid
      !         - Y         : local Y coordinates in grid
      !         - D         : Differentiation matrix given
      !                       by Bengt Fornberg's weights.f
      !         - n         : number of intervals given by 
      !                       xnodes (i.e. # of nodes - 1) 
      !
      ! Outputs: Metric information, in particular r_x, s_x,
      !          r_y, s_y, and the Jacobian of the transformation
      ! ========================================================
      use type_defs
      implicit none
      integer :: n
      real(kind=dp), dimension(0:n,0:n) :: r_x,s_x,r_y,s_y,jac,X,Y,D
      real(kind=dp), dimension(0:n,0:n) :: x_r, x_s, y_r, y_s

      integer :: i
      !% Compute the derivatives w.r.t r & s
      do i = 0,n
       x_r(:,i) = matmul(D,X(:,i))
       y_r(:,i) = matmul(D,Y(:,i))
       x_s(i,:) = matmul(D,X(i,:))
       y_s(i,:) = matmul(D,Y(i,:))
      end do
      jac = x_r*y_s-y_r*x_s
      r_x =  y_s/jac
      r_y = -x_s/jac
      s_x = -y_r/jac
      s_y =  x_r/jac

    end subroutine compute_curve_metric

    subroutine assemble(qd,nint,P,DERP,weights)
      ! ========================================================
      ! Inputs: - qd         : quad element containing the  
      !                        information of a given element
      !         - nint       : number of intervals given by 
      !                        xnodes (i.e. # of nodes - 1) 
      !         - P          : matrix containing the evaluation
      !                        of each Leg poly at the quadrature
      !                        nodes
      !         - DERP       : matrix containing the evaluation
      !                        of derivatives of each Leg poly
      !         - weights    : array containing the weights 
      !                        for gaussian quadrature
      !
      ! Output:   Build and store the mass and differentiation
      !           matrices for each element in the grid
      ! ========================================================
      use type_defs
      use quad_element
      use problemsetup, only : q
      implicit none
      integer :: nint
      type(quad) :: qd
      real(kind=dp) :: P(0:nint,0:q) !Legendre polys at quadrature nodes
      real(kind=dp) :: DERP(0:nint,0:q) ! Legendre and derivative of L. at quadrature nodes
      real(kind=dp) :: fint(0:nint),weights(0:nint),fint_x(0:nint),fint_y(0:nint)
      integer :: i,j,k,l,iy,row,col
      !
      ! This routine assembles the mass matrix M
      !
      ! We assume the modes are ordered in column major order with the "1" coordinate first:
      ! u_00,1, u_10,1, u_20,1,..., u_01,1,..., u_qq,1, u_00,2, u_10,2, u_20,2,..., u_01,2,..., u_qq,2....
      !
      ! Assemble Mass and Stiffness matrices
      ! i,k is index in r. j,l is index in s.
      ! i,j for phi
      ! k,l for u
      qd%M(:,:) = 0.0_dp
      qd%Diff_x(:,:) = 0.0_dp
      qd%Diff_y(:,:) = 0.0_dp
      do j = 0,q
       do i = 0,q
        row = i + j*(q+1)
        do l = 0,q    !track degree in y direciton
         do k = 0,q   !track degree in x direction
          col = k + l*(q+1)
          ! Integrate in r for each s
          do iy = 0,nint

           !Mass matrix quadrature in r
           fint(iy) = sum(weights*qd%jac(:,iy)&
             *P(:,i)*P(iy,j)&
             *P(:,k)*P(iy,l))  
           
           !Differentiation matrix quadrature in r
           !Here we differentiate in x
           fint_x(iy) = sum(weights*qd%jac(:,iy)&
             *P(iy,i)*P(:,j)*(DERP(iy,k)*P(:,l)*qd%rx(:,iy)+& 
             DERP(:,l)*P(iy,k)*qd%sx(:,iy))) 

           !Differentiation matrix quadrature in r
           !Here we differentiate in y
           fint_y(iy) = sum(weights*qd%jac(:,iy)&
             *P(iy,i)*P(:,j)*(DERP(iy,k)*P(:,l)*qd%ry(:,iy)+& 
             DERP(:,l)*P(iy,k)*qd%sy(:,iy)))   
          end do

          ! Then integrate in s
          qd%M(row,col) = sum(weights*fint)
          qd%Diff_x(row,col) = sum(weights*fint_x)
          qd%Diff_y(row,col) = sum(weights*fint_y)
          !Note that if we have to do an integration
          !by parts, we'll need to transpose and
          !multiply by 1.
         end do
        end do
       end do
      end do
    end subroutine assemble


end program coeff2d
