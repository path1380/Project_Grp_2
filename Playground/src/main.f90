!================================================================================================================
!
! File: main.f90
! Brief description: Main program of the Project.
!
! Detailed description: This program will be regulated by a perl script which gives inputs to a template file. 
!                       The code is documented in a point-wise form. Please refer to the points in order to 
!                       understand what the heck is going on in the program.
!
! PLEASE MAINTAIN THE DOCUMENTATION FORMAT OF THE CODE. IT IS ABSOLUTELY CRITICAL FOR OTHERS TO UNDERSTAND WHAT
! YOU'VE DONE.
!
! Authors: Parth Thakkar, Fortino Garcia, Alex Rybchuk, William Boshell
!
!================================================================================================================

!TEMP_STEP
!Here we are trying to solve the equation ut + ux = 0 on the 2D domain [0,1]^2. Dirichlet BC u(0,y,t) = 0 and IC
!u(x,y,0) = u_init(x,y)

program main
  use problemsetup
  use legendre_module
  use quad_element
  use type_defs
  implicit none

  real(kind=dp) :: qnodes(0:nint),weights(0:nint)
  type(quad) :: qd

  integer :: i,j,n_gll

  call allocate_quad(qd,q,nint,1)
  
  n_gll = nint+1

  write(*,*) nx,ny
end program main
