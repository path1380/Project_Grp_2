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

program main
  use InputControl
  use leg_funs
  use lgl
  use quad_1dmod
  use type_defs
  implicit none

!================================================================================================================
!
! 1) Evaluating the Gauss-Legendre Quadrature nodes and weights for the given maximum degree leg_degree.
!
! We define the degree 'leg_degree', weights 'leg_weights' and nodes 'leg_nodes'. 
!
! THE 'leg_degree' WILL BE GIVEN BY THE PERL SCRIPT TO THE TEMPLATE OF MAIN. 

! Then, memory is allocated to evaluate 'leg_weights' and 'leg_nodes', the function is called and the memory is 
! deallocated at the end of the program.
!
! DO NOT REMOVE THIS PART FROM THE CODE AS ALL THE FUNCTIONS ASSUME leg_nodes and leg_weights to be 
! known.
!
!===============================================================================================================

  integer :: leg_degree
  real(dp), dimension(:), allocatable :: leg_nodes, leg_weights

  leg_degree = 5

  allocate(leg_nodes(0:leg_degree), leg_weights(0:leg_degree))

  call lglnodes(leg_nodes, leg_weights, leg_degree)

!===============================================================================================================
!
! 2) The arrays are deallocated. (INCLUDING 'leg_weights' AND 'leg_nodes'). 
!
!===============================================================================================================


  deallocate(leg_nodes, leg_weights)


  write(*,*) -4*ATAN(-1.0_dp)
  write(*,*) exp(1.0_dp)
end program main
