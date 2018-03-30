module coeff
!==================================================================================================================
!
! File: coeff.f90
! Brief description: Evaluates the coefficients of the function approximation as a linear combination of Legendre
!                    polynomials.
!
! Detailed description: This module defines the function "element" which computes the coefficients of the L2 
!                       projection of a given function onto the space of Legendre polynomials for a given interval
!
! Authors: Parth Thakkar, Fortino Garcia, Alex Rybchuk, William Boshell      
!==================================================================================================================
  use type_defs
  use InputControl
  use lgl
  use quad_1dmod
  use leg_funs
  implicit none

contains


  function element(lt_endpt,rt_endpt)
!==================================================================================================================
! Inputs: - lt_endpt   : location of the left endpoint of 
!                        the given sub-interval
!         - rt_endpt   : location of the right endpoint
!                        of the given sub-interval
!         - leg_degree : highest degree of Legendre poly
!                        used   
!
! Output:   quad_1d data type containing the coefficients
!           of the Legendre polynomial L2 projection in 
!           the given subinterval
!==================================================================================================================
    real(dp), intent(in) :: lt_endpt, rt_endpt
    !integer, intent(in) :: leg_degree
    real(dp), dimension(0:leg_degree) :: temp_array, leg_nodes, leg_weights, fun_vals 
    integer :: i

    !Declare our quad_1d element and allocate its memory
    type(quad_1d) :: element
    
    element%q=leg_degree
    element%lt_endpt = lt_endpt
    element%rt_endpt = rt_endpt
    element%nvars=1    !NOTE: We'll need to account for multiple dimensions later
    call allocate_quad1d(element)
    
    !Pre-assign a value of zero to each entry in our quad_1d element
    do i=0, leg_degree
      element%a(i,element%nvars) = 0
    end do

    !Generate quadrature weights and nodes
    !call lglnodes(leg_nodes,leg_weights,leg_degree)

    !Evaluate the given function at each quadrature node (here the nodes are mapped to the current
    !interval and passed into function_eval to save memory). The function can be modified in InputControl.f90
    fun_vals = function_eval(leg_degree+1, 0.5_dp*(leg_nodes*(rt_endpt - lt_endpt) + (rt_endpt + lt_endpt))) 

    !Loop over each node and build our integral approximation 
    do i=0,leg_degree
      temp_array = leg(leg_degree, leg_nodes(i))
      element%a(:,element%nvars) = element%a(:,element%nvars) + temp_array*leg_weights(i)*fun_vals(i)
    end do

    !Solve for coefficients of L2 projection of the given
    !function to the space of Legendre polynomials
    do i=0,leg_degree
      element%a(i,element%nvars) = 0.5_dp*(2.0_dp*dble(i) + 1.0_dp)*element%a(i,element%nvars)
    end do

  end function element

end module coeff
!--------------------------------------------------------------------------------------------------------------
!
!Potential Improvements : 1. Find a way to use and call q(i). 2. Find a way to dynamically deallocate memory
!Aim : To be able to run the code for any case without editing any file other than InputControl.f90
!
!--------------------------------------------------------------------------------------------------------------
