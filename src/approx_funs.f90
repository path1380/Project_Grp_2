module approx_funs
! ================================================================================================
! 
! File: approx_funs.f90
! Brief description: Contains function for evaluating approximate solution 
!
! Detailed description: This module contains a function for evaluating the approximate solution
!                       on a given grid, given the coefficients of the L2 projection of the 
!                       function onto the space of Legendre polynomials.
!
!
! Authors: Parth Thakkar, Fortino Garcia, Alex Rybchuk, William Boshell
!
! ================================================================================================
use type_defs
use leg_funs
implicit none

contains

    function approx_eval(lt_endpt, rt_endpt, num_grdpts,grdpts,leg_degree,leg_coeffs)
!=================================================================================================
!
! Inputs: - num_grdpts : Number of gridpoints appoximate
!                        solution is to be evaluated on
!         - grdpts     : array containing gridpoints
!         - leg_degree : highest degree of Legendre poly
!                        used   
!         - leg_coeffs : computed coefficients of 
!                        Legendre poly projection
!
! Output:   1D Array of approximate function values on the 
!           given grid
!
!=================================================================================================
        integer, intent(in) :: num_grdpts, leg_degree
        real(dp), intent(in) :: lt_endpt, rt_endpt
        real(dp), dimension(num_grdpts), intent(in) :: grdpts
        real(dp), dimension(0:leg_degree), intent(in) :: leg_coeffs
        real(dp), dimension(0:leg_degree) :: temp_array
        real(dp), dimension(num_grdpts) :: approx_eval
        real(dp) :: scaling, shift
        integer :: n
        
        !Calculate scaling and shift needed by affine map
        !Here we must map the given interval back to [-1,1]
        scaling = 1.0_dp/(rt_endpt - lt_endpt)
        shift = (rt_endpt + lt_endpt)

        do n = 1, num_grdpts
            !evaluate each Legendre poly at a given gridpoint
            temp_array = leg(leg_degree, scaling*(2.0_dp*grdpts(n) - shift))

            !build linear combination of Legendre polynomials
            approx_eval(n) = DOT_PRODUCT(leg_coeffs, temp_array)
        end do

    end function approx_eval
end module approx_funs
