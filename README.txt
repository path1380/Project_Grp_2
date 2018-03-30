This README file accompanies the solver code for the project. PLEASE ADHERE TO THE CODING AND DOCUMENTING STANDARDS.
The standard to be used is:

!===================================================================================================================
!
! File: <Filename>
! Brief description: <Any dependencies on perl scripts to be mentioned here along with other things.>
!
! Detailed description: 
!
! Inputs: 
! Outputs:
!
! Authors:
!
!===================================================================================================================

Points to be noted:

-> The 'Makefile_test' file is used to test the different features of the code. We have 4 .PHONY commands:
(a) make compile -f <Makefile_name> : This command compiles the code.
(b) make build -f <Makefile_name> : This command builds the code.
(c) make run -f <Makefile_name> : This command runs the code.
(d) make clean -f <Makefile_name> : Self-explanatory!

-> The perl script will talk directly to the code, hence all inputs shall be given to the perl script.
   The perl script will talk to the files 'InputControl.f90' and 'main.f90'. The critical parameters 
   of the problem will be given to 'main.f90' whereas parameters needed for certain functions to run,
   i.e. those which will define certain functions which in turn would be inputs for the problem
   would be given to 'InputControl.f90'. Thus, templates for these 2 files will be made.

   (The problem is that certain parameters cannot be directly defined in module files.)

-> The 'Playground' directory is used for testing the code. PLEASE DO NOT TAMPER WITH THE CODE OUTSIDE THE 'Playground'
   DIRECTORY UNLESS YOU ARE ABSOLUTELY THOROUGHLY CONFIRMEDLY CONVINCED THAT IT WORKS. PLEASE DO THE SMALLEST OF 
   INCREMENTAL TESTING IN THIS DIRECTORY ONLY. USE THE MAIN REPO DIRECTORY FOR MAIN PROGRAM EXECUTION, POST-PROCESSING 
   AND STUFF OTHER THAN TESTING. Studying the results is an important part of our project and we need this provision
   to do this without interruption. If we find that something is wrong with the code during post-processing, Git's
   got our back!
