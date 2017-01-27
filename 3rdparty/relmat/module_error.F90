!--------------------------------------------------------------------------------------------------------------------
module module_error
!--------------------------------------------------------------------------------------------------------------------
! This module is used whenever one of the following errors occur:
! A) Error while the program is openning a file. 
!    Possible reasons:
!	- The file is not located in that directory.
!	- The file is somehow "corrupted".
! B) Error due to a variable-size definition.
!
    interface

        subroutine openError(path)
            implicit none
            character(100), intent (in) :: path
        end subroutine openError

        subroutine sizeError(flagE, svar, smax)
            implicit none
            integer (kind=8), intent(in) 	:: svar, smax
            character*4, intent(in) :: flagE
        end subroutine sizeError

        subroutine molnameError(molN)
            implicit none
            integer (kind=8), intent(in)    :: molN
        end subroutine molnameError

        subroutine isoAconameERROR(mol,iso)
            implicit none
            integer (kind=8), intent(in) :: mol,iso
        end subroutine isoAconameERROR

        subroutine addMolError(flagE, var)
            implicit none
            integer (kind=8), intent(in) :: var
            character*4, intent(in) :: flagE
        end subroutine addMolError

        subroutine sumRuleERROR()
            implicit none
        end subroutine sumRuleERROR

        subroutine offdiagEL_IN_ERROR(v1, v2, v3, v4, v5, v6, v7)
            implicit none
            double precision, intent(in)      :: v1, v2, v3, v4, v5, v6, v7
        end subroutine offdiagEL_IN_ERROR

        subroutine offdiagEL_ERROR(var1, var2, var3, var4)
            implicit none
            double precision, intent(in):: var1, var2, var3, var4
        end subroutine offdiagEL_ERROR

    end interface


end module module_error
!--------------------------------------------------------------------------------------------------------------------
subroutine openError(path)
!--------------------------------------------------------------------------------------------------------------------
! It displays an error message whenever an error occurs while the program is opening a file.
!
implicit none
    character*100, intent (in) :: path

    print*, "Open file Error: " // trim(path) // " does not exit."
    print*, "-----> Stop program"
    stop
end subroutine openError
!--------------------------------------------------------------------------------------------------------------------
subroutine sizeError(flagE, svar, smax)
!--------------------------------------------------------------------------------------------------------------------
! It displays an error message whenever an error occurs while the program is saving
! records from a file into a variable.
!
    implicit none
    integer (kind=8), intent(in) 	:: svar, smax
    character*4, intent(in) 	:: flagE

    if (flagE .eq. "1000") then
      write (*,1000) svar, smax
    else if (flagE .eq. "1001") then
      write (*,1001) 
    endif
      
    STOP
    
1000    Format(1x,"****************** Fortran STOP ******************",&
        1x,"The Number of points to compute (",I6,")",&
        1x,"is TOO large regarding the fortran code implementation.",&
        1x,"Raise the value of -nLmx- to at least ",I6)
	     
1001    Format(1x,"****************** Fortran STOP ******************",&
        1x,'Arrays in for Line data storage are too small', &
        1x,'raise the value of nLmx in the "Common Var. module". ')

end subroutine sizeError
!--------------------------------------------------------------------------------------------------------------------
subroutine molnameError(molN)
!--------------------------------------------------------------------------------------------------------------------
! It displays an error message whenever an error occurs while matching
! molecule Name.
!
    implicit none
    integer (kind=8), intent(in)    :: molN

      write (*,1002) molN
      
    STOP
    
1002    Format(1x,"****************** molid_char: name not-found.",&
        1x,"Molecule (HITRANid):",I4,"has not a valid register",&
        1x,"since there is no matching on the program internal database",&
        1x,"please, check the ID no, or update the subroutine molid_char.")         

end subroutine molnameError
!--------------------------------------------------------------------------------------------------------------------
subroutine addMolError(flagE, var)
!--------------------------------------------------------------------------------------------------------------------
! It displays an error message whenever an error occurs while the program is saving
! extra information of the molecule selected.
!
    implicit none
    integer (kind=8), intent(in):: var
    character*4, intent(in)     :: flagE

    if (flagE .eq. "1003") then
      write (*,1003) var
    else if (flagE .eq. "1004") then
      write (*,1004) var
    endif
      
    STOP
    
1003    Format(1x,"****************** addMolParam: B0 to O2 isotopologue missing",&
        1x,"O2 isotope #",I6," is not available in Arts line-mixing program.",&
        1x,"please update the code if possible.")
         
1004    Format(1x,"****************** addMolParam: molecule missing",&
        1x,'Sorry, no specific information about your selected molecule:',I6, &
        1x,'Please update the requiered information in the code.')

end subroutine addMolError
!--------------------------------------------------------------------------------------------------------------------
subroutine isoAconameERROR(mol, iso)
!--------------------------------------------------------------------------------------------------------------------
! It displays an error message whenever an error occurs while matching
! molecule Name.
!
    implicit none
    integer (kind=8), intent(in)    :: mol, iso

      write (*,1005) iso, mol
      
    STOP
    
1005    Format(1x,"****************** molid_char: isotopologue AFGL code not-found.",&
        1x,"The isopologue ",I4," from the molecule (HITRANid):",I4,&
        1x,"has not a valid AFGL register.",&
        1x,"please, check the isopologue no, or update the subroutine molid_char.")         

end subroutine isoAconameERROR
!--------------------------------------------------------------------------------------------------------------------
subroutine sumRuleERROR()
!--------------------------------------------------------------------------------------------------------------------
! It displays an error message whenever an the Sum-Rule is not full-filled.
! 
    implicit none

    write (*,1006) 
      
    STOP
    
1006    Format(1x,"****************** sumRule: verification error.",&
        1x,"The calculation DOES NOT verifies the sum-rule.",&
        1x,"Uncoment dialogs in Subroutine:sumRule to check affected lines.")        

end subroutine sumRuleERROR
!--------------------------------------------------------------------------------------------------------------------
subroutine offdiagEL_IN_ERROR(v1, v2, v3, v4, v5, v6, v7)
!--------------------------------------------------------------------------------------------------------------------
! It displays an error message whenever an element of series generating off-diagonal elements of the relaxation 
! matrix overflows.
! 
    implicit none
    double precision, intent(in)      :: v1, v2, v3, v4, v5, v6, v7

      write (*, *) "****************** K_jkCalc/K_jkO2: variable overflow "
      write (*, *) "Either wigner symbol:"
      write (*, *) "w3j1=",v3
      write (*, *) "w3j2", v4
      write (*, *) "w6j1=",v5
      write (*, *) "w6j2=",v6
      write (*, *) "w6j3=",v7
      write (*, *) "2nd Adiabatic Factor: v1=",v1
      write (*, *) "Or the basis-rate function: Q(L)=",v2
      write (*, *) "overflows."
      write (*, *) "Typically this error occurs when an arithmetic operation attempts"
      write (*, *) "to create a numeric value that is too large to be represented "
      write (*, *) "within the available storage space. Check the inputs and/or code."
    STOP          

end subroutine offdiagEL_IN_ERROR
!--------------------------------------------------------------------------------------------------------------------
subroutine offdiagEL_ERROR(var1, var2, var3, var4)
!--------------------------------------------------------------------------------------------------------------------
! It displays an error message if one of the resulting off-diagonal element of the relaxation 
! matrix overflows.
! 
    implicit none
    double precision, intent(in):: var1, var2, var3, var4

      write (*,*) "****************** K_jkCalc/K_jkO2: variable overflow "
      write (*, *) "Any of the following variable:"
      write (*, *) "C1=[Ni][Ni'][Nf+][Nf'][Jf][Jf'][Ji']^2=",var1
      write (*, *) "C2=(-1)^(Ji'+Ji+n)=",var2
      write (*, *) "** Note that [Ni]=2Ni+1."
      write (*, *) "1st Adiabatic Factor: AF1=",var3
      write (*, *) "Or the sum which results in the off-diag. element: suma=",var4
      write (*, *) "overflows."
      write (*, *) "Typically this error occurs when an arithmetic operation attempts"
      write (*, *) "to create a numeric value that is too large to be represented "
      write (*, *) "within the available storage space. Check the inputs and/or code."
      
    STOP

end subroutine offdiagEL_ERROR
!--------------------------------------------------------------------------------------------------------------------