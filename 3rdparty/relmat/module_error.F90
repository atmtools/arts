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
    
1000    Format(1x,"************* PROBLEM !!!! ******************",&
        1x,"The Number of points to compute (",I6,")",&
        1x,"is TOO large for Arrays ---> Program Stop",&
        1x,"Raise the value of -nSigmx- to at least ",I6,&
        "in ALL Parameter Statements where it appears. ")
	     
1001    Format(1x,'************ PROBLEM !!!! ******************', &
        1x,'Arrays in for Line data storage are too small', &
        1x,'raise the value of nLmx in the "Common Var. module". ')

end subroutine sizeError
!--------------------------------------------------------------------------------------------------------------------