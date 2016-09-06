!--------------------------------------------------------------------------------------------------------------------
module module_OC_files
!--------------------------------------------------------------------------------------------------------------------
    interface
        
        subroutine openFile(u, name, error_open)
            use module_error
            implicit none
            character(100), intent(in) 	:: name
            integer (kind=8), intent(in) 	:: u
            logical, intent(out) 	:: error_open
        end subroutine openFile

        subroutine closeFile(u, error_open)
            implicit none
            integer (kind=8), intent(in) 	:: u
            logical, intent(in) 	:: error_open
        end subroutine closeFile

    end interface


end module module_OC_files
!--------------------------------------------------------------------------------------------------------------------
subroutine openFile(u, fname, error_open)
!--------------------------------------------------------------------------------------------------------------------
! This subroutine open files. 
! General function + error checking.
!
    use module_error
    implicit none
    character(100), intent(in) 	:: fname
    integer (kind=8), intent(in) 	:: u
    logical, intent(out) 	:: error_open
    integer (kind=8)        :: file_error

    open (unit = u, file = trim(fname), status = 'OLD', action = 'READ', IOSTAT = file_error)
    !print *, file_error
    !stop
    if (file_error == 0) then
        error_open = .false.
    else
        error_open = .true.
        !print*, "Unable to open: ", trim(fname)
        !stop
        call openError(fname)
    end if

end subroutine openFile
!--------------------------------------------------------------------------------------------------------------------
subroutine closeFile(u, error_open)
!--------------------------------------------------------------------------------------------------------------------
! This subroutine close files. 
! General function.
!
    implicit none
    integer (kind=8), intent(in) :: u
    logical, intent(in) :: error_open

    if (.not.error_open) close (u) ! Si no ha habido error al abrirlo, entonces se cierra.

end subroutine closeFile