!--------------------------------------------------------------------------------------------------------------------
module module_error
!--------------------------------------------------------------------------------------------------------------------
! This module is used whenever one needs debugging.
!
    interface

        subroutine memoError(chvar,econ)
            use module_common_var
            implicit none
            character(6), intent (in)      :: chvar
            type (dta_ERR), intent (inout)   :: econ
        end subroutine memoError

        subroutine openError(path, econ)
            use module_common_var
            implicit none
            character(100), intent (in)      :: path
            type (dta_ERR), intent (inout)   :: econ
        end subroutine openError

        subroutine errorSPIN(econ)
            use module_common_var
            implicit none
            type (dta_ERR), intent (inout)   :: econ
        end subroutine errorSPIN

        subroutine sizeError(flagE, svar, smax, econ)
            use module_common_var
            implicit none
            integer (kind=8), intent(in):: svar, smax
            character*4, intent(in)          :: flagE
            type (dta_ERR), intent (inout)   :: econ
        end subroutine sizeError

        subroutine molnameError(molN,econ)
            use module_common_var
            implicit none
            integer (kind=8), intent(in)     :: molN
            type (dta_ERR), intent (inout)   :: econ
        end subroutine molnameError

        subroutine isoAconameERROR(mol,iso,econ)
            use module_common_var
            implicit none
            integer (kind=8), intent(in)     :: mol,iso
            type (dta_ERR), intent (inout)   :: econ
        end subroutine isoAconameERROR

        subroutine addMolError(flagE, var,econ)
            use module_common_var
            implicit none
            integer (kind=8), intent(in):: var
            character(4),     intent(in)     :: flagE
            type (dta_ERR), intent (inout)   :: econ
        end subroutine addMolError

        subroutine Qparam_error(econ)
            use module_common_var
            implicit none
            type (dta_ERR), intent (inout)   :: econ
        end subroutine Qparam_error

        subroutine sumRuleERROR(econ)
            use module_common_var
            implicit none
            type (dta_ERR), intent (inout)   :: econ
        end subroutine sumRuleERROR

        subroutine wignerS_ERROR(v3, v4, v5, v6, v7, econ)
            use module_common_var
            implicit none
            double precision, intent(in)     :: v3, v4, v5, v6, v7
            type (dta_ERR), intent (inout)   :: econ
        end subroutine wignerS_ERROR

        subroutine offdiagEL_IN_ERROR(v0, v1, v2, econ)
            use module_common_var
            implicit none
            double precision, intent(in)     :: v0, v1, v2
            type (dta_ERR), intent (inout)   :: econ
        end subroutine offdiagEL_IN_ERROR

        subroutine offdiagEL_ERROR(var1, var2, var3, econ)
            use module_common_var
            implicit none
            double precision, intent(in)     :: var1, var2, var3
            type (dta_ERR), intent (inout)   :: econ
        end subroutine offdiagEL_ERROR

        subroutine renorm_error(flag, n, k, W, Su, econ)
            use module_common_var
            implicit none
            integer(kind=8) , intent (in   ) :: n, k
            double precision, intent (in   ) :: W, Su
            character(4)    , intent (in   ) :: flag
            type (dta_ERR)  , intent (inout) :: econ
        end subroutine renorm_error

        subroutine W_error(flag, n, k, W, econ)
            use module_common_var
            implicit none
            integer(kind=8) , intent (in   ) :: n, k
            double precision, intent (in   ) :: W
            character(4)    , intent (in   ) :: flag
            type (dta_ERR)  , intent (inout) :: econ
        end subroutine W_error

        subroutine LLS_error(b1, b2, b3, info, econ)
            use module_common_var
            implicit none
            integer*8       , intent (in   ) :: info
            double precision, intent (in   ) :: b1, b2, b3
            type (dta_ERR)  , intent (inout) :: econ
        end subroutine LLS_error

        subroutine errorBubble(econ)
            use module_common_var
            implicit none
            type (dta_ERR), intent (inout)   :: econ
        end subroutine errorBubble

        subroutine errorPType(pty,econ)
            use module_common_var
            implicit none
            character(3)  , intent (in   )   :: pty
            type (dta_ERR), intent (inout)   :: econ
        end subroutine errorPType

    end interface


end module module_error
!--------------------------------------------------------------------------------------------------------------------
subroutine memoError(chvar, econ)
!--------------------------------------------------------------------------------------------------------------------
! It displays an error message whenever an error occurs while the program is opening a file.
!
    use module_common_var
    implicit none
    character*6, intent (in) :: chvar
    type (dta_ERR), intent (inout)   :: econ

    if (econ % e(1) .ge. 1) then
      print*, "Not enough memory to allocate var:" // trim(chvar) // "."
    endif
    econ % e(2) = econ % e(2) + 1
end subroutine memoError
!--------------------------------------------------------------------------------------------------------------------
subroutine openError(path, econ)
!--------------------------------------------------------------------------------------------------------------------
! It displays an error message whenever an error occurs while the program is opening a file.
!
    use module_common_var
    implicit none
    character*100, intent (in) :: path
    type (dta_ERR), intent (inout)   :: econ

    if (econ % e(1) .ge. 1) then
      print*, "Open file Error: " // trim(path) // " does not exit."
    endif
    econ % e(2) = econ % e(2) + 1
end subroutine openError
!--------------------------------------------------------------------------------------------------------------------
subroutine errorSPIN(econ)
!--------------------------------------------------------------------------------------------------------------------
! It displays an error message whenever an error occurs while the program is opening a file.
!
    use module_common_var
    implicit none
    type (dta_ERR), intent (inout)   :: econ

    if (econ % e(1) .ge. 1) then
      print*, "Not SPIN-data available for O2"
    endif
    econ % e(2) = econ % e(2) + 1
end subroutine errorSPIN
!--------------------------------------------------------------------------------------------------------------------
subroutine sizeError(flagE, svar, smax, econ)
!--------------------------------------------------------------------------------------------------------------------
! It displays an error message whenever an error occurs while the program is saving
! records from a file into a variable.
!
    use module_common_var
    implicit none
    integer (kind=8), intent(in) 	   :: svar, smax
    character(4)    , intent(in)     :: flagE
    type (dta_ERR), intent (inout)   :: econ

    if (econ % e(1) .ge. 1) then
      if (flagE .eq. "1000") then
        write (*,1000) svar, smax
      else if (flagE .eq. "1001") then
        write (*,1001) 
      endif
    endif
    econ % e(2) = econ % e(2) + 1
    
1000    Format(1x,"****************** Hit2DTA: memory access violation.",&
        1x,"The Number of points to compute (",I6,")",&
        1x,"is TOO large regarding the fortran code implementation.",&
        1x,"Raise the value of -nLmx- to at least ",I6)
	     
1001    Format(1x,"****************** DipCal: memory access violation.",&
        1x,'Arrays in for Line data storage are too small', &
        1x,'raise the value of nLmx in the "module_common_var". ')

end subroutine sizeError
!--------------------------------------------------------------------------------------------------------------------
subroutine molnameError(molN, econ)
!--------------------------------------------------------------------------------------------------------------------
! It displays an error message whenever an error occurs while matching
! molecule Name.
!
    use module_common_var
    implicit none
    integer (kind=8), intent(in)    :: molN
    type (dta_ERR), intent (inout)   :: econ

    if (econ % e(1) .ge. 1) then
      write (*,1002) molN
    endif
    econ % e(2) = econ % e(2) + 1
    
1002    Format(1x,"****************** molid_char: name not-found.",&
        1x,"Molecule (HITRANid):",I4,"has not a valid register",&
        1x,"since there is no matching on the program internal database",&
        1x,"please, check the ID no, or update the subroutine molid_char.")         

end subroutine molnameError
!--------------------------------------------------------------------------------------------------------------------
subroutine addMolError(flagE, var, econ)
!--------------------------------------------------------------------------------------------------------------------
! It displays an error message whenever an error occurs while the program is saving
! extra information of the molecule selected.
!
    use module_common_var
    implicit none
    integer (kind=8), intent(in):: var
    character(4)    , intent(in)     :: flagE
    type (dta_ERR), intent (inout)   :: econ

    if (econ % e(1) .ge. 1) then
      if (flagE .eq. "1003") then
        write (*,1003) var
      else if (flagE .eq. "1004") then
        write (*,1004) var
      endif
    endif
    econ % e(2) = econ % e(2) + 1
    
1003    Format(1x,"****************** addMolParam: B0 to O2 isotopologue missing",&
        1x,"O2 isotope #",I6," is not available in Arts line-mixing program.",&
        1x,"please update the code if possible.")
         
1004    Format(1x,"****************** addMolParam: molecule missing",&
        1x,'Sorry, no specific information about your selected molecule:',I6, &
        1x,'Please update the requiered information in the code.')

end subroutine addMolError
!--------------------------------------------------------------------------------------------------------------------
subroutine isoAconameERROR(mol, iso, econ)
!--------------------------------------------------------------------------------------------------------------------
! It displays an error message whenever an error occurs while matching
! molecule Name.
!
    use module_common_var
    implicit none
    integer (kind=8), intent(in)     :: mol, iso
    type (dta_ERR), intent (inout)   :: econ

    if (econ % e(1) .ge. 1) then
      write (*,1005) iso, mol
    endif
    econ % e(2) = econ % e(2) + 1
    
1005    Format(1x,"****************** molid_char: isotopologue AFGL code not-found.",&
        1x,"The isopologue ",I4," from the molecule (HITRANid):",I4,&
        1x,"has not a valid AFGL register.",&
        1x,"please, check the isopologue no, or update the subroutine molid_char.")         

end subroutine isoAconameERROR
!--------------------------------------------------------------------------------------------------------------------
subroutine Qparam_error(econ)
!--------------------------------------------------------------------------------------------------------------------
! It displays an error message whenever there are not available parameters for the basis rate function.
! 
    use module_common_var
    implicit none
    type (dta_ERR), intent (inout)   :: econ

    if (econ % e(1) .ge. 1) then
      write (*,1006) 
    endif
    econ % e(2) = econ % e(2) + 1 
    
1006    Format(1x,"****************** systemQParam: missing a1, a2, a3",&
        1x,"No adjusted parameter for the Q basis rate avaible.",&
        1x,"Please, update the data base.")         

end subroutine Qparam_error
!--------------------------------------------------------------------------------------------------------------------
subroutine sumRuleERROR(econ)
!--------------------------------------------------------------------------------------------------------------------
! It displays an error message whenever an the Sum-Rule is not full-filled.
! 
    use module_common_var
    implicit none
    type (dta_ERR), intent (inout)   :: econ

    if (econ % e(1) .ge. 1) then
      write (*,1007) 
    endif
    econ % e(2) = econ % e(2) + 1 
    
1007    Format(1x,"****************** sumRule: verification error.",&
        1x,"The calculation DOES NOT verifies the sum-rule.",&
        1x,"Uncoment dialogs in Subroutine:sumRule to check affected lines.")        

end subroutine sumRuleERROR
!--------------------------------------------------------------------------------------------------------------------
subroutine wignerS_ERROR(v3, v4, v5, v6, v7, econ)
!--------------------------------------------------------------------------------------------------------------------
! It displays an error message whenever an element of series generating off-diagonal elements of the relaxation 
! matrix overflows.
! 
    use module_common_var
    implicit none
    double precision, intent(in)      :: v3, v4, v5, v6, v7
    type (dta_ERR), intent (inout)    :: econ

    if (econ % e(1) .ge. 1) then
      write (*, *) "****************** Kpart2/Kpart2_O2: variable overflow "
      write (*, *) "Either wigner symbol:"
      write (*, *) "w3j1=",v3
      write (*, *) "w3j2", v4
      write (*, *) "w6j1=",v5
      write (*, *) "w6j2=",v6
      write (*, *) "w6j3=",v7
      !write (*, *) "2nd Adiabatic Factor: v1=",v1
      !write (*, *) "Or the basis-rate function: Q(L)=",v2
      write (*, *) "overflows."
      write (*, *) "Typically this error occurs when an arithmetic operation attempts"
      write (*, *) "to create a numeric value that is too large to be represented "
      write (*, *) "within the available storage space. Check the inputs and/or code."
    endif
    econ % e(2) = econ % e(2) + 1          

end subroutine wignerS_ERROR
!--------------------------------------------------------------------------------------------------------------------
subroutine offdiagEL_IN_ERROR(v0, v1, v2, econ)
!--------------------------------------------------------------------------------------------------------------------
! It displays an error message whenever an element of series generating off-diagonal elements of the relaxation 
! matrix overflows.
! 
    use module_common_var
    implicit none
    double precision, intent(in)     :: v0, v1, v2
    type (dta_ERR), intent (inout)   :: econ

    if (econ % e(1) .ge. 1) then
      write (*, *) "****************** Kpart2/Kpart2_O2: variable overflow "
      write (*, *) "2nd Adiabatic Factor: v1=",v0
      write (*, *) "Or the basis-rate function: Q(L)=",v1
      write (*, *) "Inner-part of the infinite sum: suma(L)=",v2
      write (*, *) "overflows."
      write (*, *) "Typically this error occurs when an arithmetic operation attempts"
      write (*, *) "to create a numeric value that is too large to be represented "
      write (*, *) "within the available storage space. Check the inputs and/or code."
    endif
    econ % e(2) = econ % e(2) + 1          

end subroutine offdiagEL_IN_ERROR
!--------------------------------------------------------------------------------------------------------------------
subroutine offdiagEL_ERROR(var1, var2, var3, econ)
!--------------------------------------------------------------------------------------------------------------------
! It displays an error message if one of the resulting off-diagonal element of the relaxation 
! matrix overflows.
! 
    use module_common_var
    implicit none
    double precision, intent(in)     :: var1, var2, var3
    type (dta_ERR), intent (inout)   :: econ

    if (econ % e(1) .ge. 1) then
      write (*,*) "****************** K_jkCalc/K_jkO2: variable overflow "
      write (*, *) "Any of the following variable:"
      write (*, *) "C1=(-1)^(Ji'+Ji+n)[Ni][Ni'][Nf+][Nf']*"
      write (*, *) "[Jf][Jf'][Ji']^2 = ",var1
      write (*, *) "** Note that [Ni]=2Ni+1."
      write (*, *) "1st Adiabatic Factor: AF1=",var2
      write (*, *) "Or the sum which results in the off-diag. element: suma=",var3
      write (*, *) "overflows."
      write (*, *) "Typically this error occurs when an arithmetic operation attempts"
      write (*, *) "to create a numeric value that is too large to be represented "
      write (*, *) "within the available storage space. Check the inputs and/or code."
    endif      
    econ % e(2) = econ % e(2) + 1

end subroutine offdiagEL_ERROR
!--------------------------------------------------------------------------------------------------------------------
subroutine renorm_error(flag, n, k, W, Su, econ)
!--------------------------------------------------------------------------------------------------------------------
! It displays an error message whenever an error occurs during
! the relaxation matrix' renormalization.
!
    use module_common_var
    implicit none
    integer(kind=8) , intent (in   ) :: n, k
    double precision, intent (in   ) :: W, Su
    character(4)    , intent (in   ) :: flag
    type (dta_ERR)  , intent (inout) :: econ

    if (econ % e(1) .ge. 1) then
      if (flag .eq. "1008") then
        write (*,1008) Su,n,k,W
      else if (flag .eq. "1009") then
        write (*,1009) Su,n,k,W
      endif
    endif
    econ % e(2) = econ % e(2) + 1
    
1008    Format(1x,"****************** RN_Wmat: upper Sum-Renormatization process- overflows",&
        1x,"Upper sum: Sup=",e12.2 &
        1x,"Last elemet if the sum is: W_rn(",i3,",",i3,")=",e12.2)
         
1009    Format(1x,"****************** RN_Wmat: Lower Sum-Renormatization process- overflows",&
        1x,"Lower sum: Slow=",e12.2 &
        1x,"Last elemet if the sum is: W_rn(",i3,",",i3,")=",e12.2)

end subroutine renorm_error
!--------------------------------------------------------------------------------------------------------------------
subroutine W_error(flag, n, k, W, econ)
!--------------------------------------------------------------------------------------------------------------------
! It displays an error message whenever an error occurs during
! the relaxation matrix calculation.
!
    use module_common_var
    implicit none
    integer(kind=8) , intent (in   ) :: n, k
    double precision, intent (in   ) :: W
    character(4)    , intent (in   ) :: flag
    type (dta_ERR)  , intent (inout) :: econ

    if (econ % e(1) .ge. 1) then
      if (flag .eq. "1010") then
        write (*,1010) n,k,W
      else if (flag .eq. "1011") then
        write (*,1011) n,k,W
      endif
    endif
    econ % e(2) = econ % e(2) + 1
    
1010    Format(1x,"****************** WelCAL: Relaxation Matrix element overflow",&
        1x,"Either element W(n,k) or W(n,k) of the matrix overflowed." &
        1x,"Please check the inputs of the code, calculations will be not valid.", &
        1x,">> Downwards transition: W(",i3,",",i3,")=",e12.2)
         
1011    Format(1x,"****************** WelCAL: Relaxation Matrix element overflow",&
        1x,"Either element W(n,k) or W(n,k) of the matrix overflowed." &
        1x,"Please check the inputs of the code, calculations will be not valid.", &
        1x,">> Upwards transition: W(",i3,",",i3,")=",e12.2)

end subroutine W_error
!--------------------------------------------------------------------------------------------------------------------
subroutine LLS_error(b1, b2, b3, info, econ)
!--------------------------------------------------------------------------------------------------------------------
! It displays an error message if an error occured at
! the Linear-Least Square method (Lapack) and its Solution is not valid.
!
    use module_common_var
    implicit none
    integer*8       , intent (in   ) :: info
    double precision, intent (in   ) :: b1, b2, b3
    type (dta_ERR)  , intent (inout) :: econ

    if (econ % e(1) .ge. 1) then
        write (*,1012) 
        write (*,1013) info
        write (*,*), "< 0:  if INFO = -i, the i-th argument had an illegal value."
        write (*,*), "> 0:  if INFO =  i, the i-th diagonal element of the"
        write (*,*), "      triangular factor of A is zero, so that A does not have"
        write (*,*), "      full rank; the least squares solution could not be computed."
        write (*,1014) b1,b2,b3
    endif
    econ % e(2) = econ % e(2) + 1
    
1012    Format(1x,"****************** calc_QParam: Lapack internal error",&
        1x,"The Least squares solution could not be computed." &
        1x,"Please check the input matrix element in the code, calculations are not valid.")

1013    Format(1x,"Lapack internal info:",i8)

1014    Format(1x,">> Solution: A1=",e12.2,", A2=",e12.2,", A3=",e12.2)

end subroutine LLS_error
!--------------------------------------------------------------------------------------------------------------------
subroutine errorBranch(pos, econ)
!--------------------------------------------------------------------------------------------------------------------
! It displays an error message whenever any transition does not follow 
! ∆J=0, +-1, +-2
!------------------------------
! NOTATION: ∆J = Jupper-Jlower
!------------------------------
! 
    use module_common_var
    implicit none
    double precision, intent(in)      :: pos
    type (dta_ERR), intent (inout)    :: econ

    if (econ % e(1) .ge. 1) then
        write(*,*), "****************** delta2branch/branch2delta:"
        write(*,*), "Transition in position", pos,"does not follows selection rules."
        write(*,*), "NOTE: HITRAN is an empirical DB it should not contain any line"
        write(*,*), "that does not follow the selection rules:"
        write(*,*), "∆J=0, +-1, +-2"
    endif
    econ % e(2) = econ % e(2) + 1          

end subroutine errorBranch
!--------------------------------------------------------------------------------------------------------------------
subroutine errorBubble(econ)
!--------------------------------------------------------------------------------------------------------------------
! It displays an error message whenever an error occurs while the program is opening a file.
!
    use module_common_var
    implicit none
    type (dta_ERR), intent (inout)   :: econ

    if (econ % e(1) .ge. 1) then
        write(*,*), "****************** bubble_index:"
        write(*,*), 'Not supported kind of the order option, use (a) or (d) instead.'
    endif
    econ % e(2) = econ % e(2) + 1
end subroutine errorBubble
!--------------------------------------------------------------------------------------------------------------------
subroutine errorPType(pty,econ)
!--------------------------------------------------------------------------------------------------------------------
    use module_common_var
    implicit none
    character(3)  , intent (in   )   :: pty
    type (dta_ERR), intent (inout)   :: econ

    if (econ % e(1) .ge. 1) then
        write(*,*), "****************** PopuCAL:"
        write(*,*), "Non-Valid Population calculation type:", pty
    endif
    econ % e(2) = econ % e(2) + 1
end subroutine errorPType
!--------------------------------------------------------------------------------------------------------------------