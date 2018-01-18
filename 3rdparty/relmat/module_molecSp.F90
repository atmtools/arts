MODULE module_molecSp
!--------------------------------------------------------------------------------------------------------------------
! This module contains all subroutines related to molecular symmetry and behavior
!
    interface

        subroutine r_arts_LocalQ(dta1,pos,Q0,Q00,econ)
            use module_common_var
            use module_error
            use module_maths
            implicit none
            integer*8              , intent(inout) :: pos
            integer*8              , intent(in   ) :: Q0(4), Q00(4)
            type (dta_SDF), pointer, intent(inout) :: dta1
            type (dta_ERR)         , intent(inout) :: econ
        end subroutine r_arts_LocalQ

        character*1 function delta2branch(delta,pos,econ)
            use module_common_var
            use module_error
            implicit none
            integer (kind=8), intent(in)    :: delta, pos
            type (dta_ERR)  , intent(inout) :: econ
        end function delta2branch

        integer*8 function branch2delta(branch,pos,econ)
            use module_common_var
            use module_error
            implicit none
            character(1)    , intent(in)    :: branch
            integer (kind=8), intent(in)    :: pos
            type (dta_ERR)  , intent(inout) :: econ
        end function branch2delta

        subroutine systemQParam(sys,molP,econ)
            use module_common_var
            use module_error
            implicit none
            integer*8, intent(in)         :: sys(2)
            type (dta_MOL), intent(inout) :: molP
            type (dta_ERR), intent(inout) :: econ
        end subroutine systemQParam

        subroutine systemQParam_LLS(sys,molP)
            use module_common_var
            implicit none
            integer*8, intent(in)         :: sys(2)
            type (dta_MOL), intent(inout) :: molP
        end subroutine systemQParam_LLS

        subroutine Ai_zeros(molP)
            use module_common_var
            implicit none
            type (dta_MOL), intent(inout) :: molP
        end subroutine Ai_zeros

    end interface

END module module_molecSp
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE r_arts_LocalQ(dta1,pos,Q0,Q00,econ) 
!--------------------------------------------------------------------------------------------------------------------
! r_arts_LocalQ : identifies localQ from ARTS' input.
!
! Detailed info:
! --------------
! This subroutine indentifies the upper and lower state local quanta for linear molecules 
! (definitions would be found in HITRAN04).
!
! Last Update: Teresa Mendaza 18/01/2018
!--------------------------------------------------------------------------------------------------------------------
        use module_common_var
        use module_error
        use module_maths
        implicit none
        integer*8              , intent(inout) :: pos 
        integer*8              , intent(in)    :: Q0(4), Q00(4) ! == artsUpp, artsLow
        type (dta_SDF), pointer, intent(inout) :: dta1
        type (dta_ERR)         , intent(inout) :: econ
        character( 1) :: delta2branch
        integer*8     :: branch2delta, iaux
        integer*8     :: my_mol, delta
        !--------------------------------------------------
        my_mol = dta1%M
        iaux = 0
        !  NOTES: 
        !  (1) Prime and double primes refer, respectively, 
        !  to upper and lower states, respectively, i.e.
        !       * upper := 2; 0 ; or ' 
        !       * lower := 1;00 ; or " 
        !--------------------------------------------------
        IF ((my_mol .eq. 2).or.(my_mol .eq. 4).or.(my_mol .eq. 5).or.(my_mol .eq. 7).or. &
            (my_mol .eq. 8).or.(my_mol .eq.13).or.(my_mol .eq.14).or.(my_mol .eq.15).or. &
            (my_mol .eq.16).or.(my_mol .eq.17).or.(my_mol .eq.18).or.(my_mol .eq.19).or. &
            (my_mol .eq.22).or.(my_mol .eq.23)) then
        !
        ! Diatomic: 
        ! CO(5), HF(14), HCl(15), HBr(16), HI(17), N2(22), NO(8), OH(13)
        ! ClO(18), O2(7)
        !
        ! Linear Triatomic: 
        ! N2O(4), OCS(19), HCN(23), CO2(2)
        ! --------------------------------------------------------------
        ! Variables in use:
        !            J: Resultant total angular momentum quantum number, 
        !               excluding nuclear spins.(initial/lower level)
        !            N: Rotational angular momentum quantum number
        !               (excluding electron and nuclear spins, 
        !               in the case where electron spin is present).
        !               J = N +- 1/2.
        !            S: Electronic spin (espin). 
        !           l2: Vibrational angular momentum.
        !
        !      Q0:                   Q00:                
        !               l2 J N S          l2 J N S     
        !
        ! LOWER LEVEL:
        ! -----------
        ! J:
            if (Q00(2) .eq. -1) then
                dta1%J(pos,1)     = 0.0d0
                iaux = -1
            else
                dta1%J(pos,1)     = real(Q00(2),dp)
            endif
        !
        ! N:
            if (Q00(3) .eq. -1) then
                dta1%N(pos,1)     = dta1%J(pos,1)
            else
                dta1%N(pos,1)     = Q00(3)
            endif
        !
        ! S:
            if (Q00(4) .eq. -1) then
                dta1%espin(pos,1) = 1.0d0
                if (my_mol .eq. 7) then
                !IF O2 and no SPIN then ERROR!!!!
                    call errorSPIN(econ)
                endif
            else
                dta1%espin(pos,1) = real(Q00(4),dp)
            endif
        !
        ! UPPER LEVEL:
        !
        ! J:
            if (Q0(2) .eq. -1) then
                dta1%J(pos,2)     = 0.0d0
                iaux = -1
            else
                dta1%J(pos,2)     = real(Q0(2),dp)
            endif
        !
        ! N:
            if (Q0(3) .eq. -1) then
                dta1%N(pos,2)     = dta1%J(pos,2)
            else
                dta1%N(pos,2)     = Q0(3)
            endif
        !
        ! S:
            if (Q0(4) .eq. -1) then
                dta1%espin(pos,2) = 1.0d0
                if (my_mol .eq. 7) then
                !! IF O2 and no SPIN then ERROR!!!!
                    call errorSPIN(econ)
                endif
            else
                dta1%espin(pos,2) = real(Q0(4),dp)
            endif
        !
        ! BRANCH:
        !
            delta = int(dta1%J(pos,2)-dta1%J(pos,1))
            dta1%br(pos) = delta2branch(delta,pos,econ)
        !
        pos = pos + iaux
        ELSE
            if (econ % e(1) .ge. 1) then
                write(*,*) "****************** r_arts_LocalQ:"
                write(*,*) " No vibrational band information or not speficied Format "
                write(*,*) " of your selected molecule (HITRANid):", my_mol  
            endif
            econ % e(2) = econ % e(2) + 1
        ENDIF

  END SUBROUTINE r_arts_LocalQ
!--------------------------------------------------------------------------------------------------------------------
  character*1 function delta2branch(delta, pos,econ)
!--------------------------------------------------------------------------------------------------------------------
! Delta2branch : Delta to branch
!
! Detailed info:
! --------------
! This function gives you the part of the branch to which transition in position 
! "pos" belongs given the difference between the upper J and the lower J.
! Possible branches are: O-, P-, Q-, R-, or S-.
!
! Author: Teresa Mendaza 02/03/2016
!--------------------------------------------------------------------------------------------------------------------
        use module_common_var
        use module_error
        implicit none
        integer*8     , intent(in)    :: delta, pos
        type (dta_ERR), intent(inout) :: econ
        !------------------------------------------------
        if ( delta .eq. 1 ) then
            ! Ju-J == +1
                delta2branch = "r"
        else if (delta .gt. 1 ) then
            ! Ju-J == 2
                delta2branch = "s"
                ! Does this make sense?
                !
        else if (delta .eq. -1 ) then
            ! Ju-J == -1
                delta2branch = "p"
        else if (delta .lt. -1 ) then
            ! Ju-J == -2
                delta2branch = "o"
                ! Does this make sense?
                !
        else if (delta .eq. 0) then 
            ! Ju-J == 0
                delta2branch = "q"
        else
            ! since HITRAN is an empirical DB it should not contain any line
            ! that does not follow the selection rules:
            ! \delta(J)=0, +-1, +-2
            ! \delta(l)= +-1
            ! However... the code checks
            call errorBranch(pos, econ)
        endif
  end function delta2branch
!--------------------------------------------------------------------------------------------------------------------
  integer*8 function branch2delta(branch,pos,econ)
!--------------------------------------------------------------------------------------------------------------------
! branch2delta : branch to delta
!
! Detailed info:
! --------------
! This function gives you the differenc in J given
! the part of the branch to which transition in position 
! "pos" belongs.
! Possible branches are: O-, P-, Q-, R-, or S-. 
! that give the following delta (respectively):
! -2, -1, 0, +1, +2
!
! NOTE! User has to watch out SELECTION RULES 
!       (for example: not every molecule has O- or S- branches)
!
! Author: Teresa Mendaza 02/03/2016
!--------------------------------------------------------------------------------------------------------------------
        use module_common_var
        use module_error
        implicit none
        character(1), intent(in)              :: branch
        integer (kind=8), intent(in)          :: pos
        type (dta_ERR)  , intent(inout)       :: econ
        !---------------------------------------------------
        ! NOTATION:
        ! ∆J = Jupper-Jlower
        !---------------------------------------------------
        if ( (branch .eq. "r") .or. (branch .eq. "R") ) then 
            ! *r = R-branch (when ∆J = +1) 
            branch2delta = 1
        else if ( branch .eq. "s" .or. (branch .eq. "S") ) then
            ! *s = S-branch (when ∆J = +2) 
            branch2delta = 2  
                ! Does this make sense?
                !
        else if ( branch .eq. "p" .or. (branch .eq. "P")) then
            ! *p = P-branch (when ∆J = -1).
            branch2delta = -1
        else if ( branch .eq. "o" .or. (branch .eq. "O") ) then
            ! *o = O-branch (when ∆J = -2).
            branch2delta = -2
        else if ( branch .eq. "q" .or. (branch .eq. "Q") .or. (branch .eq. "") ) then 
            ! *q = Q-branch (when ∆J =  0)
            branch2delta = 0
        else
            call errorBranch(pos, econ)
        endif
  end function branch2delta
!--------------------------------------------------------------------------------------------------------------------
  subroutine systemQParam(sys,molP,econ)
!--------------------------------------------------------------------------------------------------------------------
! systemQParam : Add the system (molecule-perturber) adjusted parameters 
!
! Detailed info:
! --------------
! This subroutine acts as a "look-up table" of the adjusted a1, a2, a3, dc
! parameters that are requeried for the Basis Rates (function 'Q_Mol_X').
! 
! Parameters written here come either from previous literature in the field or
! have been calculated.
!
! Author: Teresa Mendaza 17/02/2017
!--------------------------------------------------------------------------------------------------------------------
            use module_common_var
            use module_error
            implicit none
            integer*8, intent(in)         :: sys(2)
            type (dta_MOL), intent(inout) :: molP
            type (dta_ERR), intent(inout) :: econ
            double precision              :: T
            !-----------------------------------------
            T = molP % Temp
            ! NOTE: T0 needs to be considered as the 
            ! reference temperature to use this parameters
            ! T0 = 296 K
            !-----------------------------------------
            if ((sys(1) .eq. 2) .and. (sys(2) .eq. 22)) then
            ! CO2(2)-N2(4)
            ! -------------------------------
            ! Rodriguez et al. 1999; CO2 - N2
                molP%a1 = 0.0181D0 !cm-1/atm
                molP%a2 = 0.81D0   !cm-1
                molP%a3 = 0.008D0  !cm-1
                molP%dc = 2.2      !Å (amstrong)
                molP%ex1 = 0.85   
                molP%ex2 = 0.0152
            
                if (T .ne. 296d0) then
                    molP%a1 = molP%a1*(296d0/T)**molP%ex1 
                    molP%a2 = molP%a2*(296d0/T)**molP%ex2
                endif
            ! -------------------------------
            else if ((sys(1) .eq. 2) .and. (sys(2) .eq. 7)) then
            ! CO2(2)-O2(7)
            ! -------------------------------
            ! Rodriguez et al. 1999; CO2 - O2
                molP%a1 = 0.0168D0 !cm-1/atm
                molP%a2 = 0.82D0   !cm-1
                molP%a3 = 0.007D0  !cm-1
                molP%dc = 2.4      !Å (amstrong)
                molP%ex1 = 0.50   
                molP%ex2 = -0.091
                if (T .ne. 296d0) then
                    molP%a1 = molP%a1*(296d0/T)**molP%ex1 
                    molP%a2 = molP%a2*(296d0/T)**molP%ex2
                endif
            ! -------------------------------
            else if ((sys(1) .eq. 7) .and. (sys(2) .eq. 7)) then
            ! Tran et al. 2006; O2 - O2
                molP%a1 = 0.0275D0 !cm-1/atm
                molP%a2 = 0.935D0  !cm-1
                molP%a3 = 1.01D0   !cm-1
                molP%dc = 1.05     !Å (amstrong)
                molP%ex1 = 1.0  
                molP%ex2 = 1.0
                if (T .ne. 296d0) then
                    molP%a1 = molP%a1*(296d0/T)**molP%ex1 
                    molP%a2 = molP%a2*(296d0/T)**molP%ex2
                endif
            else if ((sys(1) .eq. 7) .and. (sys(2) .eq. 22)) then
            ! Tran et al. 2006; O2 - N2
                molP%a1 = 0.0285D0 !cm-1/atm
                molP%a2 = 0.950D0  !cm-1
                molP%a3 = 1.03D0   !cm-1
                molP%dc = 1.0      !Å (amstrong)
                molP%ex1 = 1.0  
                molP%ex2 = 1.0
                if (T .ne. 296d0) then
                    molP%a1 = molP%a1*(296d0/T)**molP%ex1 
                    molP%a2 = molP%a2*(296d0/T)**molP%ex2
                endif
            else if ((sys(1) .eq. 4) .and. (sys(2) .eq. 7)) then
            ! Hartmann et al. 1999; N2O - O2
                molP%a1 = 0.0147D0 !cm-1/atm
                molP%a2 = 0.77D0   !cm-1
                molP%a3 = 0.025D0  !cm-1
                molP%dc = 2.9      !Å (amstrong)
                molP%ex1 = 0.85   
                molP%ex2 = 1.0
                if (T .ne. 296d0) then
                    molP%a1 = molP%a1*(296d0/T)**molP%ex1 
                    molP%a2 = molP%a2*(296d0/T)**molP%ex2
                endif
            else if ((sys(1) .eq. 4) .and. (sys(2) .eq. 22)) then
            ! Hartmann et al. 1999; N2O - N2
                molP%a1 = 0.0174D0 !cm-1/atm
                molP%a2 = 0.77D0   !cm-1
                molP%a3 = 0.025D0  !cm-1
                molP%dc = 2.9      !Å (amstrong)
                molP%ex1 = 0.85   
                molP%ex2 = 1.0
                if (T .ne. 296d0) then
                    molP%a1 = molP%a1*(296d0/T)**molP%ex1 
                    molP%a2 = molP%a2*(296d0/T)**molP%ex2
                endif
            else
                call Qparam_error(econ)
            endif
    end subroutine systemQParam
!--------------------------------------------------------------------------------------------------------------------
  subroutine systemQParam_LLS(sys,molP)
!--------------------------------------------------------------------------------------------------------------------
! systemQParam : Add the system (molecule-perturber) adjusted parameters 
!
! Detailed info:
! --------------
! This subroutine acts as a "dataBase" of the adjusted a1, a2, a3, dc
! parameters that are requeried for the Basis Rates (function 'Q_Mol_X').
! 
! Parameters written here come either from previous literature in the field or
! have been calculated.
!
! Author: Teresa Mendaza 16/02/2017
!--------------------------------------------------------------------------------------------------------------------
            use module_common_var
            implicit none
            integer*8, intent(in)              :: sys(2)
            type (dta_MOL), intent(inout)      :: molP
            double precision                   :: T
            !-----------------------------------------
            T = molP % Temp
            !-----------------------------------------
            if ((sys(1) .eq. 2) .and. (sys(2) .eq. 22)) then
            ! CO2(2)-N2(4)
            ! -------------------------------
            ! Rodriguez et al. 1999; CO2 - N2
                molP%dc = 2.2      !Å (amstrong)
            ! -------------------------------
            else if ((sys(1) .eq. 2) .and. (sys(2) .eq. 7)) then
            ! CO2(2)-O2(7)
            ! -------------------------------
            ! Rodriguez et al. 1999; CO2 - O2
                molP%dc = 2.4      !Å (amstrong)
            ! -------------------------------
            else if ((sys(1) .eq. 7) .and. (sys(2) .eq. 7)) then
            ! Tran et al. 2006; O2 - O2
                molP%dc = 1.01     !Å (amstrong)
            else if ((sys(1) .eq. 7) .and. (sys(2) .eq. 22)) then
            ! Tran et al. 2006; O2 - N2
                molP%dc = 1.0      !Å (amstrong)
            else if ((sys(1) .eq. 4) .and. (sys(2) .eq. 7)) then
            ! Hartmann et al. 1999; N2O - O2
                molP%dc = 2.9      !Å (amstrong)
            else if ((sys(1) .eq. 4) .and. (sys(2) .eq. 22)) then
            ! Hartmann et al. 1999; N2O - N2
                molP%dc = 2.9      !Å (amstrong)
            else
                molP%availableParam = .false.
                molP%AF_ON = .false.
            endif
            call Ai_zeros(molP)

    end subroutine systemQParam_LLS
!--------------------------------------------------------------------------------------------------------------------
  subroutine Ai_zeros(molP)
!--------------------------------------------------------------------------------------------------------------------
! Ai_zeros : Are Ai parameters Zero? 
!
! Detailed info:
! --------------
! This subroutine checks whether the adjusted a1, a2, a3
! parameters are stored in the internal database.
!
! Author: Teresa Mendaza 16/02/2017
!--------------------------------------------------------------------------------------------------------------------
            use module_common_var
            implicit none
            type (dta_MOL), intent(inout)      :: molP
            !-----------------------------------------
            if ((molP%a1 .eq. 0.0_dp) .and. &
                (molP%a2 .eq. 0.0_dp) .and. &
                (molP%a1 .eq. 0.0_dp) ) then
                molP%availableParam = .false.
            endif

    end subroutine Ai_zeros
!--------------------------------------------------------------------------------------------------------------------
