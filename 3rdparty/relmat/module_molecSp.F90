MODULE module_molecSp
!--------------------------------------------------------------------------------------------------------------------
! This module contains all subroutines related to molecular symmetry and behavior
!
    interface

        subroutine moleculeID(my_mol,isotope,molP)
            use module_common_var
            use module_OC_files
            implicit none
            integer*8, intent (in)       :: my_mol , isotope
            type (dta_MOL), intent(inout)      :: molP
        end subroutine moleculeID

        subroutine addMolParam(molP)
            use module_common_var
            use module_OC_files
            implicit none
            type (dta_MOL), intent(inout)      :: molP
        end subroutine addMolParam

        subroutine molid_PF(molP,iso)
            use module_common_var
            use module_OC_files
            implicit none
            integer*8, intent(in)        :: iso
            type (dta_MOL), intent(inout)      :: molP
        end subroutine molid_PF

        subroutine r_arts_LocalQ(dta1,pos,Q0,Q00)
            use module_common_var
            use module_maths
            implicit none
            integer*8, intent(in)        :: pos
            integer*8, intent(in)        :: Q0(4), Q00(4)
            type (dta_SDF), intent(inout)      :: dta1
        end subroutine r_arts_LocalQ

        character function delta2branch(delta,pos)
            implicit none
            integer (kind=8), intent(in)       :: delta, pos
        end function delta2branch

        integer*8 function branch2delta(branch,pos)
            implicit none
            character, intent(in)              :: branch
            integer (kind=8), intent(in)       :: pos
        end function branch2delta

        subroutine systemQParam(sys,molP)
            use module_common_var
            implicit none
            integer*8, intent(in)         :: sys(2)
            type (dta_MOL), intent(inout)       :: molP
        end subroutine systemQParam

    end interface

END module module_molecSp
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE moleculeID(my_mol,isotope,molP)
!--------------------------------------------------------------------------------------------------------------------
! "moleculeID": Molecule HITRAN's ID
! 
! Detailed Description:
! ---------------------
! This subroutine check if the molecule specified in "module_common_var" by it's sc.formulation
! belongs to HITRAN records (user has to check whether the list of molecules is updated or not).
!
! Variables:
!
! Input/Output Parameters of Routine (Arguments or Common)
! ----------------------------------
! molecule: string containing the molecule's formulation.
! my_mol  : HITRAN's ID number for that molecule.
!
! Accessed Files:  'HITRAN_molList.txt'
! --------------
!
! Called Routines: None
! ---------------  
!
! Called By: Main Program
! ---------
!
!
! T.Mendaza, last change 15 February 2016
!--------------------------------------------------------------------------------------------------------------------
!
    use module_common_var
    use module_OC_files
    implicit none
    integer*8, intent (in)  :: my_mol, isotope
    type (dta_MOL), intent(inout) :: molP
    integer*8               :: i,j,k, mP_size
    integer*8               :: g_j
    double precision              :: Qmo
    character(  6)                :: auxmol
    character(100)                :: fname
    logical                       :: molfound
    integer*8               :: error_read
    logical                       :: error_open

!----------
! 
    i = 1
    auxmol = ""
    molfound = .false.
    fname = trim(in_file_path)//trim(in_file_molp)
    call openFile(u5, fname, error_open)
    do
      read (u5, 1001, iostat = error_read, err = 100, end = 101), auxmol, i 
      !stop
        if( i .eq. my_mol ) then
            !print*, auxmol,i
            molP % chmol = auxmol
            molP % iso_m = isotope
            molP % M = my_mol 
            molfound = .true. 
            do
                read (u5, *, iostat = error_read, err = 101, end = 101), &
                molP % Aco(mP_size), molP % IAb(mP_size), Qmo, &
                g_j, molP % mms(mP_size) 
                !print*, my_mol, "iso#",molP%Aco(mp_size), molP % IAb(mP_size)!, molP % mms(mP_size)
                mP_size = mP_size + 1
            enddo 
            !---------
            ! Completing molecular parameters.
            !
            ! print*, 'Adding miscelaneous molecule information...'
            CALL addMolParam(molP)
            exit      
        endif    
        
      100 if (error_read .ne. 0) cycle ! Reading Error => skip to next available data
      101 if (error_read .ne. 0) exit  ! End of file   => stops "do-loop"
    
    end do
    molP%m_size = mP_size - 1
    call closeFile(u5, error_open)
    if (molfound) then
      print*, "Molecule:",molP % chmol,"has HITRAN's ID:",my_mol
    else
      print*, "Molecule:",molP % chmol,"does not belong to HITRAN DB."
      stop
    endif 
!
1001  Format(A6,2x,i2)
! 
!   Data from HITRAN molparam.txt file
!   ----------------------------------
!   Dry air mixing ratio from
!   GLOBALVIEW-CH4: Cooperative Atmospheric Data Integration
!   Project - Methane. CD-ROM, NOAA ESRL, Boulder, Colorado
!   [Also available on Internet via anonymous FTP to ftp.cmdl.noaa.gov,
!   Path: ccg/ch4/GLOBALVIEW], 2009.
    Return
  END SUBROUTINE moleculeID
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE addMolParam(molP)
!--------------------------------------------------------------------------------------------------------------------
! "addMolParam": Add molecular paramters
! 
! Detailed Description:
! ---------------------
! This subroutine add miscelaneous molecules' information required by this software. 
! If there is no default information about the user's molecule in study, the later will have to 
! added manualy. 
!
! Variables:
!
! Input/Output Parameters of Routine (Arguments or Common)
! ----------------------------------
! molP    : dta type containing molecule's information.
!
! Accessed Files:  None
! --------------
!
! Called Routines: None
! ---------------  
!
! Called By: Main Program
! ---------
!
!
! T.Mendaza, last change 15 February 2016
!--------------------------------------------------------------------------------------------------------------------
!
    use module_common_var
    implicit none
    type (dta_MOL), intent(inout) :: molP
    
!----------
! 
    call molid_PF(molP, molP % Aco(molP%iso_m))
    IF (molP % M .eq. 1) then
    ! water, H2O
        print*, "no water available"
        stop
    ELSEIF (molP % M .eq. 2) then
    ! Carbon dioxide, CO2
        molP % Nmcon = 0.42896E-3! CO2 molar concentration at 296K (mol·cm3), 1MPa
        molP % B0    = 0.39021   ! CO2 Rotational constant B0 (cm-1) 
    ELSEIF (molP % M .eq. 6) then
    ! Methane, CH4
        molP % Nmcon = 0.41245E-4! CH4 molar concentration at 296K (mol·cm3), 1MPa
        molP % B0    = 5.2       ! CH4 Rotational constant B0 (cm-1) 
                                 ! [Brown et al. 2003] 
    ELSEIF (molP % M .eq. 7) then
    ! Oxygen, O2
        molP % Nmcon = 0.41245E-4! O2 molar concentration at 296K (mol·cm3), 1MPa
        !Rotational constant dependent on the vibrational state and isotope.
        if (molP % Aco(molP%iso_m) .eq. 66) then
            ! AFGL code = 66 => 16O2 (isotope 1 in HITRAN)
            ! If nu1 = 0 : B0 = 43100.430 MHz 
            !                 = 1.437 cm-1 (MHz*10^6/c)  
            ! If nu1 = 1 : B0 = 42626.398 MHz 
            !                 = 1.421 cm-1 
            molP % B0    = 1.43  ! O2 Rotational constant B0 (cm-1) 
                                 ! NIST 
        elseif (molP % Aco(molP%iso_m) .eq. 67) then
            ! AFGL code = 67 => 16O17O (isotope 3 in HITRAN)
            ! If nu1 = 0 : B0 = 40561.35 MHz 
            !                 = 1.353 cm-1 
            molP % B0    = 1.35  ! O2 Rotational constant B0 (cm-1) 
                                 ! NIST 
        elseif (molP % Aco(molP%iso_m) .eq. 68) then
            ! AFGL code = 67 => 16O18O (isotope 2 in HITRAN)
            ! If nu1 = 0 : B0 = 38313.761 MHz 
            !                 = 1.278 cm-1 
            ! If nu1 = 1 : B0 = 37916.618 MHz 
            !                 = 1.265 cm-1 
            molP % B0    = 1.27  ! O2 Rotational constant B0 (cm-1) 
                                 ! NIST 
        else 
            print*,"O2 isotope #", molP%iso_m, " not available in Arts data base."
        endif
        
    ELSE 
        print*, "addMolParam: sorry, no specific information about your selected molecule"
        stop
    ENDIF
! 
      Return
  END SUBROUTINE addMolParam
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE molid_PF(molP,iso)
!--------------------------------------------------------------------------------------------------------------------
! "molid_PF": Partition function coefficient per molecule ID
! 
! Detailed Description:
! ---------------------
! This subroutine check if the molecule-isotope specified has associatted 
! Partition function.
!
! Variables:
!
! Input/Output Parameters of Routine (Arguments or Common)
! ----------------------------------
! molP  : Molecule's basic information.
!
! Accessed Files:  'PF_coeff.txt'
! --------------
!
! Called Routines: None
! ---------------  
!
! Called By: addMolParam
! ---------
!
!
! T.Mendaza, last change 07 March 2016
!--------------------------------------------------------------------------------------------------------------------
!
    use module_common_var
    use module_OC_files
    implicit none
    integer*8, intent (in )  :: iso !== Molecules isotope(AFGL code)
    type (dta_MOL), intent(inout) :: molP
    integer*8        :: i,j,k
    integer*8        :: aux_iso
    double precision              :: Qc(4)
    double precision              :: T
    character(  6)                :: auxmol
    character(100)                :: fname
    logical                       :: molfound
    integer*8        :: error_read
    logical                       :: error_open

!-----------------------------------------
    T = molP % Temp
!-----------------------------------------
! 
    i = 1
    auxmol = ""
    molfound = .false.
    fname = trim(in_file_path)//trim(in_file_PFco)
    call openFile(u7, fname, error_open)
    do
      read (u7, 1001, iostat = error_read, err = 100, end = 101), auxmol, i 
      !stop
        if( i .eq. molP%M ) then
            !print*, auxmol,i
            do
                if ( (T .gt. 69.0d0) .and. (T .lt. 401.0d0) ) then
                    read (u7, *, iostat = error_read, err = 101, end = 101), &
                    aux_iso, Qc(1), Qc(2), Qc(3), Qc(4) 
                    read (u7, *)
                elseif ( (T .gt. 400.0d0) .and. (T .lt. 2005.0d0) ) then
                        read (u7, *)
                        read (u7, *, iostat = error_read, err = 101, end = 101), &
                        aux_iso, Qc(1), Qc(2), Qc(3), Qc(4)
                else
                    print*, "molid_PF: Temperature out of range."
                    stop
                endif
                if( aux_iso .eq. iso ) then
                    molfound = .true. 
                    do j=1,4
                       molP % Qcoef(j) = Qc(j)
                    enddo
                endif 
            enddo 
            exit      
        endif  
      100 if (error_read .ne. 0) cycle ! Reading Error => skip to next available data
      101 if (error_read .ne. 0) exit  ! End of file   => stops "do-loop"
    
    end do
    call closeFile(u7, error_open)
    if (.not.(molfound)) then
      print*, "Molecule (HITRANid):",molP%M,"is not included in this PF_coeff file."
      stop
    endif 
!
1001  Format(A6,2x,i2)
! 
!   HITRAN Data PF_coeff.txt file
!   ----------------------------------
    Return
  END SUBROUTINE molid_PF
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE r_arts_LocalQ(dta1,pos,Q0,Q00) 
!--------------------------------------------------------------------------------------------------------------------
        use module_common_var
        use module_maths
        implicit none
        integer*8   , intent(in)      :: pos 
        integer*8   , intent(in)      :: Q0(4), Q00(4) ! == artsUpp, artsLow
        type (dta_SDF), intent(inout)    :: dta1
        character( 1)                    :: delta2branch
        integer*8              :: branch2delta
        integer*8              :: my_mol, delta
        !--------------------------------------------------
        my_mol = dta1%M
        !  NOTES: 
        !  (1) Prime and double primes refer, respectively, 
        !  to upper and lower states, respectively, i.e.
        !       * upper := 2; 0 ; or ' 
        !       * lower := 1;00 ; or " 
        !--------------------------------------------------
        IF ((my_mol .eq. 2).or.(my_mol .eq. 7)) then
        ! Group 2: Diatomic and Linear molecules.
        ! CO2( 2)
        ! ------------------------------
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
            dta1%J(pos,1)     = real(Q00(2),dp)
            dta1%N(pos,1)     = Q00(3)
            ! CAREFUL WITH THE SPIN!!! if it is integer that means something ASK RICHARD!
            dta1%espin(pos,1) = real(Q00(4),dp)
        !
        ! UPPER LEVEL:
            dta1%J(pos,2)     = real(Q0(2),dp)
            dta1%N(pos,2)     = Q0(3)
            dta1%espin(pos,2) = real(Q0(4),dp)
        ! BRANCH:
            delta = int(dta1%J(pos,2)-dta1%J(pos,1))
            dta1%br(pos) = delta2branch(delta,pos)
        !
        ELSE
            print*, " No vibrational band information or not speficied Format your selected molecule (HITRANid):", my_mol
        ENDIF

  END SUBROUTINE r_arts_LocalQ
!--------------------------------------------------------------------------------------------------------------------
  character function delta2branch(delta, pos)
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
        implicit none
        integer*8, intent(in)  :: delta, pos
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
                print*, "transition in position", pos,"do not follows the selection rules?"
                STOP
        endif
  end function delta2branch
!--------------------------------------------------------------------------------------------------------------------
  integer*8 function branch2delta(branch,pos)
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
!       (foe example: not every molecule has O- or S- branches)
!
! Author: Teresa Mendaza 02/03/2016
!--------------------------------------------------------------------------------------------------------------------
        implicit none
        character, intent(in)                 :: branch
        integer (kind=8), intent(in)          :: pos
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
            ! since HITRAN is an empirical DB it should not contain any line
            ! that does not follow the selection rules:
            ! \delta(J)=0, +-1, +-2
            ! \delta(l)= +-1
            print*, "transition in position", pos,"do not follows the selection rules?"
            STOP
        endif
  end function branch2delta
!--------------------------------------------------------------------------------------------------------------------
  subroutine systemQParam(sys,molP)
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
! Author: Teresa Mendaza 08/07/2016
!--------------------------------------------------------------------------------------------------------------------
            use module_common_var
            implicit none
            integer*8, intent(in)         :: sys(2)
            type (dta_MOL), intent(inout)      :: molP
            double precision                   :: T
            !-----------------------------------------
            T = molP % Temp
            !-----------------------------------------
            if ((sys(1) .eq. 2) .and. (sys(2) .eq. 22)) then
            ! CO2(2)-N2(4)
            ! Strow & Reuter 1970
            !molP%a1 = 0.02597D0,& !cm-1/atm
            !molP%a2 = 0.3307D0 ,& !cm-1
            !molP%a3 = 1.0350D0 ,& !cm-1
            !molP%dc = 2.2 ! Å (amstrong)
            ! -------------------------------
            ! Rodriguez et al. 1999; CO2 - N2
            !    molP%a1 = 0.0181D0 !cm-1/atm
            !    molP%a2 = 0.81D0   !cm-1
            !    molP%a3 = 0.008D0  !cm-1
            !    molP%dc = 2.2      !Å (amstrong)
            !    molP%ex1 = 0.85   
            !    molP%ex2 = 0.0152
            !
            !    if (T .ne. T0) then
            !        molP%a1 = molP%a1*(T0/T)**molP%ex1 
            !        molP%a2 = molP%a2*(T0/T)**molP%ex2
            !    endif
            ! -------------------------------
            ! Teresa Mendaza:
                molP%a1 = 3.678794E-01 !cm-1/atm
                molP%a2 =-2.629672E-12 !cm-1
                molP%a3 = 3.907117E-14 !cm-1
                molP%dc = 2.2      !Å (amstrong)
            
            else if ((sys(1) .eq. 2) .and. (sys(2) .eq. 7)) then
            ! CO2(2)-O2(7)
            ! -------------------------------
            ! Rodriguez et al. 1999; CO2 - O2
            !    molP%a1 = 0.0168D0 !cm-1/atm
            !    molP%a2 = 0.82D0   !cm-1
            !    molP%a3 = 0.007D0  !cm-1
            !    molP%dc = 2.4      !Å (amstrong)
            !    molP%ex1 = 0.50   
            !    molP%ex2 = -0.091
            !    if (T .ne. T0) then
            !        molP%a1 = molP%a1*(T0/T)**molP%ex1 
            !        molP%a2 = molP%a2*(T0/T)**molP%ex2
            !    endif
            ! -------------------------------
            ! Teresa Mendaza:
                molP%a1 = 3.678794E-01 !cm-1/atm
                molP%a2 =-2.629672E-12 !cm-1
                molP%a3 = 3.907117E-14 !cm-1
                molP%dc = 2.4      !Å (amstrong)


            else if ((sys(1) .eq. 7) .and. (sys(2) .eq. 7)) then
            ! Tran et al. 2006; O2 - O2
                molP%a1 = 0.0275D0 !cm-1/atm
                molP%a2 = 0.935D0  !cm-1
                molP%a3 = 1.01D0   !cm-1
                molP%dc = 1.01     !Å (amstrong)
                molP%ex1 = 1.0  
                molP%ex2 = 1.0
                if (T .ne. T0) then
                    molP%a1 = molP%a1*(T0/T)**molP%ex1 
                    molP%a2 = molP%a2*(T0/T)**molP%ex2
                endif
            else if ((sys(1) .eq. 7) .and. (sys(2) .eq. 22)) then
            ! Tran et al. 2006; O2 - N2
                molP%a1 = 0.0285D0 !cm-1/atm
                molP%a2 = 0.950D0  !cm-1
                molP%a3 = 1.03D0   !cm-1
                molP%dc = 1.0      !Å (amstrong)
                molP%ex1 = 1.0  
                molP%ex2 = 1.0
                if (T .ne. T0) then
                    molP%a1 = molP%a1*(T0/T)**molP%ex1 
                    molP%a2 = molP%a2*(T0/T)**molP%ex2
                endif
            else if ((sys(1) .eq. 4) .and. (sys(2) .eq. 7)) then
            ! Hartmann et al. 1999; N2O - O2
                molP%a1 = 0.0147D0 !cm-1/atm
                molP%a2 = 0.77D0   !cm-1
                molP%a3 = 0.025D0  !cm-1
                molP%dc = 2.9      !Å (amstrong)
                molP%ex1 = 0.85   
                molP%ex2 = 1.0
                if (T .ne. T0) then
                    molP%a1 = molP%a1*(T0/T)**molP%ex1 
                    molP%a2 = molP%a2*(T0/T)**molP%ex2
                endif
            else if ((sys(1) .eq. 4) .and. (sys(2) .eq. 22)) then
            ! Hartmann et al. 1999; N2O - N2
                molP%a1 = 0.0174D0 !cm-1/atm
                molP%a2 = 0.77D0   !cm-1
                molP%a3 = 0.025D0  !cm-1
                molP%dc = 2.9      !Å (amstrong)
                molP%ex1 = 0.85   
                molP%ex2 = 1.0
                if (T .ne. T0) then
                    molP%a1 = molP%a1*(T0/T)**molP%ex1 
                    molP%a2 = molP%a2*(T0/T)**molP%ex2
                endif
            else
                print*, "systemQParam: No adjusted parameter for the Q basis rate avaible."
                print*, "Please, update the data base."
            endif
    end subroutine systemQParam
!--------------------------------------------------------------------------------------------------------------------
