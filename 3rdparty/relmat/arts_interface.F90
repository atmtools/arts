!--------------------------------------------------------------------------------------------------------------------
module arts_interface
!--------------------------------------------------------------------------------------------------------------------
! This module contains subroutines to interact with ARTS 
!
    interface

        subroutine RM_LM_tmc_arts(nLines, sgmin, sgmax, &
                          artsM, artsI, artsWNO, &
                          artsS, artsGA, artsE00, &
                          artsNA, artsUpp, artsLow, &
                          artsg0 , artsg00, &
                          T, Ptot, QT, QT0, mass, &
                          npert, pert, i_pert, p_mass, p_vmr,&
                          W_rn, dipo, rho) 
        !
        !   INPUT VARIABLES
        !   ----------------------
        !   nLines   : These are the number of lines identified as belonging to the band.
        !   sgmin    : The minimum frequency in cm-1 [defaults to the first value of artsWNO]
        !   sgmax    : The maximum frequency in cm-1 [defaults to the last value of artsWNO]
        !   artsM    : The HITRAN molecule number (e.g., 7 is O2)
        !   artsI    : The HITRAN isotope number
        !   artsWNO  : frequency in cm-1 in vacuum of the lines.  Sorted by lowest first
        !   artsS    : Intensity of the line at the same temperature as QT0 but abundance has already been considered
        !   artsGA   : Air pressure broadening in cm-1
        !   arts E00 : Lower state energy in cm-1
        !   artsNA   : Pressure broadening constant in cm-1
        !   artsUpp  : Upper state quantum numbers. First is L2, then is J, then is N, then is S.  
        !              If the quantum number is not applicable, the position contains the number -1.
        !   artsLow  : Same as artsUpp but for lower state numbers.
        !   artsg0   : Upper state g-constant
        !   artsg00  : Lower state g-constant
        !   T        : Temperature of the system in Kelvin.
        !   Ptot     : Total pressure in Pascal.
        !   QT       : Partition function at Temperature T.
        !   QT0      : Partition function at same temperature as line intensity.
        !   mass     : Mass of the molecule.
        !   npert    :  Number of perturbers/colliders.  Initially assume this is 2 for N2 and O2.
        !   pert     : Array containing HitranID of the system perturber (colliders). Index of perturbing molecule in HITRAN. Same as artsM. 
        !              Assumed position 0 is O2-66 and position 1 is N2-44.
        !   i_pert   : Array of the colliders's isotopes. Index of isotope in HITRAN.  Same as artsI.  
        !              Assumed position 0 is O2-66 and position 1 is N2-44.
        !   p_mass   : mass array of length npert of perturbing gases.  
        !              Assumed position 0 is for O2-66 and position 1 is for N2-44.
        !   p_vmr    : volume mixing ratio (VMR) of perturbing gases.  
        !              Assumed position 0 is for O2-66 and position 1 is for N2-44.
        !
        !   OUTPUT VARIABLES
        !   -----------------
        !   W_rn     : The relaxation matrix.  For output. 
        !   dipo     : The dipole moments.  For output.
        !   rho      : Populations.  For output.
        ! 
        !!!!! MODULES IN USE:
            use module_common_var
            use module_error
            !use module_OC_files
                use module_maths
                use module_molecSp
                use module_read
                    use module_phsub
                    use module_linemixing

            Implicit none
            !   INPUT variables
            integer*8         :: nLines, npert
            integer*8         :: artsM, artsI
            integer*8         :: artsg0(nLines), artsg00(nLines)
            integer*8         :: artsLow(4,nLines), artsUpp(4,nLines)
            integer*8         :: pert(npert), i_pert(npert)
            Double Precision  :: sgmin, sgmax, T, Ptot
            Double Precision  :: QT, QT0, mass
            Double Precision  :: p_vmr(npert), p_mass(npert)
            Double Precision  :: artsWNO(nLines),artsS(nLines), &
                                 artsGA(nLines), artsE00(nLines), &
                                 artsNA(nLines)
            !   OUTPUT variables
            Double Precision  :: rho(nLines), dipo(nLines) 
            Double Precision  :: W_rn(nLines,nLines)
        end subroutine RM_LM_tmc_arts

        subroutine InitW(n,W)
            implicit none
            integer*8 , intent(in ) :: n !== Molecules isotope(AFGL code)
            Double Precision, intent(out) :: W(n,n)
        end subroutine InitW

        subroutine add2Wfinal(n,Wfinal,Wadd,xMol)
            implicit none
            integer*8, intent (in)   :: n 
            Double Precision, intent (in)   :: xMOL
            Double Precision, intent (in)   :: Wadd(n,n)
            Double Precision, intent(inout) :: Wfinal(n,n)
        end subroutine add2Wfinal
        
        subroutine show_W(n,Wfinal)
            use module_common_var
            implicit none
            integer*8, intent (in)   :: n 
            Double Precision, intent (in)   :: Wfinal(n,n)
        end subroutine show_W


    end interface

END module arts_interface
!--------------------------------------------------------------------------------------------------------------------
! ARTS driver SUBROUTINES -------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------------
SUBROUTINE RM_LM_tmc_arts(nLines, sgmin, sgmax, &
                          artsM, artsI, artsWNO, &
                          artsS, artsGA, artsE00, &
                          artsNA, artsUpp, artsLow, &
                          artsg0 , artsg00, &
                          T, Ptot, QT, QT0, mass, &
                          npert, pert, i_pert, p_mass, p_vmr,&
                          W_rn, dipo, rho) bind(C, name='arts_relmat_interface')
!--------------------------------------------------------------------------------------------------------------------
!
! This SUBROUTINE is used for computing the following variables:
!	Dipole and e-level population	(Spectroscopy Data)
!   Relaxation Matrix for a given T	(Relaxation Matrix)
! LINEAR MOLECULES.
!
!
!   INPUT VARIABLES
!   ----------------------
!   nLines   : These are the number of lines identified as belonging to the band.
!   sgmin    : The minimum frequency in cm-1 [defaults to the first value of artsWNO]
!   sgmax    : The maximum frequency in cm-1 [defaults to the last value of artsWNO]
!   artsM    : The HITRAN molecule number (e.g., 7 is O2)
!   artsI    : The HITRAN isotope number
!   artsWNO  : frequency in cm-1 in vacuum of the lines.  Sorted by lowest first
!   artsS    : Intensity of the line at the same temperature as QT0 but abundance has already been considered
!   artsGA   : Air pressure broadening in cm-1
!   arts E00 : Lower state energy in cm-1
!   artsNA   : Pressure broadening constant in cm-1
!   artsUpp  : Upper state quantum numbers. First is L2, then is J, then is N, then is S.  
!              If the quantum number is not applicable, the position contains the number -1.
!   artsLow  : Same as artsUpp but for lower state numbers.
!   artsg0   : Upper state g-constant
!   artsg00  : Lower state g-constant
!   T        : Temperature of the system in Kelvin.
!   Ptot     : Total pressure in Pascal.
!   QT       : Partition function at Temperature T.
!   QT0      : Partition function at same temperature as line intensity.
!   mass     : Mass of the molecule.
!   npert    : Number of perturbers/colliders.  Initially assume this is 2 for N2 and O2.
!   pert     : Array containing HitranID of the system perturber (colliders). Index of perturbing molecule in HITRAN. Same as artsM. 
!              Assumed position 0 is O2-66 and position 1 is N2-44.
!   i_pert   : Array of the colliders's isotopes. Index of isotope in HITRAN.  Same as artsI.  
!              Assumed position 0 is O2-66 and position 1 is N2-44.
!   p_vmr    : volume mixing ratio (VMR) of perturbing gases.  
!              Assumed position 0 is for O2-66 and position 1 is for N2-44.
!   p_mass   : mass array of length npert of perturbing gases.  
!              Assumed position 0 is for O2-66 and position 1 is for N2-44.
!
!   OUTPUT VARIABLES
!   -----------------
!   W_rn     : Renormalized Relaxation Matrix acording to the renormalization procedure described in
!              Niro et al. (2004). 
!   dipo     : Dipole transition Moments of the Lines.
!   rho      : Populations of the Lower Levels of the Lines at 296 K.  
! 
!
!	Accessed Files:	 'none'
!	---------------
!
!	Called Routines: 'VarInit'  (VARiable INITialization)
!	---------------  'MoleculeID' (Molecule IDentification)
!                    'Readline' (Read HITRAN12 CH4 file)
!                    'DipCAL'   (DIPole elements CALculation)
!                    'systemQParam' (External Parameters of the system)
!                    'WelCAL'   (W elements CALculation)
!                    'RN_Wmat' (ReNormalization of W matrix)
!
!
!	T. Mendaza last change 22 Jun 2016
!-------------------------------------------------------------------
!
! MODULES IN USE:
    use module_common_var
    use module_error
    !use module_OC_files
        use module_maths
        use module_molecSp
        use module_read
           use module_phsub
            use module_linemixing

    Implicit none
!   INPUT variables
    integer*8 :: nLines, npert
    integer*8 :: artsM, artsI
    integer*8 :: artsg0(nLines), artsg00(nLines)
    integer*8 :: artsLow(4,nLines), artsUpp(4,nLines)
    integer*8 :: pert(npert), i_pert(npert)
    Double Precision  :: sgmin, sgmax, T, Ptot
    Double Precision  :: QT, QT0, mass
    Double Precision  :: p_vmr(npert), p_mass(npert)
    Double Precision  :: artsWNO(nLines),artsS(nLines), &
                         artsGA(nLines), artsE00(nLines), &
                         artsNA(nLines)
    integer*8              :: dta_size1, dta_size2
    integer*8              :: iLine, i, j, k
!   OUTPUT variables
    Double Precision       :: rho(nLines), dipo(nLines) 
    Double Precision       :: W_rn(nLines,nLines)
    Double Precision       :: Wmat(nWmax,nWmax),&
                              Wper(nWmax,nWmax)
!   Double Precision       :: Y_RosT(nLmx)
    Double Precision       :: xMOLp(npert)
    type (dta_SDF)         :: dta1
    type (dta_RMF)         :: dta2
    type (dta_MOL)         :: molP
    type (dta_MOL)         :: PerM
    integer*8              :: sys(2)  
    character*6            :: perturber(2)
    logical                :: enough_Lines
!
    write(*,2016), T
2016 Format("Starting Linemixing Relaxation Matrix software. T=",f5.0,"K")
!----------
! Band quantities specification
    print*, 'Init. Variables...'
    CALL VarInit(dta1,dta2,molP)
    molP % Temp = T
    molP % Ptot = Ptot
!
!----------
! Obtainig the molecule ID from the Formula specified in "module_common_var"
    print*, 'Identifying molecule and loading its parameters...'
    !CALL moleculeID(artsM(1), artsI(1),molP) 
    CALL moleculeID(artsM, artsI, mass, QT, QT0, .true., molP)
    dta1%M = molP%M
!---------
! Call for reading HITRAN spectroscopy data.
!
    print*, 'Reading HITRAN12 File...'
    !CALL readline(dta1, dta_size1)
    CALL Hit2DTA(dta1, dta_size1, nLines, enough_Lines, &
                                            artsWNO, &
                                            artsS, &
                                            artsGA, &
                                            artsE00, &
                                            artsNA , &
                                            artsUpp, artsLow, &
                                            artsg0 , artsg00, &
                                            sgmin  , sgmax)
    nLines = dta_size1 !nLines = 856 #raw in CH4_2nu3_HIT12.dat
    dta_size2=nLines**2 !== dta_size1*dta_size1
    !print*, "nLines:", nLines, dta_size1
    !stop 
!---------
! Compute the relative population of the lower state
! at Temperature T0
    call PopuCAL(dta1,nLines, molP)
    do j = 1, nLines
        if (T .eq. T0) then
            rho(j) = dta1 % PopuT0(j)
        else
            rho(j) = dta1 % PopuT(j)
        endif
    enddo
! NOTE: we use 'tra' mode (see 'PopuT0') because we are producing 
! Input files to Ha Tran Line-mixing code.
!---------
! Calculate Dipole element for each line.
!
    print*, 'Calculating Dipole moment'
    CALL DipCAL(dta1,nLines,molP)
    do j = 1, nLines
        dipo(j) = dta1 % D0(j)
    enddo
!---------
! Write SDF file
!
    print*,"Looping over system of perturbers..."
    print*,"--------------------------------->"

!
!   Init W rel-mat.
    CALL InitW(nWmax,Wmat)
!   Looping:
    DO i = 1,npert
        !
        ! Identifying Perturbers Molecule and
        ! its ATM concentration (renorm to 100%):
        xMOLp(i) = p_vmr(i)*100/sum(p_vmr)
        !
        CALL InitW(nWmax,Wper)
        !---------
        ! Identifying perturber molecule.
        !
        ! Perturber Molecule: This molecule has to be 
        ! significatively faster than the molecule in study.
        print*, '>Identifying perturber molecule...'
        !CALL moleculeID(pert(i),i_pert(i), PerM)
        CALL moleculeID(pert(i), i_pert(i), p_mass(i), 0.0_dp, 0.0_dp, .false. ,PerM)
        !----------
        ! let's take th proper a1, a2, a3, dc adjust parameters for the
        ! system:
        ! sys = "CO2-N2"
        ! or
        ! sys = "CO2-O2"
        !
        sys(1) = molP%M
        sys(2) = PerM%M
        CALL systemQParam(sys,molP)
        perturber(i)=PerM%chmol
        !---------
        ! Obtain Relaxation matrix elements for each line.
        !
        print*, '>Building Relaxation Matrix...'
        CALL WelCAL(dta1, nLines, molP, PerM, Wper)
        !---------
        ! Adding the corresponding perturber-molecule 
        ! contribution to the relaxation matrix.
        !
        call add2Wfinal(nLines,Wmat,Wper,xMolp(i))
    ENDDO
!---------
    print*,"<--------------------------Finished loop"
    !CALL show_W(nLines,Wmat)
    !stop
!---------
! Renormalization of the Relaxation matrix.
!
    CALL InitW(nLines,W_rn)
    print*, 'Renormalization procedure of the RM...'
    CALL RN_Wmat(nLines, dta1, Wmat, W_rn) 
!
!
    PRINT *, "Successful run!"
!   
    !STOP
!
  END SUBROUTINE RM_LM_tmc_arts
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE InitW(n,W)
!--------------------------------------------------------------------------------------------------------------------
! 
!
    implicit none
    integer*8, intent(in ) :: n !== Molecules isotope(AFGL code)
    Double Precision, intent(out)       :: W(n,n)
    integer*8              :: i,j,k

!----------
! 
    DO i = 1, n
        DO j = 1,n 
            W(i,j) = 0.0
        ENDDO
    ENDDO

  END SUBROUTINE InitW
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE add2Wfinal(n,Wfinal,Wadd,xMol)
!--------------------------------------------------------------------------------------------------------------------
! 
!
    implicit none
    integer*8, intent (in)   :: n 
    Double Precision, intent (in)         :: xMOL
    Double Precision, intent (in)         :: Wadd(n,n)
    Double Precision, intent(inout)       :: Wfinal(n,n)
    integer*8                :: i,j,k

!----------
! 
    DO i = 1, n
        DO j = 1,n 
            Wfinal(i,j) = Wfinal(i,j) + xMOL*Wadd(i,j)
        ENDDO
    ENDDO

  END SUBROUTINE add2Wfinal
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE show_W(n,Wfinal)
!--------------------------------------------------------------------------------------------------------------------
! 
!
    use module_common_var
    implicit none
    integer*8, intent (in)   :: n 
    Double Precision, intent (in)         :: Wfinal(n,n)
    integer*8                :: i,j,k

!----------
! 
    DO i = 1, n
        DO j = 1,n 
            write(*,1000), "W(",i,",",j,")=",Wfinal(i,j) 
1000 Format(a2,i3,a1,i3,a3,E12.3)
        ENDDO
    ENDDO

  END SUBROUTINE show_W
!--------------------------------------------------------------------------------------------------------------------
