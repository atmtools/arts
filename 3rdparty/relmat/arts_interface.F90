!--------------------------------------------------------------------------------------------------------------------
module arts_interface
!--------------------------------------------------------------------------------------------------------------------
! This module contains the "RM_LM_tmc_arts" and "RM_LM_LLS_tmc_arts" subroutines which interface with ARTS. 
! It also contains subroutines to do internal checks on models' outputs. 
!
    interface

        subroutine RM_LM_tmc_arts(nLines, sgmin, sgmax, &
                          artsM, artsI, artsWNO, &
                          artsS, artsGA, artsE00, &
                          artsNA, artsUpp, artsLow, &
                          artsg0 , artsg00, &
                          T, Ptot, QT, QT0, mass, &
                          npert, pert, i_pert, p_mass, p_vmr,&
                          runE_deb,ordered,&
                          W_rn, dipo, rho, &
                          Y1,Y2,Y3) 
        !   MODEL:
        !   -----
        !   This subroutine contains a version of the HARTMANN-NIRO Linemixing model.
        !
        !   INPUT VARIABLES
        !   ----------------------
        !   nLines   : This is the number of lines identified as belonging to the band (INTEGER).
        !   sgmin    : The minimum frequency in cm-1 [defaults to the first value of artsWNO; REAL*8]
        !   sgmax    : The maximum frequency in cm-1 [defaults to the last value of artsWNO; REAL*8]
        !   artsM    : The HITRAN molecule number (e.g., 7 is O2). (INTEGER).
        !   artsI    : The HITRAN isotope number (INTEGER).
        !   artsWNO  : frequency in cm-1 in vacuum of the lines.  Sorted by lowest first
        !   artsS    : Intensity of the line at the same temperature as QT0 but abundance has already been considered.
        !              (REAL*8 Array of length NLINES).
        !   artsGA   : Air pressure broadening in cm-1/atm halfwidth (REAL*8 Array of length NLINES).
        !   artsE00  : Lower state energy in cm-1 (REAL*8 Array of length NLINES).
        !   artsNA   : Pressure broadening constant in cm-1 (REAL*8 Array of length NLINES).
        !   artsUpp  : REAL*8 Array containing the UPPER STATE quantum numbers (length=NLINES). 
        !              First is L2, then is J, then is N, then is S.  
        !              If a quantum number is not applicable, the position contains the number -1.
        !   artsLow  : Same as artsUpp but for LOWER STATE numbers (REAL*8 Array of length NLINES).
        !   artsg0   : REAL*8 Array of the Upper state g-constant per band-line (length NLINES).
        !   artsg00  : REAL*8 Array of the Lower state g-constant per band-line (length NLINES).
        !   T        : REAL*8 Temperature of the system in Kelvin.
        !   Ptot     : REAL*8 Total pressure in Pascal.
        !   QT       : REAL*8 Partition function at Temperature T.
        !   QT0      : REAL*8 Partition function at same temperature as line intensity.
        !   mass     : REAL*8 Mass of the molecule.
        !   npert    : INTEGER Number of perturbers/colliders.  Initially assume this is 2 for N2 and O2.
        !   pert     : INTEGER Array (length= NPERT) containing HitranID of the system perturber (colliders). 
        !              Index of perturbing molecule in HITRAN. Same as artsM. 
        !              *Default case*: Assumed position 0 is O2-66 and position 1 is N2-44.
        !   i_pert   : INTEGER Array of the colliders's isotopes (length= NPERT). Index of isotope in HITRAN.  
        !              Same as artsI.  
        !              *Default case*: Assumed position 0 is O2-66 and position 1 is N2-44.
        !   p_mass   : REAL*8 mass array of perturbing gases (length=NPERT).  
        !              *Default case*: Assumed position 0 is for O2-66 and position 1 is for N2-44.
        !   p_vmr    : REAL*8 volume mixing ratio (VMR) of perturbing gases (length=NPERT).  
        !              *Default case*: Assumed position 0 is for O2-66 and position 1 is for N2-44.
        !   runE_deb : INTEGER that allows to print on screen (1) or not (0). It is used on debugging.
        !   ordered  : INTEGER which includes the selection in output:
        !              -1 -> Returns W, but no Y_Ro is returned.
        !               0 -> W Diagonal, and no Y_Ro is returned.
        !               1 -> Returns W, and just the first order coefficient Y1 is returned = 1st order linemixing.
        !               2 -> Returns W, and all Y1,Y2,Y3 arrays = 2nd order linemixing.
        !                    1st column = Rosenkranz parameter/line
        !                    2nd column = g (2nd ord-parameter)/line
        !                    3rd column = dv(2nd ord-shift)/line
        !
        !
        !   OUTPUT VARIABLES
        !   -----------------
        !   W_rn     : The relaxation matrix.  
        !   Y1       : 1st order linemixing parameter. 
        !   Y2       : 2nd order linemixing parameter. 
        !   Y3       : 3rd order linemixing parameter. 
        !   dipo     : The dipole moments.  
        !   rho      : Populations.  
        ! 
        !---MODULES IN USE:
            use module_common_var
            use module_error
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
            integer*8         :: ordered
            integer*8         :: runE_deb
            Double Precision  :: sgmin, sgmax, T, Ptot
            Double Precision  :: QT, QT0, mass
            Double Precision  :: p_vmr(npert), p_mass(npert)
            Double Precision  :: artsWNO(nLines),artsS(nLines), &
                                 artsGA(nLines), artsE00(nLines), &
                                 artsNA(nLines)
            !   OUTPUT variables
            Double Precision  :: rho(nLines), dipo(nLines)
            Double Precision  :: Y1(nLines),Y2(nLines),Y3(nLines) 
            Double Precision  :: W_rn(nLines,nLines)
        end subroutine RM_LM_tmc_arts

        subroutine RM_LM_LLS_tmc_arts(nLines, sgmin, sgmax, &
                          artsM, artsI, artsWNO, &
                          artsS, artsGA, artsE00, &
                          artsNA, artsUpp, artsLow, &
                          artsg0 , artsg00, &
                          T, Ptot, QT, QT0, mass, &
                          npert, pert, i_pert, p_mass, p_vmr,&
                          runE_deb,ordered,&
                          W_rn, dipo, rho, &
                          Y1,Y2,Y3) 
        !   MODEL:
        !   -----
        !   This subroutine contains the MENDAZA's Linemixing model which computes 
        !   the Relaxation Matrix, Rosenkranz parameter and 2nd order linemixing 
        !   parameters within the ECS approach for a given linear molecule 
        !   in a system of colliders.
        !
        !   INPUT VARIABLES
        !   ----------------------
        !   Same as in "RM_LM_tmc_arts"
        !
        !   OUTPUT VARIABLES
        !   -----------------
        !   Same as in "RM_LM_tmc_arts"
        ! 
        !---MODULES IN USE:
            use module_common_var
            use module_error
                use module_maths
                use module_molecSp
                use module_read
                    use module_phsub
                    use module_LLS
                    use module_linemixing

            Implicit none
            !   INPUT variables
            integer*8, intent (in) :: nLines, npert
            integer*8, intent (in) :: artsM, artsI
            integer*8, intent (in) :: artsg0(nLines), artsg00(nLines)
            integer*8, intent (in) :: artsLow(4,nLines), artsUpp(4,nLines)
            integer*8, intent (in) :: pert(npert), i_pert(npert)
            integer*8, intent (in) :: ordered
            integer*8, intent (inout) :: runE_deb
            Double Precision, intent (in) :: sgmin, sgmax, T, Ptot
            Double Precision, intent (in) :: QT, QT0, mass
            Double Precision, intent (in) :: p_vmr(npert), p_mass(npert)
            Double Precision, intent (in) :: artsWNO(nLines),artsS(nLines), &
                                             artsGA(nLines), artsE00(nLines), &
                                             artsNA(nLines)
            !   OUTPUT variables
            Double Precision, intent (out):: rho(nLines), dipo(nLines)
            Double Precision, intent (out):: Y1(nLines),Y2(nLines),Y3(nLines)
            Double Precision, intent (out):: W_rn(nLines,nLines)
        end subroutine RM_LM_LLS_tmc_arts

        subroutine alloSDF(n,dta1)
            use module_common_var
            implicit none
            integer*8, intent(in )        :: n 
            type (dta_SDF),intent(out)    :: dta1
        end subroutine alloSDF

        subroutine alloRMF(n,dta2)
            use module_common_var
            implicit none
            integer*8     , intent(in )   :: n 
            type (dta_RMF), intent(out)   :: dta2
        end subroutine alloRMF

        subroutine VarInit(molP,econ,runE)
            use module_common_var
            implicit none
            integer*8     , intent(in)    :: runE
            type (dta_MOL), intent(inout) :: molP
            type (dta_ERR), intent(inout) :: econ
        end subroutine VarInit

        subroutine mol_Init(molP)
            use module_common_var
            implicit none
            type (dta_MOL), intent(inout) :: molP
        end subroutine mol_Init

        subroutine InitM(n,m,W)
            implicit none
            integer*8 , intent(in ) :: n,m !== Molecules isotope(AFGL code)
            Double Precision, intent(out) :: W(n,m)
        end subroutine InitM

        subroutine includeW(n,indx,W_rn, NA, GA, rTT0, P, n0, Wrno) 
            implicit none
            integer*8, intent(in)         :: n, n0 
            integer*8, intent(in)         :: indx(n) 
            Double Precision, intent(in ) :: NA(n), GA(n), rTT0, P
            Double Precision, intent(in ) :: Wrno(n0,n0)
            Double Precision, intent(out) :: W_rn(n,n)
        end subroutine includeW

        subroutine includeY(n,indx,Yf, n0, Yc)
            implicit none
            integer*8, intent(in )        :: n !number of lines input from ARTS 
            integer*8, intent(in )        :: n0!number of lines with proper quantum numbers
            integer*8, intent(in )        :: indx(n) 
            Double Precision, intent(in ) :: Yc(n0)
            Double Precision, intent(out) :: Yf(n)
        end subroutine includeY

        subroutine just_fill_DiagWRn(n, NA, GA, rTT0, P, Wrn) 
            implicit none
            integer*8, intent(in )        :: n 
            Double Precision, intent(in)  :: NA(n), GA(n), rTT0, P
            Double Precision, intent(out) :: Wrn(n,n)
        end subroutine just_fill_DiagWRn

        subroutine add2Wfinal(n,Wfinal,Wadd,xMol)
            implicit none
            integer*8, intent (in)        :: n 
            Double Precision,intent (in)  :: xMOL
            Double Precision,intent (in)  :: Wadd(n,n)
            Double Precision,intent(inout):: Wfinal(n,n)
        end subroutine add2Wfinal
        
        subroutine show_W(n,Wfinal)
            use module_common_var
            implicit none
            integer*8, intent (in)        :: n 
            Double Precision, intent (in) :: Wfinal(n,n)
        end subroutine show_W

        subroutine show_PD(n,wno,p,d)
            use module_common_var
            implicit none
            integer*8, intent (in)        :: n 
            Double Precision, intent (in) :: wno(n), p(n), d(n)
        end subroutine show_PD

        subroutine save_W2plot(n, dta1, dta2 ,molP, npert, pert, econ, model)
            use module_common_var
            use module_molecSP
            use module_error
            implicit none
            integer*8, intent (in)        :: n, npert !n=dta_size1
            type (dta_MOL), intent (in)   :: molP
            type (dta_SDF), intent (in)   :: dta1
            type (dta_RMF), intent (in)   :: dta2
            type (dta_ERR), intent (in)   :: econ
            character (6) , intent (in)   :: pert(npert)
            character (3) , intent (in)   :: model
        end subroutine save_W2plot

        subroutine save_Yrp(dta1, n, molP, Y_RosenP, model)
            use module_common_var

            implicit none
            integer*8  , intent (in)      :: n
            double precision, intent (in) :: Y_RosenP(n)
            type (dta_SDF)  , intent (in) :: dta1
            type (dta_MOL)  , intent (in) :: molP
            character(3)    , intent (in) :: model
        end subroutine save_Yrp


    end interface

END module arts_interface
!--------------------------------------------------------------------------------------------------------------------
! ARTS driver SUBROUTINES -------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------------
SUBROUTINE RM_LM_tmc_arts(nLines, sgmin, sgmax, &
                          artsM, artsI, artsWNO, &
                          artsS, artsGA, artsE00, &
                          artsNA, artsUpp, artsLow, &
                          artsg0 , artsg00, &
                          T, Ptot, QT, QT0, mass, &
                          npert, pert, i_pert, p_mass, p_vmr,&
                          runE_deb,ordered,&
                          W_rn, dipo, rho, &
                          Y1,Y2,Y3) bind(C, name='arts_relmat_interface__hartmann_and_niro_type')
!--------------------------------------------------------------------------------------------------------------------
!
! This SUBROUTINE is used to compute the following variables:
!	*Dipole and e-level population	    (Spectroscopy Data)
!   *Relaxation Matrix for a given T	(Relaxation Matrix)
!   *Rosenkranz parameter for a 
!              given pair (P,T)         (Rosenkranz's Para)
!   *Second order linemixing P. 
!              on given (P,T)           (2nd order Linemix)
!   on LINEAR MOLECULES.
!
!   NOTE: Check variables up in the interface.
!   ----
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
!	T. Mendaza last change 20 Feb 2017
!-------------------------------------------------------------------
!
! MODULES IN USE:
    use module_common_var
    use module_error
        use module_maths
        use module_molecSp
        use module_read
           use module_phsub
            use module_linemixing

    Implicit none
!   INPUT variables
    integer*8, intent(in) :: nLines, npert
    integer*8, intent(in) :: artsM, artsI
    integer*8, intent(in) :: artsg0(nLines), artsg00(nLines)
    integer*8, intent(in) :: artsLow(4,nLines), artsUpp(4,nLines)
    integer*8, intent(in) :: pert(npert), i_pert(npert)
    integer*8, intent(in) :: ordered
    integer*8, intent(inout) :: runE_deb
    Double Precision, intent(in)  :: sgmin, sgmax, T, Ptot
    Double Precision, intent(in)  :: QT, QT0, mass
    Double Precision, intent(in)  :: p_vmr(npert), p_mass(npert)
    Double Precision, intent(in)  :: artsWNO(nLines),artsS(nLines), &
                                     artsGA(nLines), artsE00(nLines), &
                                     artsNA(nLines)
!
!   OUTPUT variables
    Double Precision, intent(out) :: rho(nLines), dipo(nLines)
    Double Precision, intent(out) :: Y1(nLines),Y2(nLines),Y3(nLines)
    Double Precision, intent(out) :: W_rn(nLines,nLines)
!OTHER VARIABLES
    integer*8              :: dta_size1, IERR1, IERR2, IERR3, IERR4
    integer*8              :: iLine, i, j, k
    Double Precision, ALLOCATABLE :: Wmat(:,:),&
                                     Wper(:,:),&
                                     Wrno(:,:)
    Double Precision, ALLOCATABLE :: Y_RosT(:),Y_G(:),Y_DV(:)
    Double Precision       :: xMOLp(npert)
    integer*8              :: vLines_Indx(nLines)
    Double Precision       :: faH, rT
    type (dta_SDF), target :: dta1
    type (dta_SDF), pointer:: pd1
    type (dta_RMF), target :: dta2
    type (dta_RMF), pointer:: pd2
    type (dta_MOL)         :: molP
    type (dta_MOL)         :: PerM
    type (dta_ERR)         :: econtrol
    integer*8              :: sys(2), lM, lP,aupp,alow
    character*6            :: perturber(2)
    character*2            :: fmt1, fmt2
    logical                :: enough_Lines, my_print
!   
! ROUTINE:
!
!----------
! Inital common temperature-dependent constants
    rT = T/T0
!
    if (runE_deb .ge. 1) then 
        write(*,*)  "RELMAT RUN-TYPE = Verbose."
        write(*,*)  "Type: Hartmann and Niro"
        write(*,2016) T, Ptot
    endif
2016 Format("Starting Linemixing Relaxation Matrix software. T=",f5.0,"K; P=",f7.2," atm")
!
!----------
! Variable Allocation
    if (runE_deb .ge. 1) write(*,'(a33,i4,a1)') 'Allocate SDF and RMF variables to arrays (', nLines,')'
    CALL alloSDF(nLines, dta1)
    pd1 => dta1
!
!----------
! Band quantities specification
    if (runE_deb .ge. 1) write(*,*) 'Init. Variables...'
    CALL VarInit(molP,econtrol,runE_deb)
    molP % Temp = T !Kelvin
    molP % Ptot = Ptot!atm
    molP % QTy  = "REG"
    my_print = .false.
!
!----------
! Obtainig the molecule ID from the Formula specified in "module_common_var"
    if (econtrol % e(1) .ge. 1) then
        write(*,*) 'Identifying molecule and loading its parameters...'
    endif
    !
    CALL moleculeID(artsM, artsI, mass, QT, QT0, .true., molP, econtrol)
    dta1%M = molP%M
!
!---------
! Call for reading HITRAN spectroscopy data.
    if (econtrol % e(1) .ge. 1) write(*,*) 'Locating Band information...'
    CALL Hit2DTA(pd1, dta_size1, nLines, vLines_Indx, &
                                            artsWNO, &
                                            artsS, &
                                            artsGA, &
                                            artsE00, &
                                            artsNA , &
                                            artsUpp, artsLow, &
                                            artsg0 , artsg00, &
                                            sgmin  , sgmax, &
                                            econtrol)
!
!---------
! Compute the relative population of the lower state
! at Temperature T0
    if (econtrol % e(1) .ge. 1) write(*,*) 'Counting band-lines population'
    call PopuCAL(pd1,dta_size1, molP, econtrol)
    !
    do j = 1, nLines
        if (vLines_Indx(j) .eq. 0) then
            if (T .eq. T0) then
                rho(j) = artsg00(j)*dexp(-c2*artsE00(j)/T0)/QT0
            else
                rho(j) = artsg00(j)*dexp(-c2*artsE00(j)/T0)* &
                         dexp(-c2*artsE00(k)*(1.d0/T-1.d0/T0))/QT
            endif
        else
            if (T .eq. T0) then
                rho(j) = dta1 % PopuT0(vLines_Indx(j))
            else
                rho(j) = dta1 % PopuT(vLines_Indx(j))
            endif
        endif 
    enddo
! NOTE: the code uses 'tra' mode as a default option (see 'PopuT0'). 
!       One can change this mode in "module_common_var" (check options there)
!---------
! Calculate Dipole element for each line.
    if (econtrol % e(1) .ge. 1) write(*,*) 'Calculating Dipole moment'
    CALL DipCAL(pd1,dta_size1,molP,econtrol)
!
! To save an ASCII file and shop the program Uncomment the following line: 
    !call show_PD(nLines, dta1%Sig(1:dta_size1), dta1 % PopuT0(1:dta_size1), pd1 % D0(1:dta_size1))
!
    do j = 1, nLines
        if (vLines_Indx(j) .eq. 0) then
            dipo(j) = dsqrt(artsS(k)/(artsWNO(k)* &
                        rho(k)*(1.D0 - &
                        dexp(-c2*artsWNO(k)/T0))))
        else 
            dipo(j) = dta1%DipoT(vLines_Indx(j)) 
        endif 
    enddo    
!
! Uncomment the next two lines to print POPULATION & DIPOLE to the screen:
    !call show_PD(nLines, artsWNO, rho, dipo)
    !call show_PD(dta_size1, dta1%sig, dta1%PopuT, dta1%D0)
!---------
! RELAXATION MATRIX CALCULATION STARS!
! 1) if the no. of lines of the Band meet "rule1" 
! 2) and the user has not selected returning a diagonal W (ordered==0) 
! then the process starts
    if (econtrol % e(1) .ge. 1) write(*,*) 'Relaxation Matrix calculation...'
! 
    if (rule1(dta_size1) .and. .not.(ordered .eq. 0)) then
        if (econtrol % e(1) .ge. 1) then
            write(*,*)"Looping over system of perturbers..."
            write(*,*)"----------------------------------->"
        endif
    !
    ! Allocate variables according to the valid number of lines:
    allocate(Wmat(dta_size1,dta_size1),STAT = IERR1)
    IF (IERR1 .ne. 0) call memoError(" Wmat ",econtrol)
    allocate(Wper(dta_size1,dta_size1),STAT = IERR2)
    IF (IERR2 .ne. 0) call memoError(" Wper ",econtrol)
    !
    !   Init W rel-mat.
        CALL InitM(dta_size1,dta_size1,Wmat)
    !   Looping:
        DO i = 1,npert
        !
        ! Identifying Perturbers Molecule and
        ! its ATM concentration (value must be from 0-1):
            xMOLp(i) = (p_vmr(i)/sum(p_vmr))
        !
            CALL mol_Init(PerM)
        !---------
        ! Identifying perturber molecule.
        !
        ! Perturber Molecule: This molecule has to be 
        ! significatively faster than the molecule in study.
            if (econtrol % e(1) .ge. 1) write(*,*) '>Identifying perturber molecule...'
            CALL moleculeID(pert(i), i_pert(i), p_mass(i), 0.0_dp, 0.0_dp, .false., &
                            PerM, econtrol)
        !----------
        ! let's take th proper a1, a2, a3, dc adjust parameters for the
        ! system:
            if (econtrol % e(1) .ge. 1) then
                write(*,*)  ">>System: {",trim(molP%chmol)," - ",trim(PerM%chmol),"}"
            endif
        ! sys = "CO2-N2"
        ! or
        ! sys = "CO2-O2"
        !
            sys(1) = molP%M
            sys(2) = PerM%M
            CALL systemQParam(sys,molP,econtrol)
            perturber(i)=PerM%chmol
        !---------
        ! Obtain Relaxation matrix elements for each line.
        !
            if (econtrol % e(1) .ge. 1) write(*,*) '>Building Relaxation Matrix...'
            CALL WelCAL(pd1, dta_size1, molP, PerM, Wper, econtrol)
        !---------
        ! Adding the corresponding perturber-molecule 
        ! contribution to the relaxation matrix.
        !
            call add2Wfinal(dta_size1,Wmat,Wper,xMolp(i))
        ENDDO
    !---------
        if (econtrol % e(1) .ge. 1) write(*,*)"<----------------------Finished loop"
        !
        ! Uncomment the following command to print RELAXAION MATRIX ELEMENTS to the screen:
        !
        !CALL show_W(dta_size1,Wmat,int(dta1%J(:,1),8))
        !---------
        ! Renormalization of the Relaxation matrix.
        !
        !
        if (econtrol % e(1) .ge. 1) write(*,*) 'Renormalization procedure of the RM...'
        if (rule2(Ptot,dta_size1,Wmat,dta1%Sig(1:nLines)) .or. (molP%M .eq. 7) ) then
            allocate(Wrno(dta_size1,dta_size1),STAT = IERR3)
            if (IERR3 .ne. 0) call memoError(" Wrno ",econtrol)
            CALL RN_Wmat(dta_size1, pd1, Wmat, Wrno, T, Ptot, econtrol)
            CALL includeW(nLines,vLines_Indx,W_rn, &
                          artsNA, artsGA, rT,Ptot, &
                          dta_size1, Wrno)
            !
            ! ---------
            ! Allocate and Copy RM-data to final struct and file
            !
            if (my_print) then
                CALL alloRMF(nLines, dta2)
                pd2 => dta2
            !
                if (econtrol % e(1) .ge. 1) write(*,*) 'Copying data to final struct...'
                !!CALL W2dta2(nLines, pd1, pd2, W_rn) 
                CALL W2dta2(dta_size1, pd1, pd2, Wrno) 
            !--------
            ! Write RMF file
            ! 
                if (econtrol % e(1) .ge. 1) write(*,*) 'Saving Relaxation Matrix File'
                call save_W2plot(nLines, pd1, pd2, molP, npert, perturber, econtrol, 'htm')
            !
                NULLIFY( pd2 )
            endif
            !
            !
            CALL InitM(nLines,1, Y1)
            CALL InitM(nLines,1, Y2)
            CALL InitM(nLines,1, Y3)
            if (ordered .ge. 1) then
            !---------
            ! Linemixing first order coeff. calculation.
            !
                allocate(Y_RosT(dta_size1),STAT = IERR4)
                if (IERR4 .ne. 0) call memoError("Y_RosT",econtrol)
                !
                if (econtrol % e(1) .ge. 1) write(*,*) 'Linemixing first order coeff...'
                call LM_Rosen(molP,dta_size1,pd1,Wrno,Y_RosT)
                call includeY(nLines,vLines_Indx,Y1,dta_size1,Y_RosT)

                if (ordered .eq. 2) then
                !---------
                ! Linemixing second order coeff. calculation.
                !
                    allocate(Y_G(dta_size1),STAT = IERR4)
                    if (IERR4 .ne. 0) call memoError("Y2   :",econtrol)
                    allocate(Y_DV(dta_size1),STAT = IERR4)
                    if (IERR4 .ne. 0) call memoError("Y3   :",econtrol)
                    !
                    if (econtrol % e(1) .ge. 1) write(*,*) 'Linemixing second order coeffs...'
                    call LM_2ord(molP,dta_size1,pd1,Wrno,Y_G,Y_DV)
                    !call show_PD(nLines,dta1%sig,Y_G,Y_DV)
                    CALL includeY(nLines,vLines_Indx,Y2,dta_size1,Y_G)
                    CALL includeY(nLines,vLines_Indx,Y3,dta_size1,Y_DV)
                    !call show_PD(nLines,artsWNO,Y2,Y3)
                endif
                !--------
                ! Write Y parameter file
                !
                if (my_print) then
                    if (econtrol % e(1) .ge. 1) write(*,*) 'Saving Rosenkranz parameter Y...'
                    CALL save_Yrp(pd1, dta_size1, molP, Y_RosT,'htm')
                endif
            !
            endif
        
        else
            if (econtrol % e(1) .ge. 1) write(*,*) "Rule 2 failed, RM(diagonal matrix) no OFF-diagonal elements are returned."
            econtrol%solu = 0
            CALL just_fill_DiagWRn(nLines,artsNA, artsGA, rT, Ptot,W_rn)
            CALL InitM(nLines,1, Y1)
            CALL InitM(nLines,1, Y2)
            CALL InitM(nLines,1, Y3)
        endif
        
    else
        if (econtrol % e(1) .ge. 1) then
            write(*,*) "Rule 1: Not enough Lines to calculate Relaxation Matrix"
            write(*,*) "        Diagonal matrix sent back in return."
        endif
        econtrol%solu=0
        !
        CALL just_fill_DiagWRn(nLines,artsNA, artsGA, rT, Ptot,W_rn)
        CALL InitM(nLines,1, Y1)
        CALL InitM(nLines,1, Y2)
        CALL InitM(nLines,1, Y3)
        !
        ! Uncomment the following command to print RELAXAION MATRIX ELEMENTS to the screen:
        !CALL show_W(nLines,W_rn,int(dta1%J(:,1),8))
        !---------
    endif
    !
    !
    ! Uncomment the next 3lines to save the Q-basis rates in an ASCII file:
    !aupp=artsUpp(1,1)
    !alow=artsLow(1,1)
    !call save_Q(dta1, dta_size1, molP, 'htm', aupp, alow)
    !
    NULLIFY( pd1 )
!
    if (econtrol % e(1) .ge. 1) PRINT *, "END OF RELMAT SUBROUTINE"
!
    if (econtrol % e(2) .ge. 1) then
        runE_deb = 1
    else
        runE_deb = 0
    endif
!   
!    STOP
!
  END SUBROUTINE RM_LM_tmc_arts
!--------------------------------------------------------------------------------------------------------------------
SUBROUTINE RM_LM_LLS_tmc_arts(nLines, sgmin, sgmax, &
                          artsM, artsI, artsWNO, &
                          artsS, artsGA, artsE00, &
                          artsNA, artsUpp, artsLow, &
                          artsg0 , artsg00, &
                          T, Ptot, QT, QT0, mass, &
                          npert, pert, i_pert, p_mass, p_vmr,&
                          runE_deb,ordered,&
                          W_rn, dipo, rho, &
                          Y1,Y2,Y3) bind(C, name='arts_relmat_interface__linear_type')
!--------------------------------------------------------------------------------------------------------------------
!
! This SUBROUTINE is used to compute the following variables:
!   *Dipole and e-level population      (Spectroscopy Data)
!   *Relaxation Matrix for a given T    (Relaxation Matrix)
!   *Rosenkranz parameter for a 
!              given pair (P,T)         (Rosenkranz's Para)
!   *Second order linemixing P. 
!              on given (P,T)           (2nd order Linemix)
!   on LINEAR MOLECULES.
!
!   NOTE: Check variables up in the interface.
!   ----
!
!   Accessed Files:  'none'
!   ---------------
!
!   Called Routines: 'VarInit'  (VARiable INITialization)
!   ---------------  'MoleculeID' (Molecule IDentification)
!                    'Hit2DTA' (Read HITRAN12 Linear-bands file)
!                    'PopuCal'  (POPUlation CALculation)
!                    'DipCAL'   (DIPole elements CALculation)
!                    'systemQParam_LLS' (External Parameters of the system)
!                    'WelCAL'   (W elements CALculation)
!                    'RN_Wmat' (ReNormalization of W matrix)
!
!
!   T. Mendaza last change 20 Feb 2017
!-------------------------------------------------------------------
!
! MODULES IN USE:
    use module_common_var
    use module_error
        use module_maths
        use module_molecSp
        use module_read
           use module_phsub
            use module_LLS
            use module_linemixing

    Implicit none
!   INPUT variables
    integer*8, intent(in) :: nLines, npert
    integer*8, intent(in) :: artsM, artsI
    integer*8, intent(in) :: artsg0(nLines), artsg00(nLines)
    integer*8, intent(in) :: artsLow(4,nLines), artsUpp(4,nLines)
    integer*8, intent(in) :: pert(npert), i_pert(npert)
    integer*8, intent(in) :: ordered
    integer*8, intent(inout) :: runE_deb
    Double Precision, intent(in)  :: sgmin, sgmax, T, Ptot
    Double Precision, intent(in)  :: QT, QT0, mass
    Double Precision, intent(in)  :: p_vmr(npert), p_mass(npert)
    Double Precision, intent(in)  :: artsWNO(nLines),artsS(nLines), &
                                     artsGA(nLines), artsE00(nLines), &
                                     artsNA(nLines)
!
!   OUTPUT variables
    Double Precision, intent(out) :: rho(nLines), dipo(nLines)
    Double Precision, intent(out) :: Y1(nLines),Y2(nLines),Y3(nLines)  
    Double Precision, intent(out) :: W_rn(nLines,nLines)
!   OTHER VARIABLES
    integer*8              :: dta_size1, IERR1, IERR2, IERR3, IERR4
    integer*8              :: iLine, i, j, k
    integer*8              :: vLines_Indx(nLines) 
    Double Precision, ALLOCATABLE :: Wmat(:,:),&
                                     Wper(:,:),&
                                     Wrno(:,:)
    Double Precision, ALLOCATABLE :: Y_RosT(:),Y_G(:),Y_DV(:)
    Double Precision       :: xMOLp(npert)
    Double Precision       :: faH, rT, deltaV
    Double Precision       :: maxInten, maxGA, maxE00, maxNA
    type (dta_SDF), target :: dta1
    type (dta_SDF), pointer:: pd1
    type (dta_RMF), target :: dta2
    type (dta_RMF), pointer:: pd2
    type (dta_MOL)         :: molP
    type (dta_MOL)         :: PerM
    type (dta_ERR)         :: econtrol
    integer*8              :: sys(2), lM, lP, np,nq,nr, aupp,alow
    character*6            :: perturber(2)
    character*2            :: fmt1, fmt2
    logical                :: enough_Lines, my_print
!   
! ROUTINE:
!
!----------
! Inital common temperature-dependent constants
    rT = T/T0
!
    if (runE_deb .ge. 1) then 
        write(*,*)  "RELMAT RUN-TYPE = Verbose."
        write(*,*)  "Type: Mendaza"
        write(*,2016) T, Ptot
    endif
2016 Format("Starting Linemixing Relaxation Matrix software. T=",f5.0,"K; P=",f7.2," atm")
!----------
! Allocation
    if (runE_deb .ge. 1) write(*,'(a33,i4,a1)') 'Allocate SDF and RMF variables to arrays (', nLines,')'
    CALL alloSDF(nLines, dta1)
    pd1 => dta1
!
!----------
! Band quantities specification
    if (runE_deb .ge. 1) write(*,*) 'Init. Variables...'
    CALL VarInit(molP,econtrol,runE_deb)
    molP % Temp = T !Kelvin
    molP % Ptot = Ptot!atm
    molP % QTy  = "TMC"
    molP % LLSty= "Linear"
    !molP % LLSty= "Model1"
    !molP % LLSty= "Model2"
    !molP % LLSty= "Model3"
    !molP % LLSty= "Model4"
    !molP % LLSty= "Li--AF"

    molP % v0   = meanV0(nLines,artsWNO)
    !
    !Disable the Adiabatic factor option because it requires from external inputs
    molP % AF_ON= .false.
    my_print = .false. 
!
!----------
! Obtainig the molecule ID from the Formula specified in "module_common_var"
    if (econtrol % e(1) .ge. 1) write(*,*) 'Identifying molecule and loading its parameters...'
    !
    CALL moleculeID(artsM, artsI, mass, QT, QT0, .true., molP,&
                    econtrol)
    dta1%M = molP%M
!---------
! Call for reading HITRAN spectroscopy data.
!
    if (econtrol % e(1) .ge. 1) write(*,*) 'Locating Band information...'
    !
    CALL Hit2DTA(pd1, dta_size1, nLines, vLines_Indx, &
                 artsWNO, artsS  , artsGA , artsE00, &
                 artsNA , artsUpp, artsLow, &
                 artsg0 , artsg00, sgmin  , sgmax, &
                 econtrol)
!---------
! Compute the relative population of the lower state
! at Temperature T0
    call PopuCAL(pd1,dta_size1, molP, econtrol)
    !
    do j = 1, nLines
        if (vLines_Indx(j) .eq. 0) then
            if (T .eq. T0) then
                rho(j) = artsg00(j)*dexp(-c2*artsE00(j)/T0)/QT0
            else
                rho(j) = artsg00(j)*dexp(-c2*artsE00(j)/T0)* &
                         dexp(-c2*artsE00(k)*(1.d0/T-1.d0/T0))/QT
            endif
        else
            if (T .eq. T0) then
                rho(j) = dta1 % PopuT0(vLines_Indx(j))
            else
                rho(j) = dta1 % PopuT(vLines_Indx(j))
            endif
        endif 
    enddo
! NOTE: we use 'tra' mode (see 'PopuT0') because we are producing 
! Input files to Ha Tran Line-mixing code.
!---------
! Calculate Dipole element for each line.
!
    if (econtrol % e(1) .ge. 1) write(*,*) 'Calculating Dipole moment'
    CALL DipCAL(pd1,dta_size1,molP,econtrol)
!
    do j = 1, nLines
        if (vLines_Indx(j) .eq. 0) then
            dipo(j) = dsqrt(artsS(k)/(artsWNO(k)* &
                        rho(k)*(1.D0 - &
                        dexp(-c2*artsWNO(k)/T0))))
        else 
            dipo(j) = dta1%DipoT(vLines_Indx(j)) 
        endif 
    enddo
!
! Uncomment the following command to print POPULATION & DIPOLE to the screen:
!
    !call show_PD(nLines, dta1%Sig(1:dta_size1), dta1 % PopuT0(1:dta_size1), dipo)
    !call show_PD(nLines, dta1%Sig(1:dta_size1), dta1 % PopuT0(1:dta_size1), dta1 % D0(1:dta_size1))
!---------
! Write SDF file
!
    if (rule1(dta_size1) .and. .not.(ordered .eq. 0)) then
        if (econtrol % e(1) .ge. 1) then
            write(*,*)"Looping over system of perturbers..."
            write(*,*)"----------------------------------->"
        endif
    !
    ! Allocate variables according to the valid number of lines:
        allocate(Wmat(dta_size1,dta_size1),STAT = IERR1)
        IF (IERR1 .ne. 0) call memoError(" Wmat ",econtrol)
        allocate(Wper(dta_size1,dta_size1),STAT = IERR2)
        IF (IERR2 .ne. 0) call memoError(" Wper ",econtrol)
    !
    !   Init final NO-RN W (relaxation matrix):
        CALL InitM(dta_size1,dta_size1,Wmat)
    !   Looping:
        DO i = 1,npert
            if ( (i .eq. 1) .or. (molP%AF_ON) ) then
            !
            ! Identifying Perturbers Molecule and
            ! its ATM concentration (value must be from 0-1):
                if (molP%AF_ON) then
                    xMOLp(i) = (p_vmr(i)/sum(p_vmr))
                else
                    xMOLp(i) = 1.0_dp
                endif
            !
            !   Init W rel-mat/perturber and collider structure:
                CALL mol_Init(PerM)
            !---------
            ! Identifying perturber molecule.
            !
            ! Perturber Molecule: This molecule has to be 
            ! significatively faster than the molecule in study.
                if (econtrol % e(1) .ge. 1) write(*,*) '>Identifying perturber molecule...'
                CALL moleculeID(pert(i), i_pert(i), p_mass(i), &
                                0.0_dp, 0.0_dp, .false., &
                                PerM, econtrol)
            !----------
            ! let's take th proper a1, a2, a3, dc adjust parameters for the
            ! system:
                if ( molP%AF_ON .and. (econtrol % e(1) .ge. 1)) then
                    write(*,*)  ">>System: {",trim(molP%chmol)," - ",trim(PerM%chmol),"}"
                    ! sys = "CO2-N2" or "CO2-O2"
                endif
            !
                sys(1) = molP%M
                sys(2) = PerM%M
                CALL systemQParam_LLS(sys,molP)
            !
                if (.not.(molP%availableParam)) then
                    if (molP % LLSty .eq. "Li--AF") then
                        CALL calc_QParam_AF(dta_size1, pd1, molP, PerM, econtrol)
                        if (econtrol % e(1) .ge. 1) write(*,*) "A1 = ", molP%a1,";A2= ", molP%a2,";A3= ", molP%a3
                        if (econtrol % e(1) .ge. 1) write(*,*) "A4 = ", molP%a4,";A5= ", molP%a5,";A6= ", molP%a6
                        if (econtrol % e(1) .ge. 1) write(*,*) "A7 = ", molP%a7,";A8= ", molP%a8,";A9= ", molP%a9
                    else !if "Linear" or "Model1/2/3/4"
                        !CALL calc_QParam(dta_size1, pd1, molP, PerM, econtrol)
                        !CALL calc_QPar_DGELSY(dta_size1, pd1, molP, PerM, econtrol)
                        !CALL calc_QPar_DGELSS(dta_size1, pd1, molP, PerM, econtrol)
                        CALL calc_QPar_DGELSD(dta_size1, pd1, molP, PerM, econtrol)
                        if (econtrol % e(1) .ge. 1) write(*,*) "A1 = ", molP%a1,";A2= ", molP%a2,";A3= ", molP%a3
                    endif
                endif
            !
                perturber(i)=PerM%chmol
            !---------
            ! Obtain Relaxation matrix elements for each line.
            !
                if (econtrol % e(1) .ge. 1) write(*,*) '>Building Relaxation Matrix...'
                CALL WelCAL(pd1, dta_size1, molP, PerM, Wper, econtrol)
            !---------
            ! Adding the corresponding perturber-molecule 
            ! contribution to the relaxation matrix.
            !
                CALL add2Wfinal(dta_size1,Wmat,Wper,xMolp(i))
            endif
        ENDDO
    !---------
        if (econtrol % e(1) .ge. 1) write(*,*)"<----------------------------Finished loop"
    !
    ! Uncomment the following command to print RELAXAION MATRIX ELEMENTS to the screen:
    !
        !CALL show_W(nLines,Wmat,int(pd1%J(:,1),8))
    !---------
    ! Renormalization of the Relaxation matrix.
    !
        if (econtrol % e(1) .ge. 1) write(*,*) 'Renormalization procedure of the RM...'
        allocate(Wrno(dta_size1,dta_size1),STAT = IERR3)
        if (IERR3 .ne. 0) call memoError(" Wrno ",econtrol)
        CALL RN_Wmat(dta_size1, pd1, Wmat, Wrno, T, Ptot, econtrol)

        CALL includeW(nLines,vLines_Indx,W_rn, &
                          artsNA, artsGA, rT, Ptot,&
                          dta_size1, Wrno)       
    ! ---------
    ! Allocate and Copy RM-data to final struct and file
    !
        if (my_print) then
            CALL alloRMF(nLines, dta2)
            pd2 => dta2
    !
            if (econtrol % e(1) .ge. 1) write(*,*) 'Copying data to final struct...'
            !!CALL W2dta2(nLines, pd1, pd2, W_rn) 
            CALL W2dta2(dta_size1, pd1, pd2, Wmat) 
    !--------
    ! Write RMF file
    ! 
            if (econtrol % e(1) .ge. 1) write(*,*) 'Saving Relaxation Matrix File'
            call save_W2plot(nLines, pd1, pd2, molP, npert, perturber, econtrol, 'tmc')
    !
            NULLIFY( pd2 )
        endif
    !
        CALL InitM(nLines,1, Y1)
        CALL InitM(nLines,1, Y2)
        CALL InitM(nLines,1, Y3)
        if (ordered .ge. 1) then
        !---------
        ! Linemixing first order coeff. calculation.
        !
            allocate(Y_RosT(dta_size1),STAT = IERR4)
            if (IERR4 .ne. 0) call memoError("Y_RosT",econtrol)
            !
            if (econtrol % e(1) .ge. 1) write(*,*) 'Linemixing first order coeff...'
            call LM_Rosen(molP,dta_size1,pd1,Wrno,Y_RosT)
            CALL includeY(nLines,vLines_Indx,Y1,dta_size1,Y_RosT)

            if (ordered .eq. 2) then
            !---------
            ! Linemixing second order coeff. calculation.
            !
                allocate(Y_G(dta_size1),STAT = IERR4)
                if (IERR4 .ne. 0) call memoError("Y2   :",econtrol)
                allocate(Y_DV(dta_size1),STAT = IERR4)
                if (IERR4 .ne. 0) call memoError("Y3   :",econtrol)
            !
                if (econtrol % e(1) .ge. 1) write(*,*) 'Linemixing second order coeffs...'
                call LM_2ord(molP,dta_size1,pd1,Wrno,Y_G,Y_DV)
                CALL includeY(nLines,vLines_Indx,Y2,dta_size1,Y_G)
                CALL includeY(nLines,vLines_Indx,Y3,dta_size1,Y_DV)
            endif
            ! --------
            ! Write Y parameter file
            !
            if (my_print) then
                if (econtrol % e(1) .ge. 1) write(*,*) 'Saving Rosenkranz parameter Y...'
                !CALL save_Yrp(pd1, dta_size1, molP, Y_RosT,'tmc')
                ! ERASE --->
                CALL save_Yrp(pd1, dta_size1, molP, Y_RosT,'tm2')
                !<--- this line and change my_print!!!
            endif
        !
        endif
    else
        if (econtrol % e(1) .ge. 1) then
            write(*,*) "Rule 1: Not enough Lines to calculate Relaxation Matrix"
            write(*,*) "        Diagonal matrix sent back in return."
            econtrol%solu = 0
        endif
        !
        !
        CALL just_fill_DiagWRn(nLines,artsNA, artsGA, rT, Ptot, W_rn)
        CALL InitM(nLines,1, Y1)
        CALL InitM(nLines,1, Y2)
        CALL InitM(nLines,1, Y3)
        !
        ! Uncomment the following command to print RELAXAION MATRIX ELEMENTS to the screen:
        !
        !CALL show_W(nLines,W_rn,int(dta1%J(:,1),8))
        !---------
    endif
    !
    np = 0
    nq = 0
    nr = 0
    do j = 1, nLines
        if (dta1%br(j) .eq. "p") then
            np = np+1
        else if (dta1%br(j) .eq. "q") then
            nq = nq+1 
        else if (dta1%br(j) .eq. "r") then
            nr = nr+1
        endif 
    enddo
    !
    !
    !    maxInten = maxf(nLines, artsS)
    !    deltaV = dabs(maxf(nLines, artsWNO) - minf(nLines, artsWNO))
    !    maxGA = maxf(nLines,artsGA) 
    !    maxE00 = maxf(nLines,artsE00)
    !    maxNA = maxf(nLines,artsNA)
    !    aupp=artsUpp(1,1)
    !    alow=artsLow(1,1)     
    !    call save_Q(dta1, dta_size1, molP, 'tmc', aupp, alow)
    !
    !
    NULLIFY( pd1 )
!
!
    if (econtrol % e(1) .ge. 1) PRINT *, "END OF RELMAT SUBROUTINE"
!
    if (econtrol % e(2) .ge. 1) then
        runE_deb = 1
    else
        runE_deb = 0
    endif
    !STOP
!
  END SUBROUTINE RM_LM_LLS_tmc_arts
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE alloSDF(n,dta1)
!--------------------------------------------------------------------------------------------------------------------
! Allocate the right block of memory for the selected band, SDF type.
! dta1    : dta type "dta_SDF". Spectrocopic parameters.
!
    use module_common_var
    implicit none
    integer*8, intent(in )             :: n 
    type (dta_SDF),intent(out)         :: dta1

!----------
! 
    ! Init. SDF data
      dta1%M   = 0
      dta1%iso = 0  
      dta1%lv2 = (/0,0/)

      ALLOCATE(dta1%J(n,2),dta1%N(n,2),dta1%nspin(n,2),dta1%espin(n,2))
      ALLOCATE(dta1%Sig(n),dta1%Str(n),dta1%E(n),dta1%HWT0(n),dta1%BHW(n),&
               dta1%SHIFT(n),dta1%swei0(n),dta1%swei00(n))
      ALLOCATE(dta1%PopuT0(n),dta1%PopuT(n))
      ALLOCATE(dta1%D0(n),dta1%DipoT(n))
      !ALLOCATE(dta1%Drigrotor(n),dta1%DipoT0(n))
      ALLOCATE(dta1%F(n,2),dta1%br(n),dta1%br_N(n))
      !ALLOCATE XXXXXX
      ALLOCATE(dta1%Qlt(n,500))

  END SUBROUTINE alloSDF
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE alloRMF(n,dta2)
!--------------------------------------------------------------------------------------------------------------------
! Allocate the right block of memory for the selected band, RMF type.
! dta2    : dta type "dta_RMF". Relaxation Matrix parameters.

!
    use module_common_var
    implicit none
    integer*8     , intent(in )        :: n 
    type (dta_RMF), intent(out)        :: dta2

!----------
! 
    ALLOCATE (dta2%WT0(n*n))

  END SUBROUTINE alloRMF
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE VarInit(molP,econ,runE)
!--------------------------------------------------------------------------------------------------------------------
! "VarInit": Variables initialization
! 
! Detailed Description:
! ---------------------
! This subroutine starts every variable type in this program and set it to zero. 
!
! Variables:
!
! Input/Output Parameters of Routine (Arguments or Common)
! ----------------------------------
! molP    : dta type "dta_MOL". Molecular structure parameters.
! econ    : dta type "dta_ERR". Error structure parameters.
! runE    : type of run selected by the user.
!
! Accessed Files:  None
! --------------
!
! Called Routines: "mol_Init"
! ---------------  
!
! Called By: "RM_LM_LLS_tmc_arts" and "RM_LM_tmc_arts"
! ---------
!
!
! T.Mendaza, last change 17 February 2017
!--------------------------------------------------------------------------------------------------------------------
!
    use module_common_var
    implicit none
    integer*8     , intent(in)       :: runE
    type (dta_MOL), intent(inout)    :: molP
    type (dta_ERR), intent(inout)    :: econ
    integer*8 :: i,j,k

!----------
! Init. MOLECULE data
      call mol_Init(molP)
! 
!
! ERROR CONTROL:
            econ % e(1) = runE
            econ % e(2) = 0
            econ % solu = 1
      Return
  END SUBROUTINE VarInit
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE mol_Init(molP)
!--------------------------------------------------------------------------------------------------------------------
! "mol_Init": Molecular-structure type initialization
!
! Variables:
!
! Input/Output Parameters of Routine (Arguments or Common)
! ----------------------------------
! molP    : dta type "dta_MOL". Molecular structure parameters.
!
! Accessed Files:  None
! --------------
!
! Called Routines: None
! ---------------  
!
! Called By: "VarInit"
! ---------
!
! T.Mendaza, last change 16 February 2017
!--------------------------------------------------------------------------------------------------------------------
!
    use module_common_var
    implicit none
    type (dta_MOL), intent(inout)    :: molP
!----------
! Init. MOLECULE data
      !Integer kind
      molP%M     = 0
      molP%iso_m = 0
      molP%Aco   = 0 
      !Double precision
      molP % mms  = 0.0_dp
      molP % Temp = 0.0_dp
      molP % Ptot = 0.0_dp
      molP % Nmcon= 0.0_dp
      molP % B0   = 0.0_dp
      molP % QT   = 0.0_dp
      molP % QT0  = 0.0_dp
      molP % a1   = 0.0_dp; molP % a2 = 0.0_dp; molP % a3 = 0.0_dp
      molP % a4   = 0.0_dp; molP % a5 = 0.0_dp; molP % a6 = 0.0_dp
      molP % a7   = 0.0_dp; molP % a8 = 0.0_dp; molP % a9 = 0.0_dp
      molP % a10  = 0.0_dp; molP % a11= 0.0_dp
      molP % dc   = 0.0_dp
      molP % ex1  = 0.0_dp
      molP % ex2  = 0.0_dp
      !character
      molP%chmol = ""
      molP%QTy   = "REG"
      molP%LLSty = "Linear"
      !logical
      molP%availableParam = .true.
      molP%AF_ON = .true.
! 
      Return
  END SUBROUTINE mol_Init
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE InitM(n,m,W)
!--------------------------------------------------------------------------------------------------------------------
! 
!
    implicit none
    integer*8, intent(in)               :: n,m 
    Double Precision, intent(out)       :: W(n,m)
    integer*8                           :: i,j,k

!----------
! 
    DO i = 1, n
        DO j = 1,m 
            W(i,j) = 0.0
        ENDDO
    ENDDO

  END SUBROUTINE InitM
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE includeW(n,indx,W_rn, NA, GA, rTT0, P, n0, Wrno)
!--------------------------------------------------------------------------------------------------------------------
! 
    implicit none
    integer*8, intent(in )              :: n, n0 
    integer*8, intent(in )              :: indx(n) 
    Double Precision, intent(in)        :: NA(n), GA(n), rTT0, P
    Double Precision, intent(in )       :: Wrno(n0,n0)
    Double Precision, intent(out)       :: W_rn(n,n)
    integer*8                           :: i,j,k
    real*8                              :: faH

!----------
! 
    DO i = 1, n
        DO j = 1,n
            if ((indx(i).eq.0) .OR. (indx(j).eq.0)) then
                ! Diagonal levels of the Matrix
                if(i .eq. j) then
                    faH = rTT0**NA(j) 
                    W_rn(j,j) = 2*P*GA(j)*faH !+ i*(-0.008)
                else
                    W_rn(i,j) = 0.0d0
                endif 
            else
                W_rn(i,j) = Wrno(indx(i),indx(j))
            endif
        ENDDO
    ENDDO

  END SUBROUTINE includeW
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE includeY(n,indx,Yf, n0, Yc)
!--------------------------------------------------------------------------------------------------------------------
! 
    implicit none
    integer*8, intent(in )              :: n !number of lines input from ARTS 
    integer*8, intent(in )              :: n0!number of lines with proper quantum numbers
    integer*8, intent(in )              :: indx(n) 
    Double Precision, intent(in )       :: Yc(n0)
    Double Precision, intent(out)       :: Yf(n)
    integer*8                           :: i,j,k

!----------
! 
    DO i = 1, n
        if (indx(i).eq.0) then
            Yf(i) = 0.0d0
        else
            Yf(i) = Yc(indx(i))
        endif
    ENDDO

  END SUBROUTINE includeY
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE just_fill_DiagWRn(n, NA, GA, rTT0, P, Wrn)
!--------------------------------------------------------------------------------------------------------------------
! 
    implicit none
    integer*8, intent(in )              :: n 
    Double Precision, intent(in)        :: NA(n), GA(n), rTT0, P
    Double Precision, intent(out)       :: Wrn(n,n)
    integer*8                           :: i,j,k
    real*8                              :: faH

!----------
! 
    CALL InitM(n,n,Wrn)
    !
    do j =1, n
        faH = (rTT0**NA(j)) 
        Wrn(j,j) = 2*P*GA(j)*faH !+ i*(-0.008)
        ! NOTE: 
        ! in case one would like to include water vapour broadening
        ! check the procedure in "WelCAL" subroutine  
    enddo

  END SUBROUTINE just_fill_DiagWRn
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE add2Wfinal(n,Wfinal,Wadd,xMol)
!--------------------------------------------------------------------------------------------------------------------
! 
!
    implicit none
    integer*8, intent (in)                :: n 
    Double Precision, intent (in)         :: xMOL
    Double Precision, intent (in)         :: Wadd(n,n)
    Double Precision, intent(inout)       :: Wfinal(n,n)
    integer*8                             :: i,j,k

!----------
! 
    DO i = 1, n
        DO j = 1,n 
            Wfinal(i,j) = Wfinal(i,j) + xMOL*Wadd(i,j)
        ENDDO
    ENDDO

  END SUBROUTINE add2Wfinal
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE show_W(n,Wfinal,Ji)
!--------------------------------------------------------------------------------------------------------------------
! 
!
    use module_common_var
    implicit none
    integer*8, intent (in)         :: n, Ji(n)
    Double Precision, intent (in)  :: Wfinal(n,n)
    integer*8                      :: i,j,k

!----------
! 
    DO i = 1, n
        DO j = 1,n 
            write(*,1002) "W(",i,",",j,")=",Wfinal(i,j),";Ji =",Ji(i),";Jip=",Ji(j) 
1002 Format(a2,i3,a1,i3,a3,E12.3,a5,i3,a5,i3)
        ENDDO
    ENDDO

    STOP

  END SUBROUTINE show_W
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE show_PD(n,wno,p,d)
!--------------------------------------------------------------------------------------------------------------------
! 
!
    use module_common_var
    implicit none
    integer*8, intent (in)          :: n 
    Double Precision, intent (in)   :: wno(n), p(n), d(n)
    integer*8                       :: i,j,k

!----------
    write (*, *) "Sig,         PopuT0,    DipoT "
    DO j = 1, n
        write (*, 1003) wno(j), p(j), d(j)
1003 format(f12.6,2x,e9.3,2x,e9.3,2x,a7,2x,a7)
    ENDDO
    stop
  END SUBROUTINE show_PD
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE save_W2plot(n, dta1, dta2 ,molP, npert, pert,econ, model)
!--------------------------------------------------------------------------------------------------------------------
! 
!
    use module_common_var
    use module_error
    use module_molecSP
    use module_error

    implicit none
    integer*8, intent (in)                  :: n, npert !n=dta_size1
    type (dta_MOL), intent (in)             :: molP
    type (dta_SDF), intent (in)             :: dta1
    type (dta_RMF), intent (in)             :: dta2
    type (dta_ERR), intent (inout)          :: econ
    character (6) , intent (in)             :: pert(npert)
    character (3) , intent (in)             :: model
    !subroutine Variables
    integer*8, parameter                    :: u=10 !file unit
    character*100                           :: path
    character*60                            :: coup_levels
    character*6                             :: cTemp
    character                               :: ai, af
    integer*8                               :: imatrix
    integer*8                               :: i, j, bri, brf
    integer*8                               :: Ji, Jip
    integer*8                               :: today(3)
!
!------> T. Mendaza; last change 30 January 2017
!
! INIT. VAR.
    write(cTemp,'(f5.1)') molP%Temp
!   NO-Renormalized Matrix
    !path = "RMF2plot_"//trim(cTemp(1:3))//"K.dat"
!   ReNormalized Matrix
    !path = "RMF2plot_RN_"//trim(cTemp(1:3))//"K.dat"
    !path ="RMF2plot"//trim(molP%chmol)//"_"//model//"_RN_"//trim(cTemp(1:3))//"K.dat"
    path ="RM_paper_"//trim(molP%chmol)//"_"//model//"_"//trim(cTemp(1:3))//"K.dat"

    call idate(today)   ! today(1)=day, (2)=month, (3)=year
    open (UNIT = u, FILE = trim(path), STATUS = 'REPLACE', ACTION = 'WRITE')
! HEADER
    write (u, 1004) molP%chmol,adjustL(trim(pert(1))),adjustL(trim(pert(2)))!,"N2+O2"
    write (u, *) "Temperature: "//cTemp//"K"
!    write (u, *) "HITRAN Band:", band
    write (u, *) "Based on HITRAN 2012 database"
    write (u, *) "Author: Teresa Mendaza  "
    write (u, 1005) today(3), today(2), today(1)
    write (u, *) "J00i| deltai | J00f | deltaf |W element"
    ! DATA
    ! the length of RMF file (stored in dta2), is equal to
    ! the length of the SDF file to de power of 2.
    ! Of course, RMF is an squared matrix.
    imatrix=1
    do i = 1, n   ! do1
      do j = 1, n ! do2
        !  NOTES: 
        !  (1) Prime and double primes refer, respectively, 
        !  to upper and lower states, respectively, i.e.
        !       * upper := 2; 0 ; or ' 
        !       * lower := 1;00 ; or " 
        !--------------------------------------------------
        !W(i,j) = WT0(imatrix)
        !tr2tr = auxiQupp//auxiQlow//auxjQupp//auxjQlow
        !      Q0:                   Q00:                
        !                     F´          br J" sym" F"    
        !                 10x a5       5x a1 i3 a1   a5
!        read(dta2%tr2tr(imatrix),1006),ai,Ji,af,Jip
        Ji = dta1%J(i,1); Jip = dta1%J(j,1)
        ai = dta1%br(i) ; af  = dta1%br(j)
        write(u, 1007) Ji, branch2delta(ai,i,econ),&
                       Jip, branch2delta(af,j,econ), &
                       dta2%WT0(imatrix)
                      
       if (imatrix .gt. nMmx) then
           call sizeError('1001',imatrix,nMmx, econ)
       endif
       imatrix=imatrix+1
      end do  ! do2
    end do    ! do1
    
    close (u)

1004 format ( ' RMF file transition parameters system: ', a4, '--', a3,'+'a3 )    
1005 format ( ' Last update: ', i4.4, '/', i2.2, '/', i2.2 )
1006 format (20x,a1,i3,26x,a1,i3)
1007 format (i3,1x,i3,1x,i3,1x,i3,1x,E15.7)


  END SUBROUTINE save_W2plot
!--------------------------------------------------------------------------------------------------------------------
  subroutine save_Yrp(dta1, n, molP, Y_RosenP, model)
!--------------------------------------------------------------------------------------------------------------------
!  Write the results in 'Y_air_TTTK.dat'
!  Output_3: (Variables included)
!  ---------
!  wno
!  Str
!  Y
!  Ji
!  Jf
!  NOTE: Eventhough this is a Relaxation matrix it is stored in list-formated
!    ordered as follows:
!  Matrix:
!   a11 a12 a13 ...
!   a21 a22 a23 ...
!   a31 a32 a33 ...
!   ... ... ... ...
!  Raw:
!   a11 a12 a13 ... a21 a22 a23 ... a31 a32 a33 ...
!
!  where each element represents the rotational state-to-state cross sections 
!  within a single vibrational state.
!------> T. Mendaza; last change 30 January 2017
!
    use module_common_var

    implicit none
    integer*8  , intent (in)                :: n
    double precision, intent (in)           :: Y_RosenP(n)
    type (dta_SDF) , intent (in)            :: dta1
    type (dta_MOL), intent (in)             :: molP
    character(3), intent (in)               :: model
    !subroutine Variables
    integer*8, parameter                    :: u=11 !file unit
    character*100                           :: path
    character*6                             :: cTemp
    integer*8                               :: i, j, k
    integer*8                               :: today(3)

! INIT. VAR.
    write(cTemp,'(f5.1)') molP%Temp
    !path = trim(out_file_path)//trim(out_fil2_RMF)
    path ="Y_Test_"//trim(molP%chmol)//"_"//model//"_"//trim(cTemp(1:3))//"K.dat"
    call idate(today)   ! today(1)=day, (2)=month, (3)=year
    open (UNIT = u, FILE = trim(path), STATUS = 'REPLACE', ACTION = 'WRITE')
! HEADER
    write (u, *) "Line-mixing Rosenkranz Parameter file; Molecule:", molP%chmol
!    write (u, *) "HITRAN Band:", trim(dta1 % hitBand)
    write (u, *) "Based on HITRAN 2012 database"
    write (u, *) "Author: Teresa Mendaza  "
    write (u, 1008) today(3), today(2), today(1)
    write (u, *) "|WaveNumber(cm-1) |Intensity cm-1/(molecules·cm-2) | Y Parameter  | Jlow | Jup"
   !
    ! DATA
    do i = 1, n   ! do1
        write (u, 1009) dta1%sig(i),dta1%Str(i),&
                       Y_RosenP(i),&
                       int(dta1%J(i,1)),int(dta1%J(i,2))
    end do    ! do1
    
    close (u)
    
1008 format ( ' Last update: ', i4.4, '/', i2.2, '/', i2.2 )
1009 format (F12.6,8x,E10.3,14x,F20.4,4x,i3,4x,i3)

end subroutine save_Yrp
!--------------------------------------------------------------------------------------------------------------------
  subroutine save_Q(dta1, n, molP, model, aUpp, aLow)
!--------------------------------------------------------------------------------------------------------------------
!  ! ERASE THIS FUNCTION
!  ---------
!  (a(i,j), j=1,numcols)
!  Raw:
!   a11 a12 a13 ... a21 a22 a23 ... a31 a32 a33 ...
!
!  where each element represents the rotational state-to-state cross sections 
!  within a single vibrational state.
!------> T. Mendaza; last change 30 January 2017
!
    use module_common_var

    implicit none
    integer*8  , intent (in)                :: n, aUpp, aLow
    type (dta_SDF) , intent (in)            :: dta1
    type (dta_MOL), intent (in)             :: molP
    character(3), intent (in)               :: model
    !subroutine Variables
    integer*8, parameter                    :: u=14 !file unit
    character*100                           :: path
    character*6                             :: cTemp
    character*2                             :: band
    integer*8                               :: i, j, k
    integer*8                               :: today(3)

! INIT. VAR.
    band=''
    write(cTemp,'(f5.1)') molP%Temp
    write(band,1000) aUpp,aLow
    path ="Q000_"//trim(molP%chmol)//"_"//model//"_"//trim(cTemp(1:3))//"K_"//band//".dat"
    call idate(today)   ! today(1)=day, (2)=month, (3)=year
    open (UNIT = u, FILE = trim(path), STATUS = 'REPLACE', ACTION = 'WRITE')
! HEADER
   !
    ! DATA
    do i = 1, n   ! do1 over transitions
        write (u, *) (dta1%Qlt(i,j), j=1,500)
    end do    ! do1
    
    close (u)
1000 format (2(i1))
1001 format ( ' Last update: ', i4.4, '/', i2.2, '/', i2.2 )

end subroutine save_Q
!--------------------------------------------------------------------------------------------------------------------
