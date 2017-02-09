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
        !   artsGA   : Air pressure broadening in cm-1/atm halfwidth
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

        subroutine show_PD(n,wno,p,d)
            use module_common_var
            implicit none
            integer*8, intent (in)   :: n 
            Double Precision, intent (in)   :: wno(n), p(n), d(n)
        end subroutine show_PD

        subroutine save_W2plot(n, dta1, dta2 ,molP, npert, pert)
            use module_common_var
            use module_molecSP
            use module_error

            implicit none
            integer*8, intent (in)                  :: n, npert !n=dta_size1
            type (dta_MOL), intent (in)             :: molP
            type (dta_SDF), intent (in)             :: dta1
            type (dta_RMF), intent (in)             :: dta2
            character (6)                           :: pert(npert)
        end subroutine save_W2plot

        subroutine save_Yrp(dta1, n, molP, Y_RosenP)
            use module_common_var

            implicit none
            integer*8  , intent (in)                :: n
            double precision, intent (in)           :: Y_RosenP(n)
            type (dta_SDF) , intent (in)            :: dta1
            type (dta_MOL), intent (in)             :: molP
        end subroutine save_Yrp


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
!   artsGA   : Air pressure broadening in cm-1/atm halfwidth
!   arts E00 : Lower state energy in cm-1
!   artsNA   : Pressure broadening constant in cm-1
!   artsUpp  : Upper state quantum numbers. First is L2, then is J, then is N, then is S.  
!              If the quantum number is not applicable, the position contains the number -1.
!   artsLow  : Same as artsUpp but for lower state numbers.
!   artsg0   : Upper state g-constant
!   artsg00  : Lower state g-constant
!   T        : Temperature of the system in Kelvin.
!   Ptot     : Total pressure in atm.
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
    Double Precision       :: Wmat(nLines,nLines),&
                              Wper(nLines,nLines)
    Double Precision       :: Y_RosT(nLines)
    Double Precision       :: xMOLp(npert)
    type (dta_SDF)         :: dta1
    type (dta_RMF)         :: dta2
    type (dta_MOL)         :: molP
    type (dta_MOL)         :: PerM
    integer*8              :: sys(2), lM, lP
    character*6            :: perturber(2)
    character*2            :: fmt1, fmt2
    logical                :: enough_Lines
!
    !Ptot = Ptot*9.869200E-06
    write(*,2016), T, Ptot
2016 Format("Starting Linemixing Relaxation Matrix software. T=",f5.0,"K; P=",f7.2," atm")
!----------
! Band quantities specification
    print*, 'Init. Variables...'
    CALL VarInit(dta1,dta2,molP)
    molP % Temp = T !Kelvin
    molP % Ptot = Ptot!atm
!
!----------
! Obtainig the molecule ID from the Formula specified in "module_common_var"
    print*, 'Identifying molecule and loading its parameters...'
    !
    CALL moleculeID(artsM, artsI, mass, QT, QT0, .true., molP)
    dta1%M = molP%M
!---------
! Call for reading HITRAN spectroscopy data.
!
    print*, 'Locating Band information...'
    !
    CALL Hit2DTA(dta1, dta_size1, nLines, enough_Lines, &
                                            artsWNO, &
                                            artsS, &
                                            artsGA, &
                                            artsE00, &
                                            artsNA , &
                                            artsUpp, artsLow, &
                                            artsg0 , artsg00, &
                                            sgmin  , sgmax)
    nLines = dta_size1 
    dta_size2=nLines**2 
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
!
    do j = 1, nLines
!
!   DATA
        dipo(j) = dta1%DipoT(j)
    enddo
!
! Uncomment the following command to print POPULATION & DIPOLE to the screen:
!
    !call show_PD(nLines, dta1%Sig(1:nLines), dta1 % PopuT0(1:nLines), dipo)
    !call show_PD(nLines, dta1%Sig(1:nLines), dta1 % PopuT0(1:nLines), dta1 % D0(1:nLines))
!---------
! Write SDF file
!
    print*,"Looping over system of perturbers..."
    print*,"--------------------------------->"

!
!   Init W rel-mat.
    CALL InitW(nLines,Wmat)
!   Looping:
    DO i = 1,npert
        !
        ! Identifying Perturbers Molecule and
        ! its ATM concentration (value must be from 0-1):
        xMOLp(i) = (p_vmr(i)/sum(p_vmr))
        !
        CALL InitW(nLines,Wper)
        !---------
        ! Identifying perturber molecule.
        !
        ! Perturber Molecule: This molecule has to be 
        ! significatively faster than the molecule in study.
        print*, '>Identifying perturber molecule...'
        CALL moleculeID(pert(i), i_pert(i), p_mass(i), 0.0_dp, 0.0_dp, .false. ,PerM)
        !----------
        ! let's take th proper a1, a2, a3, dc adjust parameters for the
        ! system:
        write(*,*), ">>System: {",trim(molP%chmol)," - ",trim(PerM%chmol),"}"
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
!
! Uncomment the following command to print RELAXAION MATRIX ELEMENTS to the screen:
!
    !CALL show_W(nLines,Wmat,int(dta1%J(:,1),8))
    !stop
!---------
! Renormalization of the Relaxation matrix.
!
    CALL InitW(nLines,W_rn)
    print*, 'Renormalization procedure of the RM...'
    CALL RN_Wmat(nLines, dta1, Wmat, W_rn) 
!
!
!    print*, 'Copying data to final struct...'
!    CALL W2dta2(nLines, dta1, dta2, W_rn) 
!    CALL W2dta2(nLines, dta1, dta2, Wmat) 
!--------
! Write RMF file
! 
!    print*, 'Saving Relaxation Matrix File'
!    call save_W2plot(nLines, dta1, dta2, molP, npert, perturber)
!
!---------
! Linemixing first order coeff. calculation.
!
!    print*, 'Linemixing first order coeff...'
!    call LM_Rosen(molP,nLines,dta1,W_rn,Y_RosT)
!--------
! Write Y parameter file
!
!    print*, 'Saving Rosenkranz parameter Y...'
!    CALL save_Yrp(dta1, nLines, molP, Y_RosT)
!
    PRINT *, "Successful run! :)"
!   
!    STOP
!
  END SUBROUTINE RM_LM_tmc_arts
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE InitW(n,W)
!--------------------------------------------------------------------------------------------------------------------
! 
!
    implicit none
    integer*8, intent(in )              :: n !== Molecules isotope(AFGL code)
    Double Precision, intent(out)       :: W(n,n)
    integer*8                           :: i,j,k

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
            write(*,1000), "W(",i,",",j,")=",Wfinal(i,j),";Ji =",Ji(i),";Jip=",Ji(j) 
1000 Format(a2,i3,a1,i3,a3,E12.3,a5,i3,a5,i3)
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
        write (*, 1004) wno(j), p(j), d(j)
1004 format(f12.6,2x,e9.3,2x,e9.3,2x,a7,2x,a7)
    ENDDO
    stop
  END SUBROUTINE show_PD
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE save_W2plot(n, dta1, dta2 ,molP, npert, pert)
!--------------------------------------------------------------------------------------------------------------------
! 
!
    use module_common_var
    use module_molecSP
    use module_error

    implicit none
    integer*8, intent (in)                  :: n, npert !n=dta_size1
    type (dta_MOL), intent (in)             :: molP
    type (dta_SDF), intent (in)             :: dta1
    type (dta_RMF) , intent (in)            :: dta2
    character (6)                           :: pert(npert)
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
    write(cTemp,'(f5.1)'),molP%Temp
!   NO-Renormalized Matrix
    !path = "RMF2plot_"//trim(cTemp(1:3))//"K.dat"
!   ReNormalized Matrix
    path = "RMF2plot_RN_"//trim(cTemp(1:3))//"K.dat"
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
        write(u, 1007) Ji, branch2delta(ai,i),&
                       Jip, branch2delta(af,j), &
                       dta2%WT0(imatrix)
                      
       if (imatrix .gt. nMmx) then
           call sizeError('1001',imatrix,nMmx)
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
  subroutine save_Yrp(dta1, n, molP, Y_RosenP)
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
    !subroutine Variables
    integer*8, parameter                    :: u=11 !file unit
    character*100                           :: path
    character*6                             :: cTemp
    integer*8                               :: i, j, k
    integer*8                               :: today(3)

! INIT. VAR.
    write(cTemp,'(f5.1)'),molP%Temp
    !path = trim(out_file_path)//trim(out_fil2_RMF)
    path ="Y_air_"//trim(cTemp(1:3))//"K.dat"
    call idate(today)   ! today(1)=day, (2)=month, (3)=year
    open (UNIT = u, FILE = trim(path), STATUS = 'REPLACE', ACTION = 'WRITE')
! HEADER
    write (u, *) "Line mixin Rosenkranz Parameter file; Molecule:", molP%chmol
!    write (u, *) "HITRAN Band:", trim(dta1 % hitBand)
    write (u, *) "Based on HITRAN 2012 database"
    write (u, *) "Author: Teresa Mendaza  "
    write (u, 1005) today(3), today(2), today(1)
    write (u, *) "|WaveNumber(cm-1) |Intensity cm-1/(molecules·cm-2) | Y Patameter  | Jlow | Jup"
   !
    ! DATA
    do i = 1, n   ! do1
       write (u, 1006) dta1%sig(i),dta1%Str(i),&
                       Y_RosenP(i),&
                       int(dta1%J(i,1)),int(dta1%J(i,2))
    end do    ! do1
    
    close (u)
    
1005 format ( ' Last update: ', i4.4, '/', i2.2, '/', i2.2 )
1006 format (F12.6,1x,E10.3,1x,F7.4,1x,i3,1x,i3)

end subroutine save_Yrp
!--------------------------------------------------------------------------------------------------------------------