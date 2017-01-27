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
            ! MODULES IN USE:
            use module_common_var
            use module_error
            use module_OC_files
                use module_maths
                use module_molecSp
                use module_read
                    use module_phsub
                    use module_linemixing

            Implicit none
            !   INPUT variables
            integer*8         :: nLines, npert
            integer*8         :: artsM(nLines), artsI(nLines)
            integer*8         :: artsg0(nLines), artsg00(nLines)
            integer*8         :: artsLow(4,nLines), artsUpp(4,nLines)
            integer*8         :: pert(npert), i_pert(npert)
            Double Precision  :: sgmin, sgmax, T, Ptot
            Double Precision  :: QT, QT0, mass
            Double Precision  :: p_mass(npert), p_vmr(npert)
            Double Precision  :: artsWNO(nLines),artsS(nLines), &
                                 artsGA(nLines), artsE00(nLines), &
                                 artsNA(nLines)
            !   OUTPUT variables
            Double Precision:: rho(nLines), dipo(nLines) 
            Double Precision:: W_rn(nLines,nLines)
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
! for molecule: CH4.
!
!	INPUT VARIABLES
!	----------------------
!	nLines	 : number of Lines of the band
!   hit      : Array with hitran spectroscopy data.
!   T        : Temperature of the system.
!   molecN   : Hitran ID of the absorption molecule.
!   i_molec  : Molecule's isotope to be considered in this calculation.
!   npert    : Anumber of perturbers/colliders.
!   pert     : Array containing HitranID of the system perturber (colliders)
!   i_pert   : Array of the colliders's isotopes.
!   p_vmr    : pertubers volume mixing ratio.
!   sigmin   : Minimum wavenumber of the interval where the absorption is calculated.
!   sigmax   : Maximum wavenumber of the interval where the absorption is calculated.
!
!	OUTPUT VARIABLES
!	-----------------
!	Dipo	 : Dipole transition Moments of the Lines
!	rho 	 : Populations of the Lower Levels of the Lines
!	           at 296 K.
!   W_rn     : Renormalized Relaxation Matrix acording to the renormalization procedure described in
!              Niro et al. (2004).
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
    use module_OC_files
        use module_maths
        use module_molecSp
        use module_read
           use module_phsub
            use module_linemixing

    Implicit none
!   INPUT variables
    integer*8 :: nLines, npert
    integer*8 :: artsM(nLines), artsI(nLines)
    integer*8 :: artsg0(nLines), artsg00(nLines)
    integer*8 :: artsLow(4,nLines), artsUpp(4,nLines)
    integer*8 :: pert(npert), i_pert(npert)
    Double Precision  :: sgmin, sgmax, T, Ptot
    Double Precision  :: QT, QT0, mass
    Double Precision  :: p_mass(npert), p_vmr(npert)
    Double Precision  :: artsWNO(nLines),artsS(nLines), &
                       artsGA(nLines), artsE00(nLines), &
                       artsNA(nLines)
    integer*8 :: dta_size1, dta_size2
    integer*8 :: iLine, i, j, k
!   OUTPUT variables
    Double Precision       :: rho(nLines), dipo(nLines) 
    Double Precision       :: W_rn(nLines,nLines)
!    integer*8 :: i_pert(npert) 
    Double Precision       :: Wmat(nWmax,nWmax),&
                       Wper(nWmax,nWmax)
!   Double Precision       :: Y_RosT(nLmx)
    Double Precision       :: xMOLp(npert)
    type (dta_SDF)         :: dta1
    type (dta_RMF)         :: dta2
    type (dta_MOL)         :: molP
    type (dta_MOL)         :: PerM
    integer*8 :: sys(2)  
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
    CALL moleculeID(artsM(1), artsI(1), molP)
    dta1%M = molP%M
!---------
! Call for reading HITRAN spectroscopy data.
!
    print*, 'Reading HITRAN12 File...'
    !CALL readline(dta1, dta_size1)
    CALL Hit2DTA(dta1, dta_size1, nLines, enough_Lines, artsM, artsI, &
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
        CALL moleculeID(pert(i),i_pert(i), PerM)
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
        !print*,perturber, PerM%iso_m,PerM%mms(PerM%iso_m)
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
! Write Y parameter file
!
!    print*, 'Saving Rosenkranz parameter Y...'
!    CALL outFile_Yrp(u9, dta1, dta_size1, Y_RosT(1:nLines))
!
!
    PRINT *, "Successful run!"
!   
    STOP
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
