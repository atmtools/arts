!--------------------------------------------------------------------------------------------------------------------
module module_read
!--------------------------------------------------------------------------------------------------------------------
    interface

        subroutine Hit2DTA(dta1, dta_size, nLines, enough_Lines, &
                                            artsWNO, &
                                            artsS, &
                                            artsGA, &
                                            artsE00, &
                                            artsNA , &
                                            artsUpp, artsLow, &
                                            artsg0 , artsg00, &
                                            sgmin  , sgmax)
            use module_common_var
            use module_error
            use module_molecSp
            implicit none
            integer*8, intent (in)  :: nLines
            integer*8, intent (in)  :: artsg0(nLines), artsg00(nLines)
            integer*8, intent (in)  :: artsLow(4,nLines), artsUpp(4,nLines)
            Double Precision      , intent (in)  :: sgmin, sgmax
            Double Precision      , intent (in)  :: artsWNO(nLines),artsS(nLines), &
                                                    artsGA(nLines), artsE00(nLines), &
                                                    artsNA(nLines)
            logical               , intent (out) :: enough_Lines
            integer*8, intent (out) :: dta_size
            type (dta_SDF), intent (inout)       :: dta1
        end subroutine Hit2DTA

        subroutine moleculeID(my_mol,isotope, Mass, PFmol_T, PFmol_T0, flagON, molP)
            use module_common_var
            use module_error
            implicit none
            integer*8, intent (in)             :: my_mol , isotope
            Double Precision, intent (in)      :: Mass, PFmol_T0, PFmol_T
            logical         , intent (in)      :: flagON
            type (dta_MOL), intent(inout)      :: molP
        end subroutine moleculeID

        subroutine molid_char(molP)
            use module_common_var
            use module_error
            implicit none
            type (dta_MOL), intent(inout)      :: molP
        end subroutine molid_char

        subroutine addMolParam(molP)
            use module_common_var
            use module_error
            implicit none
            type (dta_MOL), intent(inout)      :: molP
        end subroutine addMolParam

    end interface

end module module_read
!--------------------------------------------------------------------------------------------------------------------
subroutine Hit2DTA(dta1, dta_size, nLines, enough_Lines, &
                                            artsWNO, &
                                            artsS, &
                                            artsGA, &
                                            artsE00, &
                                            artsNA , &
                                            artsUpp, artsLow, &
                                            artsg0 , artsg00, &
                                            sgmin  , sgmax)
!--------------------------------------------------------------------------------------------------------------------
! "Hit2DTA": READ LINEs data
! .....................................................
! . Subroutine to read HITRAN data from its ASCII file.
! .....................................................
!
! Input/Output Parameters of Routine (Arguments or Common)
! ----------------------------------
! nLines   : Integer Array of the number of lines (Output).
!
! Other important Output Quantities (through Common Statements)
! ---------------------------------
! Sig    : WaveNumbers of the Lines (Cm-1) 
! E    : Energies of the Lower levels of the lines (Cm-1)
! HWT0   : Air-broadened (or H2-broadened) Half-Widths (at 296 K) of the 
!           Lines (Cm-1/Atm)
! BHW    : Temperature Dependence Coefficients of HWT0
!
! Accessed Files:
! --------------
!   'CH4_2nu3_F2iso1_HIT12.dat' Spectroscopic data file of the CH4 2·nu3 band.
!
! Called By: Main Program     
! ---------
!
! Double Precision Version
!     
! T. Mendaza last change 11 May 2015
!--------------------------------------------------------------------------------------------------------------------
!
      use module_common_var
      use module_error
      use module_molecSp
      Implicit None
      integer*8, intent (in)         :: nLines
      integer*8, intent (in)         :: artsg0(nLines), artsg00(nLines)
      integer*8, intent (in)         :: artsLow(4,nLines), artsUpp(4,nLines)
      Double Precision, intent (in)  :: sgmin, sgmax
      Double Precision, intent (in)  :: artsWNO(nLines),artsS(nLines), &
                                        artsGA(nLines), artsE00(nLines), &
                                        artsNA(nLines)
      integer*8, intent (out)        :: dta_size
      logical  , intent (out)        :: enough_Lines
      integer*8                      :: i, j, k
      type (dta_SDF), intent (inout) :: dta1
      integer                        :: error_read
!
!--------------------------
! Example of line: HITRAN12
!--------------------------
! PART 1:
!C/MI/wno......../S......../A......../GamA/GamS/E''....../nAi/shift..
!C 61 4414.089492 2.960E-23 2.520E-02.05500.078  219.91340.70-.008800
! PART 2:
!C/V'............/V''.........../Q'............/Q''.........../Ierr./Iref.......*/g'..../g''...||
!C    0 0 1 1 1F1    0 0 0 0 1A1    7E  43         6E   1     334332453638 7 1 7    30.0   26.0
! 
! but all Lines from "CH4_2nu3_HIT12.dat" has:
!            V'                | V"
!            nu1 2 3 4 C alph' | nu1 2 3 4 C alph"
!         "    0 0 2 0 1F1"    |"  0 0 0 0 1A1"
!
! which is the way to recognise the band 2·nu3
!
    i = 1
    do j = 1,nLines
        if( (artsWNO(i) .ge. sgmin  ) .and. (artsWNO(i) .le. sgmax  ) ) then    
        ! & the frequency lies in the band-interval

          if (i .eq. 1) then
              dta1%lv2(1) = artsLow(1,1) ! Lower level vibrationa angular momentum
              dta1%lv2(2) = artsUpp(1,1) ! Upper level vibrationa angular momentum
          endif 
          dta1 % sig(i)    = artsWNO(i)
          dta1 % Str(i)    = artsS(i)
          dta1 % HWT0(i)   = artsGA(i)
          dta1 % BHW(i)    = artsNA(i)
          dta1 % E(i)      = artsE00(i)
          dta1 % swei0(i)  = artsg0(i)
          dta1 % swei00(i) = artsg00(i)
          call r_arts_LocalQ(dta1, i, artsUpp(:,i), artsLow(:,i))

          i=i+1
          ! Check how many lines has read the subrutine...
          if ( i.gt.nLmx ) then
           call sizeError("1001",i,nLmx)
          endif 
        endif       
        
    end do
    dta_size = i - 1
    if (dta_size .lt. 10) then
      enough_Lines = .false.
    else
      enough_Lines = .true.
    endif  
    Return
end subroutine Hit2DTA
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE moleculeID(my_mol,isotope, Mass, PFmol_T, PFmol_T0, flagON, molP)
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
    use module_error
    implicit none
    integer*8, intent (in)        :: my_mol, isotope
    double precision, intent (in) :: PFmol_T, PFmol_T0, Mass
    logical         , intent (in) :: flagON
    type (dta_MOL), intent(inout) :: molP
    integer*8                     :: i,j,k, mP_size
    integer*8                     :: g_j
    double precision              :: Qmo
    character(  6)                :: auxmol
    character(100)                :: fname
    logical                       :: molfound
    integer*8                     :: error_read
    logical                       :: error_open

!----------
! 
    molP % M = my_mol
    molP % iso_m = isotope
    call molid_char(molP)
    molP % mms = Mass
    if (flagON) then
     molP % QT0 = PFmol_T0
     molP % QT  = PFmol_T
     CALL addMolParam(molP) 
    endif
!
! 
    Return
  END SUBROUTINE moleculeID
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE molid_char(molP)
!--------------------------------------------------------------------------------------------------------------------
! "molid_char": string-name per molecule ID
! 
! Detailed Description:
! ---------------------
! This subroutine check if the molecule-isotope specified has associatted 
! string-name within HITRAN.
!
! Variables:
!
! Input/Output Parameters of Routine (Arguments or Common)
! ----------------------------------
! molP  : Molecule's basic information.
!
! Accessed Files:  'self-database'
! --------------
!
! Called Routines: None
! ---------------  
!
! Called By: moleculeID
! ---------
!
!
! T.Mendaza, last change 24 January 2017
!--------------------------------------------------------------------------------------------------------------------
!
    use module_common_var
    use module_error
    implicit none
    type (dta_MOL), intent(inout) :: molP
    character(  6)                :: auxmol
    integer*8                     :: m, n
    integer*8, DIMENSION(47,11), PARAMETER :: isotopollist = reshape( (/ 161, 181, 171, 162, 182, 172,   0,   0,   0,   0,   0, & !H20
                                                                         626, 636, 628, 627, 638, 637, 828, 827, 727, 838, 837, & !CO2
                                                                         666, 668, 686, 667, 676,   0,   0,   0,   0,   0,   0, & !O3
                                                                         446, 456, 546, 448, 447,   0,   0,   0,   0,   0,   0, & !N2
                                                                          26,  36,  28,  27,  38,  37,   0,   0,   0,   0,   0, & !CO
                                                                         211, 311, 212, 312,   0,   0,   0,   0,   0,   0,   0, & !CH4
                                                                          66,  68,  67,   0,   0,   0,   0,   0,   0,   0,   0, & !O2
                                                                          46,  56,  48,   0,   0,   0,   0,   0,   0,   0,   0, & !NO
                                                                         626, 646,   0,   0,   0,   0,   0,   0,   0,   0,   0, & !SO2
                                                                         646,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, & !NO2
                                                                        4111,5111,   0,   0,   0,   0,   0,   0,   0,   0,   0, & !NH3
                                                                         146, 156,   0,   0,   0,   0,   0,   0,   0,   0,   0, & !HNO3
                                                                          61,  81,  62,   0,   0,   0,   0,   0,   0,   0,   0, & !OH
                                                                          19,  29,   0,   0,   0,   0,   0,   0,   0,   0,   0, & !HF
                                                                          15,  17,  25,  27,   0,   0,   0,   0,   0,   0,   0, & !HCl
                                                                          19,  11,  29,  21,   0,   0,   0,   0,   0,   0,   0, & !HBr
                                                                          17,  27,   0,   0,   0,   0,   0,   0,   0,   0,   0, & !HI
                                                                          56,  76,   0,   0,   0,   0,   0,   0,   0,   0,   0, & !ClO
                                                                         622, 624, 632, 623, 822,   0,   0,   0,   0,   0,   0, & !OCS
                                                                         126, 136, 128,   0,   0,   0,   0,   0,   0,   0,   0, & !H2CO
                                                                         165, 167,   0,   0,   0,   0,   0,   0,   0,   0,   0, & !HOCl
                                                                          44,  45,   0,   0,   0,   0,   0,   0,   0,   0,   0, & !N2
                                                                         124, 134, 125,   0,   0,   0,   0,   0,   0,   0,   0, & !HCN
                                                                         215, 217,   0,   0,   0,   0,   0,   0,   0,   0,   0, & !CH3Cl
                                                                        1661,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, & !H2O2
                                                                        1221,1231,1222,   0,   0,   0,   0,   0,   0,   0,   0, & !C2H2
                                                                        1221,1231,   0,   0,   0,   0,   0,   0,   0,   0,   0, & !C2H6
                                                                        1111,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, & !PH3
                                                                         269, 369,   0,   0,   0,   0,   0,   0,   0,   0,   0, & !COF2
                                                                          29,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, & !SF6
                                                                         121, 141, 131,   0,   0,   0,   0,   0,   0,   0,   0, & !H2S
                                                                         126,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, & !HCOOH
                                                                         166,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, & !HO2
                                                                           6,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, & !O
                                                                        5646,7646,   0,   0,   0,   0,   0,   0,   0,   0,   0, & !ClONO2
                                                                          46,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, & !NO+
                                                                         169, 161,   0,   0,   0,   0,   0,   0,   0,   0,   0, & !HOBr
                                                                         221, 231,   0,   0,   0,   0,   0,   0,   0,   0,   0, & !C2H4
                                                                        2161,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, & !CH3OH
                                                                         219, 211,   0,   0,   0,   0,   0,   0,   0,   0,   0, & !CH3Br
                                                                        2124,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, & !CH3CN
                                                                          29,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, & !CF4
                                                                        2211,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, & !C4H2
                                                                        1224,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, & !HC3N
                                                                          11,  12,   0,   0,   0,   0,   0,   0,   0,   0,   0, & !H2
                                                                          22,  24,  32,  23,   0,   0,   0,   0,   0,   0,   0, & !CS
                                                                          26,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0/),&!SO3
                                                                         shape(isotopollist), order=(/2,1/) )
!-----------------------------------------
! 
    auxmol = ""
    m = molP%M
    !n = molP % iso_m - 10*m
    n = molP % iso_m
    !
    !Molecule's AFGL code
    molP % Aco = isotopollist(m,n) 
    !
    if ( m .eq. 1 ) then
          auxmol = "H2O" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !         161  .997317E+00    1.7464E+02    1     18.010565
          !         181  1.99983E-03    1.7511E+02    1     20.014811
          !         171  3.71884E-04    1.0479E+03    6     19.014780
          !         162  3.10693E-04    8.5901E+02    6     19.016740
          !         182  6.23003E-07    8.7519E+02    6     21.020985
          !         172  1.15853E-07    5.2204E+03   36     20.020956
    elseif ( m .eq. 2 ) then
          auxmol = "CO2" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !         626  .984204E+00    2.8694E+02    1     43.989830
          !         636  1.10574E-02    5.7841E+02    2     44.993185
          !         628  3.94707E-03    6.0948E+02    1     45.994076
          !         627  7.33989E-04    3.5527E+03    6     44.994045
          !         638  4.43446E-05    1.2291E+03    2     46.997431
          !         637  8.24623E-06    7.1629E+03   12     45.997400
          !         828  3.95734E-06    3.2421E+02    1     47.998322
          !         827  1.47180E-06    3.7764E+03    6     46.998291
          !         727  1.36847E-07    1.1002E+04    1     45.998262
          !         838  4.44600E-08    6.5350E+02    2     49.001675
          !         837  1.65354E-08    7.6152E+03   12     48.001646
    elseif ( m .eq. 3 ) then
          auxmol = "O3" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !         666  .992901E+00    3.4838E+03    1     47.984745
          !         668  3.98194E-03    7.4657E+03    1     49.988991
          !         686  1.99097E-03    3.6471E+03    1     49.988991
          !         667  7.40475E-04    4.3331E+04    6     48.988960
          !         676  3.70237E-04    2.1405E+04    6     48.988960
    elseif ( m .eq. 4 ) then
          auxmol = "N2O" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !         446  .990333E+00    5.0018E+03    9     44.001062
          !         456  3.64093E-03    3.3619E+03    6     44.998096
          !         546  3.64093E-03    3.4586E+03    6     44.998096
          !         448  1.98582E-03    5.3147E+03    9     46.005308
          !         447  3.69280E-04    3.0971E+04   54     45.005278
    elseif ( m .eq. 5 ) then
          auxmol = "CO" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !          26  .986544E+00    1.0712E+02    1     27.994915
          !          36  1.10836E-02    2.2408E+02    2     28.998270
          !          28  1.97822E-03    1.1247E+02    1     29.999161
          !          27  3.67867E-04    6.5934E+02    6     28.999130
          !          38  2.22250E-05    2.3582E+02    2     31.002516
          !          37  4.13292E-06    1.3809E+03   12     30.002485
    elseif ( m .eq. 6 ) then
          auxmol = "CH4" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !         211  .988274E+00    5.9052E+02    1     16.031300
          !         311  1.11031E-02    1.1808E+03    2     17.034655
          !         212  6.15751E-04    4.7954E+03    3     17.037475
          !         312  6.91785E-06    9.5990E+03    6     18.040830
    elseif ( m .eq. 7 ) then
          auxmol = "O2" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !          66  .995262E+00    2.1577E+02    1     31.989830
          !          68  3.99141E-03    4.5230E+02    1     33.994076
          !          67  7.42235E-04    2.6406E+03    6     32.994045
    elseif ( m .eq. 8 ) then
          auxmol = "NO" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !          46  .993974E+00    1.1421E+03    3     29.997989
          !          56  3.65431E-03    7.8926E+02    2     30.995023
          !          48  1.99312E-03    1.2045E+03    3     32.002234
    elseif ( m .eq. 9 ) then
          auxmol = "SO2" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !         626  .945678E+00    6.3403E+03    1     63.961901
          !         646  4.19503E-02    6.3689E+03    1     65.957695
    elseif ( m .eq. 10 ) then
          auxmol = "NO2" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !         646  .991616E+00    1.3578E+04    3     45.992904

    elseif ( m .eq. 11 ) then
          auxmol = "NH3" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !        4111  .995872E+00    1.7252E+03    3     17.026549
          !        5111  3.66129E-03    1.1527E+03    2     18.023583
    elseif ( m .eq. 12 ) then
          auxmol = "HNO3" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !         146  .989110E+00    2.1412E+05    6     62.995644
          !         156  3.63600E-03    1.4187E+05    4     63.992680
    elseif ( m .eq. 13 ) then
          auxmol = "OH" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !          61  .997473E+00    8.0362E+01    2     17.002740
          !          81  2.00014E-03    8.0882E+01    2     19.006986
          !          62  1.55371E-04    2.0931E+02    3     18.008915
    elseif ( m .eq. 14 ) then
          auxmol = "HF" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !          19  .999844E+00    4.1469E+01    4     20.006229
          !          29  1.55741E-04    1.1591E+02    6     21.012404
    elseif ( m .eq. 15 ) then
          auxmol = "HCl" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !          15  .757587E+00    1.6065E+02    8     35.976678
          !          17  .242257E+00    1.6089E+02    8     37.973729
          !          25  1.18005E-04    4.6278E+02   12     36.982853
          !          27  3.77350E-05    4.6412E+02   12     38.979904
    elseif ( m .eq. 16 ) then
          auxmol = "HBr" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !          19  .506781E+00    2.0017E+02    8     79.926160
          !          11  .493063E+00    2.0023E+02    8     81.924115
          !          29  7.89384E-05    5.8640E+02   12     80.932336
          !          21  7.68016E-05    5.8676E+02   12     82.930289
    elseif ( m .eq. 17 ) then
          auxmol = "HI" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !          17  .999844E+00    3.8899E+02   12    127.912297
          !          27  1.55741E-04    1.1470E+03   18    128.918472
    elseif ( m .eq. 18 ) then
          auxmol = "ClO" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !          56  .755908E+00    3.2746E+03    4     50.963768
          !          76  .241720E+00    3.3323E+03    4     52.960819
    elseif ( m .eq. 19 ) then
          auxmol = "OCS" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !         622  .937395E+00    1.2210E+03    1     59.966986
          !         624  4.15828E-02    1.2535E+03    1     61.962780
          !         632  1.05315E-02    2.4842E+03    2     60.970341
          !         623  7.39908E-03    4.9501E+03    4     60.966371
          !         822  1.87967E-03    1.3137E+03    1     61.971231
    elseif ( m .eq. 20 ) then
          auxmol = "H2CO" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !         126  .986237E+00    2.8467E+03    1     30.010565
          !         136  1.10802E-02    5.8376E+03    2     31.013920
          !         128  1.97761E-03    2.9864E+03    1     32.014811
    elseif ( m .eq. 21 ) then
          auxmol = "HOCl" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !         165  .755790E+00    1.9274E+04    8     51.971593
          !         167  .241683E+00    1.9616E+04    8     53.968644
    elseif ( m .eq. 22 ) then
          auxmol = "N2" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !          44  .992687E+00    4.6598E+02    1     28.006148
          !          45  7.47809E-03    3.8646E+02    6     29.003182
    elseif ( m .eq. 23 ) then
          auxmol = "HCN" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !         124  .985114E+00    8.9529E+02    6     27.010899
          !         134  1.10676E-02    1.8403E+03   12     28.014254
          !         125  3.62174E-03    6.2141E+02    4     28.007933
    elseif ( m .eq. 24 ) then
          auxmol = "CH3Cl" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !         215  .748937E+00    5.7916E+04    4     49.992328
          !         217  .239491E+00    5.8834E+04    4     51.989379
    elseif ( m .eq. 25 ) then
          auxmol = "H2O2" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !        1661  .994952E+00    9.8198E+03    1     34.005480
    elseif ( m .eq. 26 ) then
          auxmol = "C2H2" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !        1221  .977599E+00    4.1403E+02    1     26.015650
          !        1231  2.19663E-02    1.6562E+03    8     27.019005
          !        1222  3.04550E-04    1.5818E+03    6     27.021825
    elseif ( m .eq. 27 ) then
          auxmol = "C2H6" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !        1221  .976990E+00    7.0881E+04    1     30.046950
          !        1231  2.19526E-02    3.6191E+04    2     31.050305
    elseif ( m .eq. 28 ) then
          auxmol = "PH3" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !        1111  .999533E+00    3.2486E+03    2     33.997238
    elseif ( m .eq. 29 ) then
          auxmol = "COF2" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !         269  .986544E+00    7.0044E+04    1     65.991722
          !         369  1.10834E-02    3.7844E+04    2     66.995083
    elseif ( m .eq. 30 ) then
          auxmol = "SF6" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !          29  .950180E+00    1.6233E+06    1    145.962492
    elseif ( m .eq. 31 ) then
          auxmol = "H2S" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !         121  .949884E+00    5.0307E+02    1     33.987721
          !         141  4.21369E-02    5.0435E+02    1     35.983515
          !         131  7.49766E-03    2.0149E+03    4     34.987105
    elseif ( m .eq. 32 ) then
          auxmol = "HCOOH" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !         126  .983898E+00    3.9133E+04    4     46.005480
    elseif ( m .eq. 33 ) then
          auxmol = "HO2" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !         166  .995107E+00    4.3004E+03    2     32.997655
    elseif ( m .eq. 34 ) then
          auxmol = "O" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !           6  .997628E+00    6.7212E+00    1     15.994915
    elseif ( m .eq. 35 ) then
          auxmol = "ClONO2" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !         5646  .749570E+00    4.7884E+06   12     96.956672
          !         7646  .239694E+00    4.9102E+06   12     98.953723
    elseif ( m .eq. 36 ) then
          auxmol = "NO+" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !          46  .993974E+00    3.1168E+02    3     29.997989
    elseif ( m .eq. 37 ) then
          auxmol = "HOBr" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !         169  .505579E+00    2.8339E+04    8     95.921076
          !         161  .491894E+00    2.8238E+04    8     97.919027
    elseif ( m .eq. 38 ) then
          auxmol = "C2H4" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !         221  .977294E+00    1.1041E+04    1     28.031300
          !         231  2.19595E-02    4.5197E+04    2     29.034655
    elseif ( m .eq. 39 ) then
          auxmol = "CH3OH" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !        2161  .985930E+00    3.5314E+04    2     32.026215
    elseif ( m .eq. 40 ) then
          auxmol = "CH3Br" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !         219  .500995E+00    8.3051E+04    4     93.941811
          !         211  .487433E+00    8.3395E+04    4     95.939764
    elseif ( m .eq. 41 ) then
          auxmol = "CH3CN" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !        2124  .973866E+00    8.8671E+04    3     41.026549
    elseif ( m .eq. 42 ) then
          auxmol = "CF4" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !          29  .988890E+00    1.2127E+05    1     87.993616
    elseif ( m .eq. 43 ) then
          auxmol = "C4H2" 
         ! Molecule #    Iso Abundance  Q(296K)       gj    Molar Mass(g)
         !         2211  .955998E+00    9.8180E+03    1     50.015650
    elseif ( m .eq. 44 ) then
          auxmol = "HC3N" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !        1224  .963346E+00    2.4785E+04    6     51.010899
    elseif ( m .eq. 45 ) then
          auxmol = "H2" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !          11  .999688E+00    7.6724E+00    1     2.0156500
          !          12  3.11432E-04    2.9879E+01    6     3.0218250
    elseif ( m .eq. 46 ) then
          auxmol = "CS" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !          22  .939624E+00    2.5362E+02    1     43.971036
          !          24  .416817E-01    2.5777E+02    1     45.966787
          !          32  .105565E-01    5.3750E+02    2     44.974368
          !          23  .741668E-02    1.0230E+03    4     44.970399
    elseif ( m .eq. 47 ) then
          auxmol = "SO3" 
          ! Molecule #   Iso Abundance  Q(296K)       gj    Molar Mass(g)
          !          26  .943400E+00    7.9638E+03    1     79.956820
    else
          call molnameERROR(molP%M)
    endif
    molP % chmol = auxmol

    if (molP % Aco .eq. 0) call isoAconameERROR(molP%M,molP%iso_m)
!
! 
!   Data from HITRAN molparam.txt file
!   ----------------------------------
!   Dry air mixing ratio from
!   GLOBALVIEW-CH4: Cooperative Atmospheric Data Integration
!   Project - Methane. CD-ROM, NOAA ESRL, Boulder, Colorado
!   [Also available on Internet via anonymous FTP to ftp.cmdl.noaa.gov,
!   Path: ccg/ch4/GLOBALVIEW], 2009.
!
    Return
  END SUBROUTINE molid_char
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
    use module_error
    implicit none
    type (dta_MOL), intent(inout) :: molP
    
!----------
! 
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
        if (molP%Aco .eq. 66) then
            ! AFGL code = 66 => 16O2 (isotope 1 in HITRAN)
            ! If nu1 = 0 : B0 = 43100.430 MHz 
            !                 = 1.437 cm-1 (MHz*10^6/c)  
            ! If nu1 = 1 : B0 = 42626.398 MHz 
            !                 = 1.421 cm-1 
            molP % B0    = 1.43  ! O2 Rotational constant B0 (cm-1) 
                                 ! NIST 
        elseif (molP%Aco .eq. 67) then
            ! AFGL code = 67 => 16O17O (isotope 3 in HITRAN)
            ! If nu1 = 0 : B0 = 40561.35 MHz 
            !                 = 1.353 cm-1 
            molP % B0    = 1.35  ! O2 Rotational constant B0 (cm-1) 
                                 ! NIST 
        elseif (molP%Aco .eq. 68) then
            ! AFGL code = 68 => 16O18O (isotope 2 in HITRAN)
            ! If nu1 = 0 : B0 = 38313.761 MHz 
            !                 = 1.278 cm-1 
            ! If nu1 = 1 : B0 = 37916.618 MHz 
            !                 = 1.265 cm-1 
            molP % B0    = 1.27  ! O2 Rotational constant B0 (cm-1) 
                                 ! NIST 
        else 
            call addMolError("1003",  molP%iso_m)
        endif
        
    ELSE 
        call addMolError("1003",  molP%M)
    ENDIF
! 
!   Data from NIST website
!   -----------------------
!
      Return
  END SUBROUTINE addMolParam