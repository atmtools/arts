!--------------------------------------------------------------------------------------------------------------------
module module_read
!--------------------------------------------------------------------------------------------------------------------
    interface

        subroutine readMolParam(molP, mP_size, my_mol)
            use module_common_var
            use module_molecSp
            implicit none
            integer*8, intent (out) :: mP_size
            integer*8, intent (in)  :: my_mol !== dta1 % M
            type (dta_MOL), intent (inout)       :: molP
        end subroutine readMolParam

        subroutine Hit2DTA(dta1, dta_size, nLines, enough_Lines, artsM, artsI, &
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
            integer*8, intent (in)  :: artsM(nLines), artsI(nLines)
            integer*8, intent (in)  :: artsg0(nLines), artsg00(nLines)
            integer*8, intent (in)  :: artsLow(nLines,4), artsUpp(nLines,4)
            Double Precision      , intent (in)  :: sgmin, sgmax
            Double Precision      , intent (in)  :: artsWNO(nLines),artsS(nLines), &
                                                    artsGA(nLines), artsE00(nLines), &
                                                    artsNA(nLines)
            logical               , intent (out) :: enough_Lines
            integer*8, intent (out) :: dta_size
            type (dta_SDF), intent (inout)       :: dta1
        end subroutine Hit2DTA

    end interface

end module module_read
!--------------------------------------------------------------------------------------------------------------------
subroutine readMolParam(molP, mP_size, my_mol)
!--------------------------------------------------------------------------------------------------------------------
! "readMolParam": READ Molecule Parameters
! ...............................................................
! . Subroutine to read HITRAN "molParam.txt" from its ASCII file.
! ...............................................................
!
! Input/Output Parameters of Routine (Arguments or Common)
! ----------------------------------
! molP   : Array of the number of lines (Output).
! mP_size: number of lines (Output).
!
! Other important Output Quantities (through Common Statements)
! ---------------------------------
!
! Accessed Files:
! --------------
!   'molparam.txt' HITRAN information about molecules/isotop.
!
! Called By: Main Program     
! ---------
!
! Double Precision Version
!     
! T. Mendaza last change 09 Nov 2015
!--------------------------------------------------------------------------------------------------------------------
!
      use module_common_var
      use module_molecSp
      Implicit None
      integer*8, intent (out) :: mP_size
      integer, intent (in)                 :: my_mol !== dta1 % M
      type (dta_MOL), intent (inout)       :: molP
      integer*8               :: i, j, k
      integer*8               :: aux_molno, g_j
      Double Precision                     :: Qmo
      character*58                         :: moLine
      character*100                        :: fname
      integer*8               :: error_read
      logical                              :: error_open      
!
!--------------------------
! Example of HITRAN12 line:
!--------------------------
! Molecule # Iso Abundance     Q(296K)      gj    Molar Mass(g)
!C H2O (1)
!C        161  .997317E+00    1.7464E+02    1     18.010565
!
    mP_size = 1
    error_read = 0
    fname = trim(in_file_path)//trim(in_file_molp)
    call openFile(u5, fname, error_open)
    do
      read (u5, 1001, iostat = error_read, err = 100, end = 101), moLine
      !stop
      if( ( moLine(8:8) .eq. '(' ) .and. (moLine(11:11) .eq. ')') ) then
        read(moLine(9:10), '(i2)'), aux_molno
        if (aux_molno .eq. my_mol) then
        !print*, moLine
        !print*, moLine(8:8), moLine(9:10), moLine(11:11)
        !stop
          do
            read (u5, *, iostat = error_read, err = 101, end = 101), &
            molP % Aco(mP_size), molP % IAb(mP_size), Qmo, &
            g_j, molP % mms(mP_size) 
            !print*, my_mol, "iso#",mP_size, molP % IAb(mP_size)
            ! molP % iso(mP_size) = mP_size !== row number
            mP_size = mP_size + 1
          enddo
        endif     
        
      endif       
      100 if (error_read .ne. 0) cycle ! Reading Error => skip to next available data
      101 if (error_read .ne. 0) exit  ! End of file   => stops "do-loop"

        
    end do
    !stop
    mP_size = mP_size - 1
    call closeFile(u5, error_open)
!
1001  Format(A58)
1002  Format(9x,I3,2x,E7.2,4x,E6.2,3x,I2,4x,F10.6)

    
      Return
end subroutine readMolParam
!--------------------------------------------------------------------------------------------------------------------
subroutine Hit2DTA(dta1, dta_size, nLines, enough_Lines, artsM, artsI, &
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
      integer*8, intent (in)  :: nLines
      integer*8, intent (in)  :: artsM(nLines), artsI(nLines)
      integer*8, intent (in)  :: artsg0(nLines), artsg00(nLines)
      integer*8, intent (in)  :: artsLow(4,nLines), artsUpp(4,nLines)
      Double Precision, intent (in)  :: sgmin, sgmax
      Double Precision, intent (in)  :: artsWNO(nLines),artsS(nLines), &
                                        artsGA(nLines), artsE00(nLines), &
                                        artsNA(nLines)
      integer*8, intent (out) :: dta_size
      logical  , intent (out) :: enough_Lines
      integer*8               :: i, j, k
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
!
1001   Format(I2,I1,f12.6,2E10.3,2F5.4,F10.4,F4.2,F8.6, &
             4A15,I6,I12,A1,2F7.1)
1002   Format(I5,10X)
1003   Format(//,1x,'************ PROBLEM !!!! ******************', &
       1x,'Arrays in for Line data storage are too small',&
       1x,'raise the value of nLmx in the "driver" ')
    
      Return
end subroutine Hit2DTA