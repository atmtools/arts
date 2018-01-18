MODULE module_linemixing
!--------------------------------------------------------------------------------------------------------------------
! This module contains all subroutines related to molecular symmetry and behavior
!
    interface

        subroutine WelCAL(dta1, nLines, molP, PerM, W_jk, econ)
            use module_common_var
            use module_error
            use module_maths
            use module_phsub
            implicit none
            integer*8              , intent(in   )  :: nLines
            type (dta_SDF), pointer, intent(in   )  :: dta1
            type (dta_MOL)         , intent(in   )  :: molP
            type (dta_MOL)         , intent(in   )  :: PerM
            type (dta_ERR)         , intent(inout)  :: econ
            double precision       , intent(  out)  :: W_jk(nLines,nLines)
        end subroutine WelCAL

        subroutine RN_Wmat(nLines, dta1, Wmat, W_rnO, T, P, econ)
            use module_common_var
            use module_error
            use module_maths
            implicit none
            integer*8              , intent(in   ) :: nLines
            Double precision       , intent(in   ) :: T, P
            type (dta_SDF), pointer, intent(in   ) :: dta1
            Double Precision       , intent(in   ) :: Wmat(nLines,nLines)
            Double Precision       , intent(  out) :: W_rnO(nLines,nLines)
            type (dta_ERR)         , intent(inout) :: econ
        end subroutine RN_Wmat

        subroutine LM_Rosen(molP, nLines,dta1,Wmat,Y_RosT)
            use module_common_var
            use module_maths
            use module_phsub
            implicit none
            integer*8              , intent(in   ) :: nLines 
            double precision       , intent(in   ) :: Wmat(nLines,nLines)
            type (dta_MOL)         , intent(in   ) :: molP
            type (dta_SDF), pointer, intent(inout) :: dta1
            double precision       , intent(  out) :: Y_RosT(nLines)
        end subroutine LM_Rosen

        subroutine LM_2ord(molP, nLines,dta1,Wmat,Y2,Y3)
            use module_common_var
            use module_maths
            use module_phsub
            implicit none
            integer*8              , intent(in   ) :: nLines 
            double precision       , intent(in   ) :: Wmat(nLines,nLines)
            type (dta_MOL)         , intent(in   ) :: molP
            type (dta_SDF), pointer, intent(inout) :: dta1
            double precision       , intent(  out) :: Y2(nLines),Y3(nLines)
        end subroutine LM_2ord

        subroutine calc_QParam(nLines, dta1, molP, PerM, econ)
            use module_common_var
            use module_error
            use module_maths
            use module_LLS
            Implicit None
            integer*8              , intent(in   ) :: nLines
            type (dta_SDF), pointer, intent(in   ) :: dta1
            type (dta_MOL)         , intent(inout) :: molP
            type (dta_MOL)         , intent(in   ) :: PerM
            type (dta_ERR)         , intent(inout) :: econ
        end subroutine calc_QParam

        subroutine calc_QParam_AF(nLines, dta1, molP, PerM, econ)
            use module_common_var
            use module_error
            use module_maths
            use module_LLS
            Implicit None
            integer*8              , intent(in   ) :: nLines
            type (dta_SDF), pointer, intent(in   ) :: dta1
            type (dta_MOL)         , intent(inout) :: molP
            type (dta_MOL)         , intent(in   ) :: PerM
            type (dta_ERR)         , intent(inout) :: econ
        end subroutine calc_QParam_AF

        subroutine calc_QPar_DGELSY(nLines, dta1, molP, PerM, econ)
            use module_common_var
            use module_error
            use module_maths
            use module_LLS
            Implicit None
            integer*8              , intent(in   ) :: nLines
            type (dta_SDF), pointer, intent(in   ) :: dta1
            type (dta_MOL)         , intent(inout) :: molP
            type (dta_MOL)         , intent(in   ) :: PerM
            type (dta_ERR)         , intent(inout) :: econ
        end subroutine calc_QPar_DGELSY

        subroutine calc_QPar_DGELSS(nLines, dta1, molP, PerM, econ)
            use module_common_var
            use module_error
            use module_maths
            use module_LLS
            Implicit None
            integer*8              , intent(in   ) :: nLines
            type (dta_SDF), pointer, intent(in   ) :: dta1
            type (dta_MOL)         , intent(inout) :: molP
            type (dta_MOL)         , intent(in   ) :: PerM
            type (dta_ERR)         , intent(inout) :: econ
        end subroutine calc_QPar_DGELSS

        subroutine calc_QPar_DGELSD(nLines, dta1, molP, PerM, econ)
            use module_common_var
            use module_error
            use module_maths
            use module_LLS
            Implicit None
            integer*8              , intent(in   ) :: nLines
            type (dta_SDF), pointer, intent(in   ) :: dta1
            type (dta_MOL)         , intent(inout) :: molP
            type (dta_MOL)         , intent(in   ) :: PerM
            type (dta_ERR)         , intent(inout) :: econ
        end subroutine calc_QPar_DGELSD

        subroutine W2dta2(nLines, dta1, dta2, W_rnO)
            use module_common_var
            use module_maths
            use module_phsub
            implicit none
            integer*8              , intent(in   ) :: nLines
            double precision       , intent(in   ) :: W_rnO(nLines,nLines)
            type (dta_SDF), pointer, intent(in   ) :: dta1
            type (dta_RMF), pointer, intent(inout) :: dta2
        end subroutine W2dta2

        logical function rule1(nLines)
            use module_common_var
            implicit none
            integer*8, intent(in) :: nLines
        end function rule1

        logical function rule2(P,nLines,W,v,TOL_rule2)
            use module_common_var
            implicit none
            integer*8       , intent(in) :: nLines
            double precision, intent(in) :: P, TOL_rule2
            double precision, intent(in) :: W(nLines,nLines)
            double precision, intent(in) :: v(nLines) 
        end function rule2

    end interface

END module module_linemixing
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE WelCAL(dta1, nLines, molP, PerM, W_jk, econ)
!--------------------------------------------------------------------------------------------------------------------
! WelCAL: W ELements CALculation (Relaxation matrix element)
!
! Detailed description:
! ---------------------
! The matrix stored is the real part of of the relaxation operator, based on the 
! ECS approximation and properly accounts for the coupling and relaxation of rotational angular momenta
! (i.e. where each element represents the rotational state-to-state cross sections 
! within a single vibrational state). 
! 
! Variables:
!
! Input/Output Parameters of Routine (Arguments or Common)
! --------------------------------------------------------
! nLines  : Number of lines of A band (Input).
! HWT0    : Air (or H2) Broadened HalfWidths of the Lines for
!           T0 and P0 (Cm-1). (Input).
! PopuT0  : Populations of the Lower Levels of the Lines
!           at Temperature T0. (Input).
! W:jk    : Relaxation matrix (Output) 
!
! Accessed Files:  None
! --------------
!     
! Called Routines: "W_error"              
! ---------------  "K_jkO2"/"K_jkCalc"
!     
! Called By: arts_interface : "RM_LM_LLS_tmc_arts" or "RM_LM_tmc_arts"  
! ---------
!     
!     Double Precision Version
!     
! T. Mendaza last change 08 Jan 2016
!--------------------------------------------------------------------------------------------------------------------
!
    use module_common_var
    use module_error
    use module_phsub
    use module_maths
    Implicit None
    integer*8             , intent(in)    :: nLines
    type (dta_SDF),pointer, intent(in)    :: dta1
    type (dta_MOL)        , intent(in)    :: molP
    type (dta_MOL)        , intent(in)    :: PerM
    type (dta_ERR)        , intent(inout) :: econ
    double Precision      , intent(out)   :: W_jk(nLines,nLines)
    integer*8                       :: i, j, k, n
    integer*8                       :: jBIG, jSMALL
    integer*8                       :: indexI(nLines)
    Double Precision                :: r_kj, pn, pk
    Double Precision                :: auxW, auxHW
    Double Precision                :: Wtest(nLines,nLines)
    !Auxiliar Constants
    integer*8                       :: count1, count2
    double Precision                :: Kaux, HWT, faH
    double Precision                :: T, Ptot, RT
    logical                         :: tOK
!
!-----------------------------------------
    T    = molP % Temp
    Ptot = molP % Ptot
    RT   = T0/T
!-----------------------------------------
!
    faH = 1.0_dp
!
    do j=1, nLines
      do k=1, j
        !write(*,'(a3,i5,a3,i5,a1)') '(j=',j,';k=',k,')'         
        if ( j .eq. k ) then ! Qlow(j) == Qlow(k)
        ! Diagonal levels of the Matrix (Wjj)
        ! They are defined as 
        ! the Lorenz Half-Width + a line shift correction
        !           (dependent on the Molecule)
        !>> Wjj= HWT0(j) + i*shift(j)
        ! The Halfwidth dependent on Temperature:
        ! HWT0 = HITRAN(296K) air halfwidth (molecule)
        ! BHW = air temperature exponent (molecule)
        !
          if(T.ne.T0) faH = (RT**dta1%BHW(j)) 
          HWT = dta1%HWT0(j)*faH!* (1-xH2O) &
        ! NOTE: 
        ! in case one would like to include water vapour broadening
        ! add the "!!!" terms (and provide data): 
        !!! + (xH2O*(HWT0_H2O(i)) 
        ! NOTE: HWT0_H2O depends on Temperature
        ! in case one has the temperature dependency exponent add the 
        ! following line instead:
        !!! + (xH2O*(HWT0_H2O(i)*(RT**BHW_H2O(i))))
        !
        ! where:
        ! HWT0_H2O = halfwidth broadeing by water vapour (at 296K)
        ! xH2O  = H2O atmospheric percentage 
        ! BHW_H2O = water vapour temperature exponent
        !
          W_jk(j,k) = 2*Ptot*HWT !+ i*(-0.008)  
        !
        else 
        ! OFF-diagonal matrix elements (Wjk ∞ W_jk)
        ! Wjk:= initial state <<k|W|j>>
        ! where:
        ! <<k|W|j>> := transition j->k
        ! and "j" > "k" (downwards transition)
        !
        ! It is assumed that the relaxation matrix elements have
        ! the same functional form as the rotational state-to-state 
        ! cross sections within a single vibrational state W_j,k := W_jk
        !            
          if (isJb(dta1,j,k)) then
          ! CASE:  J(j) > J(k) (downwards transition j->k)
          ! or
          ! CASE: J(j) = J(k)
            jBIG   = j
            jSMALL = k
          else
          ! CASE: J(j) < J(k)
          ! so downwards transition is (k->j)
          ! pj·<<k|W|j>> = pk·<<j|W|k>>; pk = dta1%PopuT(k); pj = dta1%PopuT(j)
            jBIG   = k
            jSMALL = j
          endif
          if (molP%M .eq. 7) then
              ! O2
                W_jk(jBIG,jSMALL) = K_jkO2(jBIG,jSMALL,dta1,nLines,molP,PerM,econ) 
          else 
                W_jk(jBIG,jSMALL) = K_jkCalc(j,jBIG,jSMALL,dta1,nLines,molP,PerM,econ)  
          endif
          !                                
          ! Downwards transition is (k->j)
          ! == downwards transition:(jBIG -> jSMALL)
          ! pjSMALL·<<jBIG|W|jSMALL>> = pjBIG·<<jSMALL|W|jBIG>>; 
          ! where:
          ! pjBIG = dta1%PopuT(jBIG); pjSMALL = dta1%PopuT(jSMALL)
          !
          r_kj = dta1%PopuT(jBIG)/dta1%PopuT(jSMALL) !pjBIG/pjSMALL
          W_jk(jSMALL,jBIG) = r_kj*W_jk(jBIG,jSMALL)  
              !
          if ( isnan( W_jk(jBIG,jSMALL) ) .or. isinf( W_jk(jBIG,jSMALL) ) ) then
              call W_error("1010", jBIG, jSMALL, W_jk(jBIG,jSMALL), econ)
              call W_error("1011", jSMALL, jBIG, W_jk(jSMALL,jBIG), econ)
          endif
        ! 
        endif
      !    
      enddo
    !
    enddo
    !
    Return
  END SUBROUTINE WelCAL
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE RN_Wmat(nLines, dta1, Wmat, W_rnO, T, P,econ)
!--------------------------------------------------------------------------------------------------------------------
! RN_Wmat: Renormalization of the Relaxation matrix element.
!
! Detailed description:
! ---------------------
! Subroutine to renormalize the relaxation matrix to the number of transitions.
! Niro et al. (2006) writes that renormalization is crucial to built a relaxation matrix that simultaneously:
! (1) respects the detailed balance. 
! (2) exactly verifies the sum rule. <---- Here is where renormalization becomes essential!
! (3) has precise diagonal elements in order to correctly model the isolated lines.
!
! The first (1) and (2) are essential to obtain the correct behavior in the wings. 
! (1) and (2) can be obtained simultaneously through a renormalization of our previous ECS model.
! 
! Variables:
!
! Input/Output Parameters of Routine (Arguments or Common)
! --------------------------------------------------------
! nLines  : Number of lines of A band (Input).
! Wmat    : initial Relaxation Matrix (Input).
! PopuT   : Populations of the Lower Levels of the Lines
!           at Temperature T (Input).
! W_rnO   : renormalized Relaxation Matrix (Output).
!
! Accessed Files:  None
! --------------
!     
! Called Routines: "isnan"/"isinf"
! ---------------  "renorm_error"
!                  "bubble_index"/"ibubble_index"
!                  "sumRule"
!     
! Called By: arts_interface : "RM_LM_LLS_tmc_arts" or "RM_LM_tmc_arts" 
! ---------
!     
!     Double Precision Version
!     
! T. Mendaza last change 17 Jan 2018
!--------------------------------------------------------------------------------------------------------------------
!
    use module_common_var
    use module_error
    use module_maths
    Implicit None
    integer*8             , intent(in   )  :: nLines
    Double precision      , intent(in   )  :: T, P
    Double precision      , intent(in   )  :: Wmat(nLines,nLines)
    Double Precision      , intent(  out)  :: W_rnO(nLines,nLines)
    type (dta_SDF),pointer, intent(in   )  :: dta1
    type (dta_ERR)        , intent(inout)  :: econ
    !Local variables
    integer*8                     :: i, j, k, n, e20
    integer*8                     :: indexS(nLines),indexI(nLines)
    double Precision              :: sortV(nLines)
    Double Precision              :: Sup, Slow, S_UL_rate
    Double Precision              :: pn, pk
    Double Precision              :: W_rn(nLines,nLines)
!---------
! Reordering process
!
      do k=1,nLines
        ! The rule to sort out the matrix can be any of the followings:
        ! 1) The spectral line intensity [cm-1/(moleculecm-1)] at "296 K
        !    (similar to HITRAN's intensity, however at that variable-HITRAN-  
        !    the natural terrestrial isotopic abundance is included)
        !
         !sortV(k) = dta1%STR(k)
        !
        ! 2) Reduced line intensity (or "strenght") that uses the rigid rotor
        !
         !sortV(k) = dta1%sig(k)*dta1%popuT(k)*dta1%D0(k)**2 
        !
        ! 3) Reduced line intensity (or "strenght") that uses the dipole (depends on temperature) 
        !
         !sortV(k) = dta1%sig(k)*dta1%popuT(k)*dta1%DipoT(k)**2
        !
        ! 4) using the population of the lines (at the given temperature) as the leading factor to
        !    measure the importance of the line.
        !
         sortV(k) = dta1%popuT(k) 
      enddo
      ! The rotational components are sorted according to their "importance" 
      ! (intensity according to [Niro et al. 2004] but this code will use the one selected above) 
      ! in decreasing order (the first one is the most "intense"). Hence, << 1 | W | 1 >> is the 
      ! broadening of the most intense line, and << 2 | W | 1 >> (or << 1 | W | 2 >>) are the terms 
      ! coupling this transition to the second most intense one.
      !
      ! SORT by intensity/Populations:
      call bubble_index(nLines,sortV,indexS,'d',econ)
      ! Here we perform the 'pullback' of indexS. 
      call ibubble_index(nLines,indexS,indexI,'a',econ)
      !NOTE: the sorting algorithm used is the classic Bubble sort
      !      Although the algorithm is simple, it is too slow and 
      !      impractical for most problems even when compared to INSERTION SORT. 
      !      This part should be improved in the future.
      !
      do i=1,nLines
          do j=1,nLines
            if (i.eq.j) then
              W_rn(i,i)=Wmat(indexS(i), indexS(i))
            else
              W_rn(i,j)=-abs(Wmat(indexS(i), indexS(j)))
            endif
          enddo
      enddo
      !
      ! Then, each column "k" of the matrix, starting from the first one is treated as follows:
      DO n = 1, nLines
        Sup = 0.0; Slow = 0.0
        pn = dta1%popuT( indexS(n) )
        DO k = 1, nLines
          ! remove assymetry
          if(dta1%iso.gt.2 .AND. dta1%iso.ne.7 &
            .AND.&           
            mod(abs(int(dta1%J(n,1))-int(dta1%J(k,1))),2).ne.0)cycle
          if (k .le. n) then 
            Sup  = Sup  + dabs( dta1%D0( indexS(k) ) )* W_rn( n , k )
            if(isnan(Sup).or.isinf(Sup)) then
              call renorm_error("1008", n, k, W_rn( n , k ), Sup, econ)
            endif
          else !if (k .gt. n) then
            Slow = Slow + dabs(dta1%D0( indexS(k) )) * W_rn( n , k ) 
            if(isnan(Slow).or.isinf(Slow))then
              call renorm_error("1009", n, k, W_rn( n , k ), Slow, econ)
            endif          
          endif
        ENDDO
        S_UL_rate=Sup/Slow
        !
        ! Scale the elements of the "lower part" of the column (k>n)
        ! So we have the W_RN (renormalize Relaxation Matrix)
        DO k = n,nLines
          if ( k .ne. n) then 
          !
            if (Slow .eq. 0.0_dp) then
              W_rn(n,k) = 0.0_dp
              W_rn(k,n) = 0.0_dp
            else
              ! A) The lower-elements of the matrix, in other words, if ( n < k ) then 
              ! "n" has a lower index-value == most intense line than "k", that mark a
              ! downwards transition and this is expressed as: W_jk= - W · Sup/Slow
              ! W(n,k) = << k | W | n >> 
              !   -->       k<------n
              !
              W_rn(n,k) = -S_UL_rate*W_rn(n,k)
              !
              ! B) The upper-elements of the matrix, in other words, the corresponding 
              ! transposed line of the matrix is built from detailed balance:
              ! pk·<< n | W | k >> = pn·<< k | W | n >>  
              ! pk·   W_rn(k, n)   = pn·   W_rn(n,k)
              !
              pk = dta1%popuT(indexS(k))
              W_rn(k,n) = W_rn(n,k)*pn/pk
            endif
          endif
        ENDDO  
        !
      ENDDO
      !
      ! Sum Rule TEST
      !
      ! The Relaxation matrix must follow the SUM-RULE
      ! However, that is acomplish regarding that 
      ! every relaxation between every theoretical transition
      ! (not just empirically observed) is taken into the sum.
      ! Since the calc. is performed using lines from a database the  
      ! will not be significant before a renormalization.
      !
      ! NOTE:
      ! use '1.0' if Wii = line-width 
      ! use '2.0' if Wii = half-width 
      e20 = econ % e(2)
      CALL sumRule(nLines,indexS,dta1%D0(1:nLines),W_rn,1.0,econ)
      !
      if ( econ % e(2) .gt. e20 ) then
        ! if the test fails, return a diagonal matrix.
        !
        econ % e(2) = econ % e(2) - 1
        CALL just_fill_DiagWRn(nLines,dta1 % BHW, dta1 % HWT0, T/T0, P,W_rnO)
      else 
        ! if the matrix passes the sumRule test, the matrix will be returned
        ! in its original order.
        !
        do i=1,nLines
          do j=1,nLines
            W_rnO(i,j)  =  W_rn( indexI(i) , indexI(j) )
          enddo
        enddo
        !
      endif
  END SUBROUTINE RN_Wmat
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE LM_Rosen(molP, nLines,dta1,Wmat,Y_RosT)
!--------------------------------------------------------------------------------------------------------------------
! "LM_Rosen": Rosenkranz parameter
! 
! Detailed Description:
! ---------------------
! This subroutine gives the First order Linemixing coefficient (a.k.a Rosenkranz Parameter).
!
! Variables:
!
! Input/Output Parameters of Routine 
! ----------------------------------
! nLines  : Number of lines of A band (Input).
! Wmat    : renormalized Relaxation Matrix  (Input).
! PopuT   : Populations of the Lower Levels of the Lines
!           at Temperature T (Input).
! DipoT   : Dipole moment of the the Line
!           at Temperature T (Input).
! sig     : line wavenumber  (Input).
! Str     : line intensity as in equation (5) [Simeckova et al. 2006]  (Input).
! Y_RosT  : 1st order LM-parameter of the line (Output).
!
! Called By: "RM_LM_LLS_tmc_arts" or "RM_LM_tmc_arts"
! ---------
!
!
! T.Mendaza, last change 29 March 2017
!--------------------------------------------------------------------------------------------------------------------
!
    use module_common_var
    use module_maths
    use module_phsub
    implicit none
    integer*8             , intent(in   ) :: nLines 
    Double Precision      , intent(in   ) :: Wmat(nLines,nLines)
    double precision      , intent(  out) :: Y_RosT(nLines)
    type (dta_SDF),pointer, intent(inout) :: dta1
    type (dta_MOL)        , intent(in   ) :: molP
    !Local Var
    integer*8                     :: i,j,k
    double precision              :: delta, sumY, T, Ptot
    double precision              :: DipoI, DipoK
!-----------------------------------------
    T    = molP % Temp
    Ptot = molP % Ptot
!-----------------------------------------
!
!   Build Ym from W elements
!     
    DO i=1,nLines
         sumY=0.d0
         !DipoI= abs(dsqrt(dta1%Str(i)/(dta1%Sig(i)*dta1%PopuT0(i))))
         DipoI  = abs(dsqrt(dta1%Str(i)/(dta1%Sig(i)*dta1%PopuT(i))))
         !DipoI= dta1%DipoT(i)
         do k=1,nLines
            !DipoK= abs(dsqrt(dta1%Str(k)/(dta1%Sig(k)*dta1%PopuT0(k))))
            DipoK = abs(dsqrt(dta1%Str(k)/(dta1%Sig(k)*dta1%PopuT(k))))
            !DipoK= dta1%DipoT(k)
            if(k.eq.i)cycle
            !
            !  Correction for asym hysothopes
            if( (dta1%iso.gt.2) .AND. (dta1%iso .ne. 7) .AND. &
                mod( abs( int(dta1%J(i,1)-dta1%J(k,1)) ),2).ne.0)cycle
            !
            !     Using detailed balance
            delta = dta1%sig(i) - dta1%sig(k)
            !
            if( dabs(delta) .lt. 1.d-4 ) delta=1.d-4
            !
            sumY=sumY + (DipoK/DipoI)* &
                ( Wmat(i,k)/delta ) ! Wmat(k,i)
        enddo
        Y_RosT(i)=sumY
    ENDDO
    Y_RosT=2.0*Ptot*Y_RosT
  END SUBROUTINE LM_Rosen
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE LM_2ord(molP, nLines,dta1,Wmat,Y2,Y3)
!--------------------------------------------------------------------------------------------------------------------
! "LM_2ord": Second and third coeff of the second order linemixing formulation
! 
! Detailed Description:
! ---------------------
! This subroutine gives the Second and Third coeff of the second order linemixing formulation.
!
! Variables:
!
! Input/Output Parameters of Routine 
! ----------------------------------
! nLines  : Number of lines of A band (Input).
! Wmat    : renormalized Relaxation Matrix  (Input).
! PopuT   : Populations of the Lower Levels of the Lines
!           at Temperature T (Input).
! DipoT   : Dipole moment of the the Line
!           at Temperature T (Input).
! sig     : line wavenumber  (Input).
! Str     : line intensity as in equation (5) [Simeckova et al. 2006]  (Input).
! Y2      : g_k, second coeffcient from the second order linemixing approx [Smith, 1981]. (Output).
! Y3      :dv_k, is the second order shift [Smith, 1981]. (Output).
!
! NOTE: on the calculation of these coeff (as well the 1st order one), 
! the first order pressure shift has been neglected, as it was done by Smith (1981).
! 
!
! Called By: "RM_LM_LLS_tmc_arts" or "RM_LM_tmc_arts"
! ---------
!
!
! T.Mendaza, last change 29 March 2017
!--------------------------------------------------------------------------------------------------------------------
!
    use module_common_var
    use module_maths
    use module_phsub
    implicit none
    integer*8             , intent(in   ) :: nLines 
    Double Precision      , intent(in   ) :: Wmat(nLines,nLines)
    double precision      , intent(  out) :: Y2(nLines),Y3(nLines)
    type (dta_SDF),pointer, intent(inout) :: dta1
    type (dta_MOL)        , intent(in   ) :: molP
    !Local Var
    integer*8                     :: i,l,k
    double precision              :: delta, deltaL, T, Ptot, Pto2
    double precision              :: sumG1, sumG2, sumG3, sumG4, sumG42
    double precision              :: sumDV
    double precision              :: DipoI, DipoK, rD_ki
!-----------------------------------------
    T    = molP % Temp
    Ptot = molP % Ptot
    Pto2 = Ptot**2
!-----------------------------------------
!
!     Build the Y2 and Y3 from the W
!     
    DO i=1,nLines
         sumG1  = 0.d0
         sumG2  = 0.d0
         sumG3  = 0.d0
         sumG4  = 0.d0
         sumG42 = 0.d0
         sumDV  = 0.d0
         !DipoI  = abs(dsqrt(dta1%Str(i)/(dta1%Sig(i)*dta1%PopuT0(i))))
         DipoI  = abs(dsqrt(dta1%Str(i)/(dta1%Sig(i)*dta1%PopuT(i))))
         !DipoI= dta1%DipoT(i)
         do k=1,nLines
            !DipoK = abs(dsqrt(dta1%Str(k)/(dta1%Sig(k)*dta1%PopuT0(k))))
            DipoK = abs(dsqrt(dta1%Str(k)/(dta1%Sig(k)*dta1%PopuT(k))))
            !DipoK= dta1%DipoT(k)
            !
            rD_ki = DipoK/DipoI
            if (isnan(rD_ki).or.isinf(rD_ki))rD_ki = 1.0d0
            if(k.eq.i)cycle
            !
            !  Correction for asym hysothopes
            if( (dta1%iso.gt.2) .AND. (dta1%iso .ne. 7) .AND. &
                mod( abs( int(dta1%J(i,1)-dta1%J(k,1)) ),2).ne.0)cycle
            !
            !     Using detailed balance
            delta = dta1%sig(k) - dta1%sig(i)
            if( dabs(delta) .lt. 1.d-4 ) delta=1.d-4
            !
            sumG1 = sumG1 + Wmat(i,k)*Wmat(k,i)/(delta**2)
            sumG2 = sumG2 + rD_ki*( Wmat(i,k)/delta ) 
            sumG3 = sumG3 + rD_ki*( Wmat(i,k)*Wmat(i,i)/(delta**2) )
            do l = 1,nLines
              if (l.eq.i) then
                cycle
              else
                deltaL= dta1%sig(l) - dta1%sig(i)
                if( dabs(deltaL) .lt. 1.d-4 ) deltaL=1.d-4
                sumG42 = sumG42 + ( Wmat(l,k)*Wmat(i,l)/(delta*deltaL) )
              endif
            enddo
            sumG4 = sumG4 + rD_ki*sumG42
            sumDV = sumDV + Wmat(i,k)*Wmat(k,i)/delta
        enddo
        Y2(i)=sumG1 - (sumG2)**2 + 2.0d0*sumG3 - 2.0d0*sumG4
        Y3(i)=sumDV
    ENDDO
    Y2=Y2*Pto2
    Y3=Y3*Pto2
  END SUBROUTINE LM_2ord  
!
!--------------------------------------------------------------------------------------------------------------------
! LINEAR LEAST SQUARES Problems
! -----------------------------
! In the most usual case m>=n and {rank}(A) = n, and in this case the solution to min||Ax-B|| is 
! unique, and the problem is also referred to as finding a least squares solution to an overdetermined system 
! of linear equations.
!
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE calc_QParam(nLines, dta1, molP, PerM, econ)
!--------------------------------------------------------------------------------------------------------------------
! "calc_QParam": gives the values of a1, a2, a3 
! after solving a Linear-Least-Square system. For that Lapack is used.
!
!------------------------------------------ 
    use module_common_var
    use module_error
    use module_maths
    use module_LLS
    Implicit None
    integer*8              , intent(in   ) :: nLines
    type (dta_SDF), pointer, intent(in   ) :: dta1
    type (dta_MOL)         , intent(inout) :: molP
    type (dta_MOL)         , intent(in   ) :: PerM
    type (dta_ERR)         , intent(inout) :: econ
    double Precision :: Mlls(nLines,4)
    double Precision :: Aux_4M(4)
    integer*8        :: i, j, k
    integer*8        :: jBIG, jSMALL
    integer*8        :: indexI(nLines)
    Double Precision :: r_kj, rD0_kj
    Double Precision :: K1,K2, sumK2,sumK1
    !Auxiliar Constants
    double Precision :: Kaux, HWT, faH
    double precision :: T, Ptot, RT
    !Auxiliar for LAPACK
    !
    ! Parameters:
      integer*8, Parameter   :: N = 3
      integer*8, Parameter   :: NRHS = 1
      integer*8              :: M
      integer*8              :: LDA, LDB
      !PARAMETER        ( LDA = M, LDB = M )
      integer*8, Parameter   :: LWMAX = 100
    ! Local Scalars:
      integer*8              :: LWORK
      integer*8              :: info1, info2
    ! Local Arrays:
      Double Precision       :: WORK( LWMAX )
      Double Precision,ALLOCATABLE,DIMENSION(:,:) :: A( :, : ), B( :, : )
!-----------------------------------------
      T    = molP % Temp
      Ptot = molP % Ptot
      RT   = T0/T
!-----------------------------------------
      M = nLines
      LDA = M
      LDB = M
      allocate ( A( LDA, N ), B( LDB, NRHS ) )
!
! LAPACK is used (installation command for mac):
! sudo port install lapack
! 
! * Compilation for "Free PGI compiler" for MAC:
! pgf90 myprog.f90 -llapack -lblas
!
! * Compilation "gfortran"
! gfortran myprog.f90 -llapack
!
!---------
! FIRST: create a zero matrix for 
    do j=1, nLines
      do k=1,4
        Mlls(j,k) = 0.0_dp
      enddo  
    enddo  
!
!---------
! Generate the Matrix for LLS:
    do j=1, nLines
      do k=1, nLines
        ! 
        if (j .eq. k) then
          faH = 1.0_dp
          if(T.ne.T0)faH = (RT**dta1%BHW(j)) 
          B(j,NRHS) = 2*molP%Ptot*dta1%HWT0(j)*faH
        else              
          if (isJb(dta1,j,k)) then
          ! CASE:  J(j) > J(k) (downwards transition j->k)
          ! or
          ! CASE: J(j) = J(k)
            jBIG   = j
            jSMALL = k
            r_kj = 1.0_dp 
          else
          ! CASE: J(j) < J(k)
          ! so downwards transition is (k->j)
          ! pj·<<k|W|j>> = pk·<<j|W|k>>; pk = dta1%PopuT(k); pj = dta1%PopuT(j)
            jBIG   = k
            jSMALL = j
            r_kj = dta1%PopuT(jBIG)/dta1%PopuT(jSMALL) !pjBIG/pjSMALL
          endif
          call LLS_Matrix(jBIG,jSMALL,dta1,molP,PerM,Aux_4M, econ)
          !
          rD0_kj = dta1%D0(k)/dta1%D0(j) 
          do i =1,4
            Mlls(j,i) = Mlls(j,i) + rD0_kj*r_kj*Aux_4M(i)
          enddo  
        endif    
      enddo
      !
    enddo
!
! **********************************************************************************
! LAPACK routine:
! ---------------
! DGELS solves overdetermined or undetermined real linear systems
! involving an M-by-N matrix A, or its transpose, using RQ or LQ factorization of A.
! Note_ it is asumed that A has full rank.
! if TRANS = 'N' and m>=n: find the least squares solution of an overdetermined system, 
! i.e., solve the least squares problem minimize || B - A*X ||
!  A * X = B,
!
!---------
! Definition of A, B:
    !
    A = -Mlls(1:nLines,1:3)
    !
    do i = 1,nLines
      B(i,NRHS) = B(i,NRHS) + Mlls(i,4)
    enddo    
!
!
!  -- LAPACK driver routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
!      CHARACTER          TRANS
!      INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS
!     ..
!     .. Array Arguments ..
!      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DGELS solves overdetermined or underdetermined real linear systems
!  involving an M-by-N matrix A, or its transpose, using a QR or LQ
!  factorization of A.  It is  that A has full rank.
!
!  The following options are provided:
!
!  1. If TRANS = 'N' and m >= n:  find the least squares solution of
!     an overdetermined system, i.e., solve the least squares problem
!                  minimize || B - A*X ||.
!
!  2. If TRANS = 'N' and m < n:  find the minimum norm solution of
!  an underdetermined system A * X = B.
!
!  3. If TRANS = 'T' and m >= n:  find the minimum norm solution of
!  an undetermined system A**T * X = B.
!
!  4. If TRANS = 'T' and m < n:  find the least squares solution of
!  an overdetermined system, i.e., solve the least squares problem
!               minimize || B - A**T * X ||.
!
!  Several right hand side vectors b and solution vectors x can be
!  handled in a single call; they are stored as the columns of the
!  M-by-NRHS right hand side matrix B and the N-by-NRHS solution
!  matrix X.
!
!  Arguments
!  =========
!
!  TRANS   (input) CHARACTER*1
!       = 'N': the linear system involves A;
!       = 'T': the linear system involves A**T.
!
!  M       (input) INTEGER
!       The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!       The number of columns of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!       The number of right hand sides, i.e., the number of
!       columns of the matrices B and X. NRHS >=0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!       On entry, the M-by-N matrix A.
!       On exit,
!         if M >= N, A is overwritten by details of its QR
!                    factorization as returned by DGEQRF;
!         if M <  N, A is overwritten by details of its LQ
!                    factorization as returned by DGELQF.
!
!  LDA     (input) INTEGER
!       The leading dimension of the array A.  LDA >= max(1,M).
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!       On entry, the matrix B of right hand side vectors, stored
!       columnwise; B is M-by-NRHS if TRANS = 'N', or N-by-NRHS
!       if TRANS = 'T'.
!       On exit, if INFO = 0, B is overwritten by the solution
!       vectors, stored columnwise:
!       if TRANS = 'N' and m >= n, rows 1 to n of B contain the least
!       squares solution vectors; the residual sum of squares for the
!       solution in each column is given by the sum of squares of
!       elements N+1 to M in that column;
!       if TRANS = 'N' and m < n, rows 1 to N of B contain the
!       minimum norm solution vectors;
!       if TRANS = 'T' and m >= n, rows 1 to M of B contain the
!       minimum norm solution vectors;
!       if TRANS = 'T' and m < n, rows 1 to M of B contain the
!       least squares solution vectors; the residual sum of squares
!       for the solution in each column is given by the sum of
!       squares of elements M+1 to N in that column.
!
!  LDB     (input) INTEGER
!       The leading dimension of the array B. LDB >= MAX(1,M,N).
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!       On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!       The dimension of the array WORK.
!       LWORK >= max( 1, MN + max( MN, NRHS ) ).
!       For optimal performance,
!       LWORK >= max( 1, MN + max( MN, NRHS )*NB ).
!       where MN = min(M,N) and NB is the optimum block size.
!
!       If LWORK = -1, then a workspace query is assumed; the routine
!       only calculates the optimal size of the WORK array, returns
!       this value as the first entry of the WORK array, and no error
!       message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!       = 0:  successful exit
!       < 0:  if INFO = -i, the i-th argument had an illegal value
!       > 0:  if INFO =  i, the i-th diagonal element of the
!             triangular factor of A is zero, so that A does not have
!             full rank; the least squares solution could not be
!             computed.
!
!  =====================================================================
!
!  Command sentence:
!  CALL dgels( TRANS, M , N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )
!
!  =====================================================================
!
!   Query the optimal workspace.
!
      LWORK = -1
      info1 = 0
      CALL DGELS( 'No transpose', M, N, NRHS, A, LDA, B, LDB, WORK,&
                 LWORK, info1 )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
!
!---------
!   Solve the equations A*X = B.
!
      info2 = 0
      CALL DGELS( 'No transpose', M, N, NRHS, A, LDA, B, LDB, WORK,&
                 LWORK, info2 )
! <---------------------------------------------
! NOTE: for further information about this subroutine 
! http://www.netlib.no/netlib/lapack/double/dgels.f
! http://www.netlib.no/netlib/lapack/double/dgelsd.f
! http://www.netlib.no/netlib/lapack/double/dgglse.f
!
! For examples:
! https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/dgels_ex.f.htm
!
! For information about LAPACK solving linear methods in general visit:
! http://www.netlib.org/lapack/lug/node26.html
!**********************************************************************************
          molP%a1 = -B(1,NRHS)
          molP%a2 = -B(2,NRHS)
          molP%a3 = -B(3,NRHS)
          if (info2 .ne. 0) then
            call LLS_error(B(1,NRHS),B(2,NRHS),B(3,NRHS), info2, econ)
          endif

  END SUBROUTINE calc_QParam
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE calc_QParam_AF(nLines, dta1, molP, PerM, econ)
!--------------------------------------------------------------------------------------------------------------------
! "calc_QParam_AF": gives the values of a1, a2, a3 
! after solving a Linear-Least-Square system. For that Lapack is used.
!
!------------------------------------------ 
    use module_common_var
    use module_error
    use module_maths
    use module_LLS
    Implicit None
    integer*8              , intent(in   ) :: nLines
    type (dta_SDF), pointer, intent(in   ) :: dta1
    type (dta_MOL)         , intent(inout) :: molP
    type (dta_MOL)         , intent(in   ) :: PerM
    type (dta_ERR)         , intent(inout) :: econ
    double Precision :: Mlls(nLines,12)
    double Precision :: Aux_M(12)
    integer*8        :: i, j, k
    integer*8        :: jBIG, jSMALL
    integer*8        :: indexI(nLines)
    Double Precision :: r_kj, rD0_kj
    !Auxiliar Constants
    double Precision :: Kaux, HWT, faH
    double precision :: T, Ptot, RT
    !Auxiliar for LAPACK
    !
    ! Parameters:
      integer*8, Parameter   :: N = 11
      integer*8, Parameter   :: NRHS = 1
      integer*8              :: M
      integer*8              :: LDA, LDB
      !PARAMETER        ( LDA = M, LDB = M )
      integer*8, Parameter   :: LWMAX = 100
    ! Local Scalars:
      integer*8              :: LWORK
      integer*8              :: info1, info2
    ! Local Arrays:
      Double Precision       :: WORK( LWMAX )
      Double Precision,ALLOCATABLE,DIMENSION(:,:) :: A( :, : ), B( :, : )
!-----------------------------------------
      T    = molP % Temp
      Ptot = molP % Ptot
      RT   = T0/T
!-----------------------------------------
      M = nLines
      LDA = M
      LDB = M
      allocate ( A( LDA, N ), B( LDB, NRHS ) )
!
! LAPACK is used (installation command for mac):
! sudo port install lapack
! 
! * Compilation for "Free PGI compiler" for MAC:
! pgf90 myprog.f90 -llapack -lblas
!
! * Compilation "gfortran"
! gfortran myprog.f90 -llapack
!
!---------
! FIRST: create a zero matrix for 
    do j=1, nLines
      do k=1,(N+1)
        Mlls(j,k) = 0.0_dp
      enddo  
    enddo  
!
!---------
! Generate the Matrix for LLS:
    do j=1, nLines
      do k=1, nLines
        ! 
        if (j .eq. k) then
          faH = 1.0_dp
          if(T.ne.T0)faH = (RT**dta1%BHW(j)) 
          B(j,NRHS) = 2*molP%Ptot*dta1%HWT0(j)*faH
        else              
          if (isJb(dta1,j,k)) then
          ! CASE:  J(j) > J(k) (downwards transition j->k)
          ! or
          ! CASE: J(j) = J(k)
            jBIG   = j
            jSMALL = k
            r_kj = 1.0_dp 
          else
          ! CASE: J(j) < J(k)
          ! so downwards transition is (k->j)
          ! pj·<<k|W|j>> = pk·<<j|W|k>>; pk = dta1%PopuT(k); pj = dta1%PopuT(j)
            jBIG   = k
            jSMALL = j
            r_kj = dta1%PopuT(jBIG)/dta1%PopuT(jSMALL) !pjBIG/pjSMALL
          endif
          call LLS_AF_Matrix(jBIG,jSMALL,dta1,molP,PerM,Aux_M, econ)
          !
          rD0_kj = dta1%D0(k)/dta1%D0(j) 
          do i =1,(N+1)
            Mlls(j,i) = Mlls(j,i) + rD0_kj*r_kj*Aux_M(i)
          enddo  
        endif    
      enddo
      !
    enddo
!
! **********************************************************************************
! LAPACK routine:
! ---------------
! DGELS solves overdetermined or undetermined real linear systems
! involving an M-by-N matrix A, or its transpose, using RQ or LQ factorization of A.
! Note_ it is asumed that A has full rank.
! if TRANS = 'N' and m>=n: find the least squares solution of an overdetermined system, 
! i.e., solve the least squares problem minimize || B - A*X ||
!  A * X = B,
!
!---------
! Definition of A, B:
    !
    A = Mlls(1:nLines,1:N)
    !
    do i = 1,nLines
      B(i,NRHS) = B(i,NRHS) - Mlls(i,N+1)
    enddo    
!
!  -- LAPACK driver routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
!      CHARACTER          TRANS
!      INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS
!     ..
!     .. Array Arguments ..
!      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DGELS solves overdetermined or underdetermined real linear systems
!  involving an M-by-N matrix A, or its transpose, using a QR or LQ
!  factorization of A.  It is  that A has full rank.
!
!  The following options are provided:
!
!  1. If TRANS = 'N' and m >= n:  find the least squares solution of
!     an overdetermined system, i.e., solve the least squares problem
!                  minimize || B - A*X ||.
!
!  2. If TRANS = 'N' and m < n:  find the minimum norm solution of
!  an underdetermined system A * X = B.
!
!  3. If TRANS = 'T' and m >= n:  find the minimum norm solution of
!  an undetermined system A**T * X = B.
!
!  4. If TRANS = 'T' and m < n:  find the least squares solution of
!  an overdetermined system, i.e., solve the least squares problem
!               minimize || B - A**T * X ||.
!
!  Several right hand side vectors b and solution vectors x can be
!  handled in a single call; they are stored as the columns of the
!  M-by-NRHS right hand side matrix B and the N-by-NRHS solution
!  matrix X.
!
!  Arguments
!  =========
!
!  TRANS   (input) CHARACTER*1
!       = 'N': the linear system involves A;
!       = 'T': the linear system involves A**T.
!
!  M       (input) INTEGER
!       The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!       The number of columns of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!       The number of right hand sides, i.e., the number of
!       columns of the matrices B and X. NRHS >=0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!       On entry, the M-by-N matrix A.
!       On exit,
!         if M >= N, A is overwritten by details of its QR
!                    factorization as returned by DGEQRF;
!         if M <  N, A is overwritten by details of its LQ
!                    factorization as returned by DGELQF.
!
!  LDA     (input) INTEGER
!       The leading dimension of the array A.  LDA >= max(1,M).
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!       On entry, the matrix B of right hand side vectors, stored
!       columnwise; B is M-by-NRHS if TRANS = 'N', or N-by-NRHS
!       if TRANS = 'T'.
!       On exit, if INFO = 0, B is overwritten by the solution
!       vectors, stored columnwise:
!       if TRANS = 'N' and m >= n, rows 1 to n of B contain the least
!       squares solution vectors; the residual sum of squares for the
!       solution in each column is given by the sum of squares of
!       elements N+1 to M in that column;
!       if TRANS = 'N' and m < n, rows 1 to N of B contain the
!       minimum norm solution vectors;
!       if TRANS = 'T' and m >= n, rows 1 to M of B contain the
!       minimum norm solution vectors;
!       if TRANS = 'T' and m < n, rows 1 to M of B contain the
!       least squares solution vectors; the residual sum of squares
!       for the solution in each column is given by the sum of
!       squares of elements M+1 to N in that column.
!
!  LDB     (input) INTEGER
!       The leading dimension of the array B. LDB >= MAX(1,M,N).
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!       On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!       The dimension of the array WORK.
!       LWORK >= max( 1, MN + max( MN, NRHS ) ).
!       For optimal performance,
!       LWORK >= max( 1, MN + max( MN, NRHS )*NB ).
!       where MN = min(M,N) and NB is the optimum block size.
!
!       If LWORK = -1, then a workspace query is assumed; the routine
!       only calculates the optimal size of the WORK array, returns
!       this value as the first entry of the WORK array, and no error
!       message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!       = 0:  successful exit
!       < 0:  if INFO = -i, the i-th argument had an illegal value
!       > 0:  if INFO =  i, the i-th diagonal element of the
!             triangular factor of A is zero, so that A does not have
!             full rank; the least squares solution could not be
!             computed.
!
!  =====================================================================
!
!  Command sentence:
!  CALL dgels( TRANS, M , N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )
!
!  =====================================================================
!
!   Query the optimal workspace.
!
      LWORK = -1
      info1 = 0
      CALL DGELS( 'No transpose', M, N, NRHS, A, LDA, B, LDB, WORK,&
                 LWORK, info1 )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
!
!   Solve the equations A*X = B.
!
      info2 = 0
      CALL DGELS( 'No transpose', M, N, NRHS, A, LDA, B, LDB, WORK,&
                 LWORK, info2 )
! <---------------------------------------------
! NOTE: for further information about this subroutine 
! http://www.netlib.no/netlib/lapack/double/dgels.f
! http://www.netlib.no/netlib/lapack/double/dgelsd.f
! http://www.netlib.no/netlib/lapack/double/dgglse.f
!
! For examples:
! https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/dgels_ex.f.htm
!
! For information about LAPACK solving linear methods in general visit:
! http://www.netlib.org/lapack/lug/node26.html
!**********************************************************************************
          molP%a1 = B(1,NRHS)
          molP%a2 = B(2,NRHS)
          molP%a3 = B(3,NRHS)
          molP%a4 = B(4,NRHS)
          molP%a5 = B(5,NRHS)
          molP%a6 = B(6,NRHS)
          molP%a7 = B(7,NRHS)
          molP%a8 = B(8,NRHS)
          molP%a9 = B(9,NRHS)
          if (info2 .ne. 0) then
            call LLS_error(B(1,NRHS),B(2,NRHS),B(3,NRHS), info2, econ)
            call LLS_error(B(4,NRHS),B(5,NRHS),B(6,NRHS), info2, econ)
            call LLS_error(B(7,NRHS),B(8,NRHS),B(9,NRHS), info2, econ)
          endif

  END SUBROUTINE calc_QParam_AF
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE W2dta2(nLines, dta1, dta2, W_rn)
!--------------------------------------------------------------------------------------------------------------------
! W2dta2: Copy the Relaxation Matrix to the struct dta2.
!
! Accessed Files:  None
! --------------
!     
! Called Routines: "isnan"/"isinf"
! ---------------   
!     
! Called By: arts_interface : "RM_LM_LLS_tmc_arts" or "RM_LM_tmc_arts"
! ---------
!     
!     Double Precision Version
!     
! T. Mendaza last change 04 Abr 2016
!--------------------------------------------------------------------------------------------------------------------
!
    use module_common_var
    use module_maths
    Implicit None
    integer*8             , intent(in   ) :: nLines
    double precision      , intent(in   ) :: W_rn(nLines,nLines)
    type (dta_SDF),pointer, intent(in   ) :: dta1
    type (dta_RMF),pointer, intent(inout) :: dta2
    !Local variables
    integer                       :: i, j, k, n
    double precision              :: str(nLines)
    !Strings
    character*15                  :: auxkQupp, auxkQlow
    character*15                  :: auxnQupp, auxnQlow
    character*60                  :: tra2tra
!---------
! Reordering by line strength
!
       DO n = 1,nLines
        !auxnQupp = dta1%Qupp(n)
        !auxnQlow = dta1%Qlow(n)
        DO k = 1,nLines
          !auxkQupp = dta1%Qupp(k)
          !auxkQlow = dta1%Qlow(k)
          IF ( k .eq. n ) THEN
            !dta2%WT0( (nLines*( n - 1 ) + k) ) = W_rn( indexI(n) , indexI(k) )
            dta2%WT0( (nLines*( n - 1 ) + k) ) = W_rn( n , k )
          ELSE
            ! 
            ! Matrix Element
            !
            if (isnan( W_rn( n , k ) ) .or. isinf( W_rn( n , k ) ) .or. &
                W_rn( n , k) .eq. 0.0_dp  ) then
              dta2%WT0( (nLines*( n - 1 ) + k) ) = 0.0_dp
            else
              !dta2%WT0( (nLines*( n - 1 ) + k) ) = W_rn( indexI(n) , indexI(k) )
              dta2%WT0( (nLines*( n - 1 ) + k) ) = -abs(W_rn( n , k ))
            endif  
          ENDIF
         ENDDO
       ENDDO
       !
  END SUBROUTINE W2dta2
!--------------------------------------------------------------------------------------------------------------------
  logical function rule1(nLines)
!--------------------------------------------------------------------------------------------------------------------
! rule1: if the number of lines of a band is not enough, 
!        then calcualting its Relaxation Matrix makes no sense.
!     
! T. Mendaza last change 20 Feb 2017
! --------------------------------------------------------
    implicit none
    integer*8, intent(in)        :: nLines
    !
    if (nLines .lt. 15) then
      rule1 = .false.
    else
      rule1 = .true.
    endif
    RETURN
  end function rule1
!--------------------------------------------------------------------------------------------------------------------
  logical function rule2(P,nLines,W,v,TOL_rule2)
!--------------------------------------------------------------------------------------------------------------------
! rule2: Perturbation Theory limitation.
!
! Detailed description:
! ---------------------
! To invert the MAtrix expresion, the perturbation theory is appropriate under the condition
!
! |  P * W_lk |
! | --------- | << 1
! |  vl - vk  |
! 
! Variables:
! P    = Pressure (atm) 
! W_lk = << k | W | l >> (cm-1/atm)
! vl   = central frequency(cm-1) of line l
! vk   = central frequency(cm-1) of line k
!
! NOTE: The Tolerance used in this function (TOL_rule2)
!       is an arbitrary selection made by the user
!     
! T. Mendaza last update 16 Jan 2018
! --------------------------------------------------------
    use module_common_var
    implicit none
    integer*8       , intent(in) :: nLines
    double precision, intent(in) :: P, TOL_rule2
    double precision, intent(in) :: W(nLines,nLines)
    double precision, intent(in) :: v(nLines) 
    integer                      :: l,k,pos(2)
    double precision             :: aux,minWlk
    !
    rule2=.true.
    !
    DO l =1,nLines
      DO k =1,nLines
      !
        if ( l .ne. k ) then
          aux = P*abs(W(l,k))/abs(v(l)-v(k))
          if (aux .gt. TOL_rule2) then
            rule2 = .false.
          endif
        endif
      ENDDO
    ENDDO

    RETURN
  end function rule2
!--------------------------------------------------------------------------------------------------------------------
! OTHER METHODS: RANK DEFICIENT PROBLEMS
! -------------
! In the general case when we may have {rank}(A) < min(m,n) -- in other words, A may be rank-deficient -- we seek 
! the minimum norm least squares solution x which minimizes both |x|2 and || b - Ax ||_2
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE calc_QPar_DGELSY(nLines, dta1, molP, PerM, econ)
!--------------------------------------------------------------------------------------------------------------------
! "calc_QPar_DGELSY": gives the values of a1, a2, a3 
! assumed A could be rank deficient so it seeks for the minimum norm least squares solution x
! which minimizes both |x|2 and || b - Ax ||_2
! using a complete orthogonal factorization of A
!
!------------------------------------------ 
    use module_common_var
    use module_error
    use module_maths
    use module_LLS
    Implicit None
    integer*8              , intent(in   ) :: nLines
    type (dta_SDF), pointer, intent(in   ) :: dta1
    type (dta_MOL)         , intent(inout) :: molP
    type (dta_MOL)         , intent(in   ) :: PerM
    type (dta_ERR)         , intent(inout) :: econ
    double Precision :: Mlls(nLines,4)
    double Precision :: Aux_4M(4)
    integer*8        :: i, j, k
    integer*8        :: jBIG, jSMALL
    integer*8        :: indexI(nLines)
    Double Precision :: r_kj, rD0_kj
    Double Precision :: K1,K2, sumK2,sumK1
    !Auxiliar Constants
    double Precision :: Kaux, HWT, faH
    double precision :: T, Ptot, RT
    !Auxiliar for LAPACK
    !
    ! Parameters:
      ! Local Scalars:
      integer*8, Parameter   :: N = 3
      integer*8, Parameter   :: NB = 64
      integer*8, Parameter   :: NRHS = 1
      integer*8, Parameter   :: LWMAX = 250
      integer*8              :: M
      integer*8              :: LDA, LDB
      integer*8              :: LWORK, RANK
      integer*8              :: info1, info2
      Double Precision       :: RCOND
      ! Local Arrays:
      integer*8,ALLOCATABLE,DIMENSION(:):: JPVT(:)
      Double Precision,ALLOCATABLE,DIMENSION(:):: WORK( : )
      Double Precision,ALLOCATABLE,DIMENSION(:,:) :: A( :, : ), B( :, : )
!-----------------------------------------
      T    = molP % Temp
      Ptot = molP % Ptot
      RT   = T0/T
!-----------------------------------------
      M = nLines
      LDA = M
      LDB = M
      LWORK = 3*N + NB*( M+N )
      allocate ( JPVT( N ) )
      !allocate ( WORK( LWORK ) )
      allocate ( WORK( LWMAX ) )
      allocate ( A( LDA, N ), B( LDB, NRHS ))
      info1 = 0
!
! LAPACK is used (installation command for mac):
! sudo port install lapack
! 
! * Compilation for "Free PGI compiler" for MAC:
! pgf90 myprog.f90 -llapack -lblas
!
! * Compilation "gfortran"
! gfortran myprog.f90 -llapack
!
!---------
! FIRST: create a zero matrix for 
    do j=1, nLines
      do k=1,4
        Mlls(j,k) = 0.0_dp
      enddo  
    enddo  
!
!---------
! Generate the Matrix for LLS:
    do j=1, nLines
      do k=1, nLines
        ! 
        if (j .eq. k) then
          faH = 1.0_dp
          if(T.ne.T0)faH = (RT**dta1%BHW(j)) 
          B(j,NRHS) = 2*molP%Ptot*dta1%HWT0(j)*faH
        else              
          if (isJb(dta1,j,k)) then
          ! CASE:  J(j) > J(k) (downwards transition j->k)
          ! or
          ! CASE: J(j) = J(k)
            jBIG   = j
            jSMALL = k
            r_kj = 1.0_dp 
          else
          ! CASE: J(j) < J(k)
          ! so downwards transition is (k->j)
          ! pj·<<k|W|j>> = pk·<<j|W|k>>; pk = dta1%PopuT(k); pj = dta1%PopuT(j)
            jBIG   = k
            jSMALL = j
            r_kj = dta1%PopuT(jBIG)/dta1%PopuT(jSMALL) !pjBIG/pjSMALL
          endif
          call LLS_Matrix(jBIG,jSMALL,dta1,molP,PerM,Aux_4M, econ)
          !
          rD0_kj = dta1%D0(k)/dta1%D0(j) 
          do i =1,4
            Mlls(j,i) = Mlls(j,i) + rD0_kj*r_kj*Aux_4M(i)
          enddo  
        endif    
      enddo
      !
    enddo
!
! **********************************************************************************
! LAPACK routine:
! ---------------
! DGELSY computes the minimum-norm solution to a real linear least
! squares problem:
!    minimize || A * X - B ||
! DGELSY solves real linear systems involving an M-by-N matrix A using a complete 
! orthogonal factorization of A. 
! allowing for the possibility that A is rank-deficient.
!===================================================================================
!
! Definition of A, B:
!
! A:
    A = -Mlls(1:nLines,1:3)
! B:
    do i = 1,nLines
      B(i,NRHS) = B(i,NRHS) + Mlls(i,4)
    enddo    
!
! Purpose
! =======
!
! DGELSY computes the minimum-norm solution to a real linear least
! squares problem:
!    minimize || A * X - B ||
! using a complete orthogonal factorization of A.  A is an M-by-N
! matrix which may be rank-deficient.
!
! Several right hand side vectors b and solution vectors x can be
! handled in a single call; they are stored as the columns of the
! M-by-NRHS right hand side matrix B and the N-by-NRHS solution
! matrix X.
!
! The routine first computes a QR factorization with column pivoting:
!    A * P = Q * [ R11 R12 ]
!                [  0  R22 ]
! with R11 defined as the largest leading submatrix whose estimated
! condition number is less than 1/RCOND.  The order of R11, RANK,
! is the effective rank of A.
!
! Then, R22 is considered to be negligible, and R12 is annihilated
! by orthogonal transformations from the right, arriving at the
! complete orthogonal factorization:
!   A * P = Q * [ T11 0 ] * Z
!               [  0  0 ]
! The minimum-norm solution is then
!   X = P * Z' [ inv(T11)*Q1'*B ]
!              [        0       ]
! where Q1 consists of the first RANK columns of Q.
!
! This routine is basically identical to the original xGELSX except
! three differences:
!  o The call to the subroutine xGEQPF has been substituted by the
!    the call to the subroutine xGEQP3. This subroutine is a Blas-3
!    version of the QR factorization with column pivoting.
!  o Matrix B (the right hand side) is updated with Blas-3.
!  o The permutation of matrix B (the right hand side) is faster and
!    more simple.
!
! Arguments
! =========
!
! M       (input) INTEGER
!         The number of rows of the matrix A.  M >= 0.
!
! N       (input) INTEGER
!         The number of columns of the matrix A.  N >= 0.
!
! NRHS    (input) INTEGER
!         The number of right hand sides, i.e., the number of
!         columns of matrices B and X. NRHS >= 0.
!
! A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!         On entry, the M-by-N matrix A.
!         On exit, A has been overwritten by details of its
!         complete orthogonal factorization.
!
! LDA     (input) INTEGER
!         The leading dimension of the array A.  LDA >= max(1,M).
!
! B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!         On entry, the M-by-NRHS right hand side matrix B.
!         On exit, the N-by-NRHS solution matrix X.
!
! LDB     (input) INTEGER
!         The leading dimension of the array B. LDB >= max(1,M,N).
!
! JPVT    (input/output) INTEGER array, dimension (N)
!         On entry, if JPVT(i) .ne. 0, the i-th column of A is permuted
!         to the front of AP, otherwise column i is a free column.
!         On exit, if JPVT(i) = k, then the i-th column of AP
!         was the k-th column of A.
!
! RCOND   (input) DOUBLE PRECISION
!         RCOND is used to determine the effective rank of A, which
!         is defined as the order of the largest leading triangular
!         submatrix R11 in the QR factorization with pivoting of A,
!         whose estimated condition number < 1/RCOND.
!
! RANK    (output) INTEGER
!         The effective rank of A, i.e., the order of the submatrix
!         R11.  This is the same as the order of the submatrix T11
!         in the complete orthogonal factorization of A.
!
! WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!         On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
! LWORK   (input) INTEGER
!         The dimension of the array WORK.
!         The unblocked strategy requires that:
!            LWORK >= MAX( MN+3*N+1, 2*MN+NRHS ),
!         where MN = min( M, N ).
!         The block algorithm requires that:
!            LWORK >= MAX( MN+2*N+NB*(N+1), 2*MN+NB*NRHS ),
!         where NB is an upper bound on the blocksize returned
!         by ILAENV for the routines DGEQP3, DTZRZF, STZRQF, DORMQR,
!         and DORMRZ.
!
!         If LWORK = -1, then a workspace query is assumed; the routine
!         only calculates the optimal size of the WORK array, returns
!         this value as the first entry of the WORK array, and no error
!         message related to LWORK is issued by XERBLA.
!
! INFO    (output) INTEGER
!         = 0: successful exit
!         < 0: If INFO = -i, the i-th argument had an illegal value.
!
! Further Details
! ===============
!  Based on contributions by
!  A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA
!  E. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain
!  G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain
!
!  =====================================================================
!
!  Command sentence:
!  CALL DGELSY( M , N, NRHS, A, LDA, B, M, JPVT, RCOND, RANK, WORK, LWORK, INFO)
!
!  =====================================================================
!  Step0: Query the optimal workspace.
    if (econ % e(1) .ge. 1)  write(*,*) 'prior-Lwork', LWORK
    LWORK = -1
    info1=0
    CALL DGELSY(M,N,NRHS,A,LDA,B,M,JPVT,RCOND,RANK,WORK,LWORK,info1)
    LWORK = INT( WORK( 1 ) )
    !LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
    if (LWORK .lt. LWMAX) then
      deallocate ( WORK )
      allocate ( WORK( LWORK ) )
    endif
    if (econ % e(1) .ge. 1)  write(*,*) 'post-Lwork', LWORK
!  Step1: Initialize JPVT to be zero so that all columns are free
!         <== Done while allocating the variable
!
!  Step2: Choose RCOND to reflect the relative accuracy of the input data
    !RCOND = 0.01D0 !<== These would be hard-imposed accuracy
    RCOND = -1.0D0  !<== < 0, machine precision is used instead.
!
!  Step3: Solve the least squares problem min( norm2(b - Ax) ) for the x
!         of minimum norm.
    info2 = 0
    CALL DGELSY(M,N,NRHS,A,LDA,B,M,JPVT,RCOND,RANK,WORK,LWORK,info2)
!         
! <---------------------------------------------
! NOTE: for further information about this subroutine 
! http://www.netlib.no/netlib/lapack/double/dgelsy.f
!
! For examples:
! https://www.nag.com/lapack-ex/examples/source/dgelsy-ex.f
!
! For information about LAPACK solving linear methods in general visit:
! http://www.netlib.org/lapack/lug/node26.html
!**********************************************************************************
    molP%a1 = -B(1,NRHS)
    molP%a2 = -B(2,NRHS)
    molP%a3 = -B(3,NRHS)
    if (info2 .eq. 0) then
        if (econ % e(1) .ge. 1) write(*,*) 'Complete Orthogonal Factorization algorithm succeded!'
        !
        ! Print the effective rank of A
        !
        if (econ % e(1) .ge. 1) write(*,*) 'Tolerance used to estimate the rank of A', RCOND
        if (econ % e(1) .ge. 1) write(*,*) 'Estimated rank of A', RANK
        !
    else
        call LLS_DGELSYSD_error(B(1,NRHS),B(2,NRHS),B(3,NRHS), info1, econ)
        if (econ % e(1) .ge. 1) write(*,*) 'Complete Orthogonal Factorization of A failed.'
    endif

  END SUBROUTINE calc_QPar_DGELSY
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE calc_QPar_DGELSS(nLines, dta1, molP, PerM, econ)
!--------------------------------------------------------------------------------------------------------------------
! "calc_QPar_DGELSS": gives the values of a1, a2, a3 
! after solving linear system using the Lapack's DGELSS routine.
!
!---------------------------------------------------------------
    use module_common_var
    use module_error
    use module_maths
    use module_LLS
    Implicit None
    integer*8              , intent(in   ) :: nLines
    type (dta_SDF), pointer, intent(in   ) :: dta1
    type (dta_MOL)         , intent(inout) :: molP
    type (dta_MOL)         , intent(in   ) :: PerM
    type (dta_ERR)         , intent(inout) :: econ
    double Precision :: Mlls(nLines,4)
    double Precision :: Aux_4M(4)
    integer*8        :: i, j, k
    integer*8        :: jBIG, jSMALL
    integer*8        :: indexI(nLines)
    Double Precision :: r_kj, rD0_kj
    Double Precision :: K1,K2, sumK2,sumK1
    !Auxiliar Constants
    double Precision :: Kaux, HWT, faH
    double precision :: T, Ptot, RT
    !Auxiliar for LAPACK
    !
    ! Parameters:
      integer*8, Parameter   :: N = 3
      integer*8, Parameter   :: NB = 10
      integer*8, Parameter   :: NRHS = 1
      integer*8, Parameter   :: LWMAX = 250
      integer*8              :: M
      integer*8              :: LDA, LDB
    ! Local Scalars:
      integer*8              :: LWORK, RANK
      integer*8              :: info1, info2
    ! Local Arrays:
      Double Precision       :: RCOND, RNORM, DNRM2
      Double Precision       :: SVal(N)
      Double Precision,ALLOCATABLE,DIMENSION(:)   :: WORK( : )
      Double Precision,ALLOCATABLE,DIMENSION(:,:) :: A( :, : ), B( :, : )
!-----------------------------------------
      T    = molP % Temp
      Ptot = molP % Ptot
      RT   = T0/T
!-----------------------------------------
      M = nLines
      LDA = M
      LDB = M
      ! LWORK >= 3*min(M,N) + max( 2*min(M,N), max(M,N), NRHS )
      ! LWORK >= 3*N + max( 2*N, M, NRHS )
      LWORK = 3*N + NB*( M+N )
      allocate ( WORK( LWMAX ) )
      allocate ( A( LDA, N ), B( LDB, NRHS ))
!
! LAPACK is used (installation command for mac):
! sudo port install lapack
! 
! * Compilation for "Free PGI compiler" for MAC:
! pgf90 myprog.f90 -llapack -lblas
!
! * Compilation "gfortran"
! gfortran myprog.f90 -llapack
!
!---------
! FIRST: create a zero matrix for 
    do j=1, nLines
      do k=1,4
        Mlls(j,k) = 0.0_dp
      enddo  
    enddo  
!
!---------
! Generate the Matrix for LLS:
    do j=1, nLines
      do k=1, nLines
        ! 
        if (j .eq. k) then
          faH = 1.0_dp
          if(T.ne.T0)faH = (RT**dta1%BHW(j)) 
          B(j,NRHS) = 2*molP%Ptot*dta1%HWT0(j)*faH
        else              
          if (isJb(dta1,j,k)) then
          ! CASE:  J(j) > J(k) (downwards transition j->k)
          ! or
          ! CASE: J(j) = J(k)
            jBIG   = j
            jSMALL = k
            r_kj = 1.0_dp 
          else
          ! CASE: J(j) < J(k)
          ! so downwards transition is (k->j)
          ! pj·<<k|W|j>> = pk·<<j|W|k>>; pk = dta1%PopuT(k); pj = dta1%PopuT(j)
            jBIG   = k
            jSMALL = j
            r_kj = dta1%PopuT(jBIG)/dta1%PopuT(jSMALL) !pjBIG/pjSMALL
          endif
          call LLS_Matrix(jBIG,jSMALL,dta1,molP,PerM,Aux_4M, econ)
          !
          rD0_kj = dta1%D0(k)/dta1%D0(j) 
          do i =1,4
            Mlls(j,i) = Mlls(j,i) + rD0_kj*r_kj*Aux_4M(i)
          enddo  
        endif    
      enddo
      !
    enddo
!
! **********************************************************************************
! LAPACK routine:
! ---------------
! DGELSS computes the minimum norm solution to a real linear least
! squares problem:
!
! Minimize 2-norm(| b - A*x |).
!
!===================================================================================
! Definition of A, B:
! A:
    A = -Mlls(1:nLines,1:3)
! B:
    do i = 1,nLines
      B(i,NRHS) = B(i,NRHS) + Mlls(i,4)
    enddo    
!
!  -- LAPACK driver routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
! Purpose
! =======
!
! DGELSS uses the singular value decomposition (SVD) of A. A is an M-by-N
! matrix which may be rank-deficient.
!
!  Several right hand side vectors b and solution vectors x can be
! handled in a single call; they are stored as the columns of the
! M-by-NRHS right hand side matrix B and the N-by-NRHS solution matrix
! X.
!
! The effective rank of A is determined by treating as zero those
! singular values which are less than RCOND times the largest singular
! value.
!
!  Arguments
!  =========
!
![in]  M 
!          M is INTEGER
!          The number of rows of the matrix A. M >= 0.
![in]  N 
!          N is INTEGER
!          The number of columns of the matrix A. N >= 0.
![in]  NRHS  
!          NRHS is INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrices B and X. NRHS >= 0.
![in]  A 
![out]     A is DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the M-by-N matrix A.
!          On exit, the first min(m,n) rows of A are overwritten with
!          its right singular vectors, stored rowwise.
![in]  LDA 
!          LDA is INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
![in]  B 
![out]     B is DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the M-by-NRHS right hand side matrix B.
!          On exit, B is overwritten by the N-by-NRHS solution
!          matrix X.  If m >= n and RANK = n, the residual
!          sum-of-squares for the solution in the i-th column is given
!          by the sum of squares of elements n+1:m in that column.
![in]  LDB 
!          LDB is INTEGER
!          The leading dimension of the array B. LDB >= max(1,max(M,N)).
![out] SVal 
!          SVal is DOUBLE PRECISION array, dimension (min(M,N))
!          The singular values of A in decreasing order.
!          The condition number of A in the 2-norm = S(1)/S(min(m,n)).
![in]  RCOND 
!          RCOND is DOUBLE PRECISION
!          RCOND is used to determine the effective rank of A.
!          Singular values S(i) <= RCOND*S(1) are treated as zero.
!          If RCOND < 0, machine precision is used instead.
![out] RANK  
!          RANK is INTEGER
!          The effective rank of A, i.e., the number of singular values
!          which are greater than RCOND*S(1).
![out] WORK  
!          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
![in]  LWORK 
!          LWORK is INTEGER
!          The dimension of the array WORK. LWORK >= 1, and also:
!          LWORK >= 3*min(M,N) + max( 2*min(M,N), max(M,N), NRHS )
!          For good performance, LWORK should generally be larger.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
![out] INFO  
!          INFO is INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  the algorithm for computing the SVD failed to converge;
!                if INFO = i, i off-diagonal elements of an intermediate
!                bidiagonal form did not converge to zero.
!
!  =====================================================================
!
!  Command sentence:
!  CALL DGELSS(M,N,1,A,LDA,B,M,S,RCOND,RANK,WORK,LWORK,INFO)
!
!  =====================================================================
!  Step0: Query the optimal workspace.
    if (econ % e(1) .ge. 1)  write(*,*) 'prior-Lwork', LWORK
    LWORK = -1
    info1 = 0
    CALL DGELSS(M,N,NRHS,A,LDA,B,M,SVal,RCOND,RANK,WORK,LWORK,info1)
    LWORK = INT( WORK(1) ) 
    if (LWORK .gt. LWMAX) then
      deallocate ( WORK )
      allocate ( WORK( LWORK ) )
    endif
!
!  Step1: Choose RCOND to reflect the relative accuracy of the input data
    !RCOND = 0.01D0 !<== These would be hard-imposed accuracy
    RCOND = -1.0D0  !<== < 0, machine precision is used instead.
!
!  Step2: Solve the least squares problem min( norm2(b - Ax) ) for the x
!         of minimum norm.
    info2 = 0
    CALL DGELSS(M,N,NRHS,A,LDA,B,M,SVal,RCOND,RANK,WORK,LWORK,info2)
!
! <---------------------------------------------
! NOTE: for further information about this subroutine 
! http://www.netlib.no/netlib/lapack/double/dgelss.f
!
! For examples:
! https://www.nag.com/lapack-ex/examples/source/dgelss-ex.f
!
! For information about LAPACK solving linear methods in general visit:
! http://www.netlib.org/lapack/lug/node26.html
!**********************************************************************************
! Least squares solution:
    molP%a1 = -B(1,NRHS)
    molP%a2 = -B(2,NRHS)
    molP%a3 = -B(3,NRHS)
    if (info2 .eq. 0) then
        if (econ % e(1) .ge. 1) write(*,*) 'SVD algorithm converged'
        !
        ! Print the effective rank of A
        !
        if (econ % e(1) .ge. 1) write(*,*) 'Tolerance used to estimate the rank of A', RCOND
        if (econ % e(1) .ge. 1) write(*,*) 'Estimated rank of A', RANK
    else
        call LLS_DGELSYSD_error(B(1,NRHS),B(2,NRHS),B(3,NRHS), info1, econ)
        if (econ % e(1) .ge. 1) write(*,*) 'The SVD algorithm failed to converge'
    endif

  END SUBROUTINE calc_QPar_DGELSS
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE calc_QPar_DGELSD(nLines, dta1, molP, PerM, econ)
!--------------------------------------------------------------------------------------------------------------------
! "calc_QPar_DGELSD": gives the values of a1, a2, a3 
! after solving linear system using the Lapack's DGELSD routine.
!
!---------------------------------------------------------------
    use module_common_var
    use module_error
    use module_maths
    use module_LLS
    Implicit None
    integer*8              , intent(in   ) :: nLines
    type (dta_SDF), pointer, intent(in   ) :: dta1
    type (dta_MOL)         , intent(inout) :: molP
    type (dta_MOL)         , intent(in   ) :: PerM
    type (dta_ERR)         , intent(inout) :: econ
    double Precision :: Mlls(nLines,4)
    double Precision :: Aux_4M(4)
    integer*8        :: i, j, k
    integer*8        :: jBIG, jSMALL
    integer*8        :: indexI(nLines)
    Double Precision :: r_kj, rD0_kj
    Double Precision :: K1,K2, sumK2,sumK1
    !Auxiliar Constants
    double Precision :: Kaux, HWT, faH
    double precision :: T, Ptot, RT
    !Auxiliar for LAPACK
    !
    ! Parameters:
      integer*8, Parameter   :: N = 3
      integer*8, Parameter   :: NB = 64
      integer*8, Parameter   :: NRHS = 1
      integer*8, Parameter   :: NLVL = 10
      integer*8, Parameter   :: LWMAX = 1000
      integer*8              :: M
      integer*8              :: LDA, LDB
    ! Local Scalars:
      integer*8              :: LWORK, LIWORK, RANK
      integer*8              :: info1, INFO
      integer*8,ALLOCATABLE,DIMENSION(:):: IWORK( : )
    ! Local Arrays:
      Double Precision       :: RCOND, RNORM, DNRM2
      Double Precision       :: SVal(N), sumD0, sumrjk
      Double Precision,ALLOCATABLE,DIMENSION(:)   :: WORK( : )
      Double Precision,ALLOCATABLE,DIMENSION(:,:) :: A( :, : ), B( :, : )
!-----------------------------------------
      T    = molP % Temp
      Ptot = molP % Ptot
      RT   = T0/T
!-----------------------------------------
      M = nLines
      LDA = M
      LDB = M
      LWORK = NB*( 2*M + N )
      LIWORK= 3*N*NLVL + 11*N +1
      allocate ( IWORK( LIWORK ) )
      allocate ( WORK( LWORK ) )
      allocate ( A( LDA, N ), B( LDB, NRHS ) )
!
! LAPACK is used (installation command for mac):
! sudo port install lapack
! 
! * Compilation for "Free PGI compiler" for MAC:
! pgf90 myprog.f90 -llapack -lblas
!
! * Compilation "gfortran"
! gfortran myprog.f90 -llapack
!
!---------
! FIRST: create a zero matrix for 
    do j=1, nLines
      do k=1,4
        Mlls(j,k) = 0.0_dp
      enddo  
    enddo  
!
!---------
! Generate the Matrix for LLS:
    do j=1, nLines
      sumrjk = 0.d0
      sumD0  = 0.d0
      sumK1  = 0.d0
      sumK2  = 0.d0
      do k=1, nLines
        ! 
        if (j .eq. k) then
          faH = 1.0_dp
          if(T.ne.T0)faH = (RT**dta1%BHW(j)) 
          B(j,NRHS) = 2*molP%Ptot*dta1%HWT0(j)*faH
        else              
          if (isJb(dta1,j,k)) then
          ! CASE:  J(j) > J(k) (downwards transition j->k)
          ! or
          ! CASE: J(j) = J(k)
            jBIG   = j
            jSMALL = k
            r_kj = 1.0_dp 
          else
          ! CASE: J(j) < J(k)
          ! so downwards transition is (k->j)
          ! pj·<<k|W|j>> = pk·<<j|W|k>>; pk = dta1%PopuT(k); pj = dta1%PopuT(j)
            jBIG   = k
            jSMALL = j
            r_kj = dta1%PopuT(jBIG)/dta1%PopuT(jSMALL) !pjBIG/pjSMALL
          endif
          call LLS_Matrix(jBIG,jSMALL,dta1,molP,PerM,Aux_4M, econ)
          !
          rD0_kj = dta1%D0(k)/dta1%D0(j) 
          do i =1,4
            Mlls(j,i) = Mlls(j,i) + rD0_kj*r_kj*Aux_4M(i)
          enddo 
          !ERASE ---->
          sumrjk = sumrjk + r_kj
          sumD0 = sumD0 + rD0_kj
          sumK1 = sumK1 + K1
          sumK2 = sumK2 + K2
          !<----- these lines
        endif    
      enddo
      !
    enddo
!
! **********************************************************************************
! LAPACK routine:
! ---------------
! DGELSD computes the minimum-norm solution to a real linear least
! squares problem:
!     minimize 2-norm(| b - A*x |)
!
! Definition of A, B:
    !
    A = -Mlls(1:nLines,1:3)
    !
    do i = 1,nLines
      B(i,NRHS) = B(i,NRHS) + Mlls(i,4)
    enddo    
!
!  -- LAPACK driver routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2017
!
!  Purpose
!  =======
!
! DGELSD computes the minimum-norm solution to a real linear least
! squares problem:
!     minimize 2-norm(| b - A*x |)
! using the singular value decomposition (SVD) of A with an algorithm 
! based on divide and conquer.
! A is an M-by-N matrix which may be rank-deficient.
!
! Note that The subroutine dGELSD is significantly faster than its 
! older counterpart dGELSS, especially for large problems, but may 
! require somewhat more workspace depending on the matrix dimensions.
!
! Several right hand side vectors b and solution vectors x can be
! handled in a single call; they are stored as the columns of the
! M-by-NRHS right hand side matrix B and the N-by-NRHS solution
! matrix X.
!
! The problem is solved in three steps:
! (1) Reduce the coefficient matrix A to bidiagonal form with
!     Householder transformations, reducing the original problem
!     into a "bidiagonal least squares problem" (BLS)
! (2) Solve the BLS using a divide and conquer approach.
! (3) Apply back all the Householder transformations to solve
!     the original least squares problem.
!
! The effective rank of A is determined by treating as zero those
! singular values which are less than RCOND times the largest singular
! value.
!
! The divide and conquer algorithm makes very mild assumptions about
! floating point arithmetic. It will work on machines with a guard
! digit in add/subtract, or on those binary machines without guard
! digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
! Cray-2. It could conceivably fail on hexadecimal or decimal machines
! without guard digits, but we know of none.
!
!  Arguments
!  =========
!
![in]  M 
!          M is INTEGER
!          The number of rows of A. M >= 0.
![in]  N 
!          N is INTEGER
!          The number of columns of A. N >= 0.
![in]  NRHS  
!          NRHS is INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrices B and X. NRHS >= 0.
![in,out]  A 
!          A is DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the M-by-N matrix A.
!          On exit, A has been destroyed.
![in]  LDA 
!          LDA is INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
![in,out]  B 
!          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the M-by-NRHS right hand side matrix B.
!          On exit, B is overwritten by the N-by-NRHS solution
!          matrix X.  If m >= n and RANK = n, the residual
!          sum-of-squares for the solution in the i-th column is given
!          by the sum of squares of elements n+1:m in that column.
![in]  LDB 
!          LDB is INTEGER
!          The leading dimension of the array B. LDB >= max(1,max(M,N)).
![out] S 
!          S is DOUBLE PRECISION array, dimension (min(M,N))
!          The singular values of A in decreasing order.
!          The condition number of A in the 2-norm = S(1)/S(min(m,n)).
![in]  RCOND 
!          RCOND is DOUBLE PRECISION
!          RCOND is used to determine the effective rank of A.
!          Singular values S(i) <= RCOND*S(1) are treated as zero.
!          If RCOND < 0, machine precision is used instead.
![out] RANK  
!          RANK is INTEGER
!          The effective rank of A, i.e., the number of singular values
!          which are greater than RCOND*S(1).
![out] WORK  
!          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
![in]  LWORK 
!          LWORK is INTEGER
!          The dimension of the array WORK. LWORK must be at least 1.
!          The exact minimum amount of workspace needed depends on M,
!          N and NRHS. As long as LWORK is at least
!              12*N + 2*N*SMLSIZ + 8*N*NLVL + N*NRHS + (SMLSIZ+1)**2,
!          if M is greater than or equal to N or
!              12*M + 2*M*SMLSIZ + 8*M*NLVL + M*NRHS + (SMLSIZ+1)**2,
!          if M is less than N, the code will execute correctly.
!          SMLSIZ is returned by ILAENV and is equal to the maximum
!          size of the subproblems at the bottom of the computation
!          tree (usually about 25), and
!             NLVL = MAX( 0, INT( LOG_2( MIN( M,N )/(SMLSIZ+1) ) ) + 1 )
!          For good performance, LWORK should generally be larger.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
![out] IWORK 
!          IWORK is INTEGER array, dimension (MAX(1,LIWORK))
!          LIWORK >= max(1, 3 * MINMN * NLVL + 11 * MINMN),
!          where MINMN = MIN( M,N ).
!          On exit, if INFO = 0, IWORK(1) returns the minimum LIWORK.
![out] INFO  
!          INFO is INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  the algorithm for computing the SVD failed to converge;
!                if INFO = i, i off-diagonal elements of an intermediate
!                bidiagonal form did not converge to zero.
!
!
!  =====================================================================
!
!  Command sentence:
!  CALL DGELSD(M , N , 1 , A , LDA , B , N , S , RCOND , RANK , WORK , LWORK , IWORK , INFO )
!
!  =====================================================================
!  Step0: Query the optimal workspace.
      LWORK = -1
      info1 = 0
      CALL DGELSD(M,N,NRHS,A,LDA,B,M,Sval,RCOND,RANK,WORK,LWORK,IWORK,info1)
      LWORK = INT( WORK(1) ) 
      if (LWORK .gt. LWMAX) then
        deallocate ( WORK )
        allocate ( WORK( LWORK ) )
      endif
      IF (info1.eq.0) LIWORK = INT(IWORK(1))
!  
!  Step1: Choose RCOND to reflect the relative accuracy of the input data
         !RCOND = 0.01D0 !<== These would be hard-imposed accuracy
         RCOND = -1.0D0  !<== < 0, machine precision is used instead.
!
!  Step2: Solve the least squares problem min( norm2(b - Ax) ) for the x
!         of minimum norm.
      INFO = 0
      RANK = 0
      CALL DGELSD(M,N,NRHS,A,LDA,B,LDB,Sval,RCOND,RANK,WORK,LWORK,IWORK,INFO)
!
! <---------------------------------------------
! NOTE: for further information about this subroutine 
! http://www.netlib.no/netlib/lapack/double/dgelsd.f
! http://www.netlib.org/lapack/explore-html/d7/d3b/group__double_g_esolve_ga94bd4a63a6dacf523e25ff617719f752.html#ga94bd4a63a6dacf523e25ff617719f752
!
! For examples:
! https://www.nag.com/lapack-ex/examples/source/dgelsd-ex.f
!
! For information about LAPACK solving linear methods in general visit:
! http://www.netlib.org/lapack/lug/node26.html
!**********************************************************************************
    molP%a1 = -B(1,NRHS)
    molP%a2 = -B(2,NRHS)
    molP%a3 = -B(3,NRHS)
    if (INFO .eq. 0) then
        if (econ % e(1) .ge. 1) write(*,*) 'SVD algorithm converged.'
        !
        ! Print the effective rank of A
        !
        if (econ % e(1) .ge. 1) write(*,*) 'Tolerance used to estimate the rank of A', RCOND
        if (econ % e(1) .ge. 1) write(*,*) 'Estimated RANK of A', RANK
    !
    else
        call LLS_DGELSYSD_error(B(1,NRHS),B(2,NRHS),B(3,NRHS), INFO, econ)
        if (econ % e(1) .ge. 1) write(*,*) 'The SVD algorithm failed to converge'
    endif

  END SUBROUTINE calc_QPar_DGELSD
!--------------------------------------------------------------------------------------------------------------------
