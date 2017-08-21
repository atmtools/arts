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

        logical function rule2(P,nLines,W,v)
            use module_common_var
            implicit none
            integer*8       , intent(in) :: nLines
            double precision, intent(in) :: P
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
! Subroutine to obtain the relaxation matrix elements.
! 
! Variables:
!
! Input/Output Parameters of Routine (Arguments or Common)
! --------------------------------------------------------
! nLines  : Number of lines of A band (Input).
! HWT0    : Air (or H2) Broadened HalfWidths of the Lines for
!           T0 and P0 (Cm-1)
! PopuT0  : Populations of the Lower Levels of the Lines
!           at Temperature T0
!
! Other important Output Quantities (through "module_common_var")
! ---------------------------------
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
    !Double Precision                :: Wtest(nLines,nLines)
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
        !write(*,'(a3,i5,a3,i5,a1)'),'(j=',j,';k=',k,')'         
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
                W_jk(jBIG,jSMALL) = K_jkCalc(jBIG,jSMALL,dta1,nLines,molP,PerM,econ)  
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
          
        endif    
      enddo

    enddo
    !
    ! Sum Rule TEST
    !
    ! The Relaxation matrix must follow the SUM-RULE
    ! However, that is acomplish regarding that 
    ! every relaxation betwen every theoretical transition
    ! (not just empirically observed) is taken into the sum.
    ! Since we are doing the matching with HITRAN from the very 
    ! begining the test will fail, but every sum value has to be
    ! under the linewidth (== 2*HalfWidth which is the value 
    ! provided by HITRAN).
    !
    !print*, "sum rule?"
    !do i = 1, nLines
    !  indexI(i) = i
    !enddo
    !do i= 1, nLines
    !  do j = 1,nLines
    !    Wtest(i,j) = W_jk(i,j)
    !  enddo
    !enddo
    !CALL sumRule(nLines,indexI,dta1%D0(1:nLines),Wtest,0.5,econ,tOK)
    
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
! 
! Variables:
!
! Input/Output Parameters of Routine (Arguments or Common)
! --------------------------------------------------------
! nLines  : Number of lines of A band (Input).
! HWT0    : Air (or H2) Broadened HalfWidths of the Lines for
!           T0 and P0 (Cm-1)
! PopuT   : Populations of the Lower Levels of the Lines
!           at Temperature T
!
! Other important Output Quantities (through "module_common_var")
! ---------------------------------
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
! T. Mendaza last change 01 Abr 2016
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
    double Precision              :: str(nLines)
    Double Precision              :: Sup, Slow, S_UL_rate
    Double Precision              :: pn, pk
    Double Precision              :: W_rn(nLines,nLines)
!---------
! Reordering by line strength
!
      do k=1,nLines
         !str(k) = dta1%STR(k)
         !str(k) = dta1%sig(k)*dta1%popuT(k)*dta1%D0(k)**2 
         !str(k) = dta1%sig(k)*dta1%popuT(k)*dta1%DipoT(k)**2
         str(k) = dta1%popuT(k) 
      enddo
      ! The rotational components are sorted according to their intensities in decreasing order
      ! (the first one is the most intense). Hence, << 1 | W | 1 >> is the broadening of the most 
      ! intense line, and << 2 | W | 1 >> (or << 1 | W | 2 >>) are the terms coupling this transition
      ! to the second most intense one.
      ! SORTENING by intensity.
      call bubble_index(nLines,str,indexS,'d',econ)
      ! Here we perform the 'pullback' of indexS. 
      call ibubble_index(nLines,indexS,indexI,'a',econ)
      !
      do i=1,nLines
          do j=1,nLines
            if (i.eq.j) then
              W_rn(i,i)=Wmat(indexS(i), indexS(i))
            else
              !W_rn(i,j)=-abs(Wmat(indexS(i), indexS(j)))
              W_rn(i,j)=Wmat(indexS(i), indexS(j))
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
      ! use '1.0' if Wii = line-width 
      ! use '2.0' if Wii = half-width 
      e20 = econ % e(2)
      CALL sumRule(nLines,indexS,dta1%D0(1:nLines),W_rn,1.0,econ)
      ! 
      ! Reordered by wno
      !
      if (( econ % e(2) .gt. e20 ) .and. ((econ % e(1) .eq. -1).or.(econ % e(1) .eq. 2) )) then
        econ % e(2) = econ % e(2) - 1
        CALL just_fill_DiagWRn(nLines,dta1 % BHW, dta1 % HWT0, T/T0, P,W_rnO)
        !do i = 1, nLines
        !  print*,'wno', dta1 % sig(i) 
        !  print*,'Str', dta1 % Str(i)
        !  print*,'HWT', dta1 % HWT0(i)
        !  print*,'BHW', dta1 % BHW(i) 
        !  print*,'E  ', dta1 % E(i)
        !  print*,'g0 ', dta1 % swei0(i) 
        !  print*,'g00', dta1 % swei00(i)
        !  print*,'UP ', dta1%J(i,1),dta1%N(i,1),dta1%espin(i,1)
        !  print*,'LO ', dta1%J(i,2),dta1%N(i,2),dta1%espin(i,2)
        !  print*,'bra', dta1%br(i)
        !enddo
      else 
        do i=1,nLines
          do j=1,nLines
            W_rnO(i,j)  =  W_rn( indexI(i) , indexI(j) )
          enddo
        enddo
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
!
! Accessed Files: 
! --------------
!
! Called Routines: 
! ---------------  
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
!     Build the Ym from the W
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
! Y2 : g_k, second coeffcient from the second order linemixing approx [Smith, 1981]. 
! Y3 :dv_k, is the second order shift [Smith, 1981].
! NOTE: on the calculation of these coeff (as well the 1st order one), 
! the first order pressure shift has been neglected, as it was done by Smith (1981).
!
! Accessed Files: 
! --------------
!
! Called Routines: 
! ---------------  
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
!
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
!
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
!  CALL dgelsd(TRANS, M , N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK, LWORK, IWORK, INFO ) 
!
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
!
! Generate the Matrix for LLS:
print*, "Generate the Matrix for LLS"
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
    print*, "LLS_AF_Matrix done"
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
!
! Definition of A, B:
    !
    A = Mlls(1:nLines,1:N)
    !
    do i = 1,nLines
      B(i,NRHS) = B(i,NRHS) - Mlls(i,N+1)
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
!  CALL dgelsd(TRANS, M , N, NRHS, A, LDA, B, LDB, S, RCOND, RANK, WORK, LWORK, IWORK, INFO ) 
!
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
! Detailed description:
! ---------------------
! 
! Variables:
!
! Input/Output Parameters of Routine (Arguments or Common)
! --------------------------------------------------------
! --      : .
!
! Other important Output Quantities 
! ---------------------------------
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
  logical function rule2(P,nLines,W,v)
!--------------------------------------------------------------------------------------------------------------------
! rule2: Perturbation Theory limitation.
!
! Detailed description:
! ---------------------
! To invert the MAtrix expresion, th eperturbation theory is appropriate under the condition
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
!       is an arbitrary selection
!     
! T. Mendaza last change 20 Feb 2017
! --------------------------------------------------------
    use module_common_var
    implicit none
    integer*8       , intent(in) :: nLines
    double precision, intent(in) :: P
    double precision, intent(in) :: W(nLines,nLines)
    double precision, intent(in) :: v(nLines) 
    integer                      :: l,k,pos(2)
    double precision,parameter   :: TOL_rule2 = 0.1_dp
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
            !print*, v(l),v(k),W(l,k)
          endif
        endif
      ENDDO
    ENDDO

    RETURN
  end function rule2
!--------------------------------------------------------------------------------------------------------------------