MODULE module_linemixing
!--------------------------------------------------------------------------------------------------------------------
! This module contains all subroutines related to molecular symmetry and behavior
!
    interface

        subroutine LM_Rosen(molP, nLines,dta1,Wmat,Y_RosT)
            use module_common_var
            use module_maths
            use module_phsub
            implicit none
            integer*8      ,intent (in)  :: nLines 
            Double Precision, intent(out) :: Wmat(nWmax,nWmax)
            double precision,intent(out)  :: Y_RosT(nLines)
            type (dta_SDF), intent(inout) :: dta1
            type (dta_MOL), intent(in)    :: molP
        end subroutine LM_Rosen

        subroutine WelCAL(dta1, nLines, molP, PerM, W_jk)
            use module_common_var
            use module_maths
            use module_phsub
            implicit none
            integer*8, intent(in)          :: nLines
            Double Precision, intent(out)   :: W_jk(nWmax,nWmax)
            type (dta_SDF), intent(in)      :: dta1
            type (dta_MOL), intent(in)      :: molP
            type (dta_MOL), intent(in)      :: PerM
        end subroutine WelCAL

        subroutine RN_Wmat(nLines, dta1, Wmat, W_rnO)
            use module_common_var
            use module_maths
            use module_phsub
            implicit none
            integer*8, intent(in)             :: nLines
            Double Precision, intent(in)    :: Wmat(nWmax,nWmax)
            Double Precision, intent(out)   :: W_rnO(nLines,nLines)
            type (dta_SDF), intent(in)      :: dta1
        end subroutine RN_Wmat

    end interface

END module module_linemixing
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
! Called By: 
! ---------
!
!
! T.Mendaza, last change 6 April 2016
!--------------------------------------------------------------------------------------------------------------------
!
    use module_common_var
    use module_maths
    use module_phsub
    implicit none
    integer*8      ,intent (in)  :: nLines 
    Double Precision, intent(out) :: Wmat(nWmax,nWmax)
    double precision,intent(out)  :: Y_RosT(nLines)
    type (dta_SDF), intent(inout) :: dta1
    type (dta_MOL), intent(in)    :: molP
    !Local Var
    integer*8                    :: i,j,k
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
         DipoI= abs(dsqrt(dta1%Str(i)/(dta1%Sig(i)*dta1%PopuT0(i))))
         do k=1,nLines
            DipoK= abs(dsqrt(dta1%Str(k)/(dta1%Sig(k)*dta1%PopuT0(k))))
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
            !
            !sumY=sumY + (dabs(dta1%DipoT0(k)) &
            !    /dabs(dta1%DipoT0(i)))* &
            !    ( Wmat(i,k)/delta ) ! Wmat(k,i)
            sumY=sumY + (DipoK/DipoI)* &
                ( Wmat(i,k)/delta ) ! Wmat(k,i)
        enddo
        Y_RosT(i)=2.0*Ptot*sumY
    ENDDO

    Return
  END SUBROUTINE LM_Rosen
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE WelCAL(dta1, nLines, molP, PerM, W_jk)
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
! Called Routines: "PFmol"    (Partition Function of CH4)
! ---------------  "Readline" (READ LINE data and relaxation matrix)
!     
! Called By: "main program"  
! ---------
!     
!     Double Precision Version
!     
! T. Mendaza last change 08 Jan 2016
!--------------------------------------------------------------------------------------------------------------------
!
    use module_common_var
    use module_phsub
    use module_maths
    Implicit None
    integer*8,        intent(in)      :: nLines
    type (dta_SDF), intent(in)      :: dta1
    type (dta_MOL), intent(in)      :: molP
    type (dta_MOL), intent(in)      :: PerM
    double precision, intent(out)   :: W_jk(nWmax,nWmax)
    !double precision, intent(out)   :: W_jk(nWmax,nWmax)
    integer*8                      :: i, j, k, n
    integer*8                      :: jBIG, jSMALL
    integer*8                      :: indexI(nLines)
    Double Precision                :: r_kj, pn, pk
    Double Precision                :: auxW, auxHW
    !Double Precision                :: Wtest(nLines,nLines)
    !Auxiliar Constants
    integer*8                      :: count1, count2
    double Precision                :: Kaux, HWT, faH
    double Precision                :: T, Ptot, RT
!
!-----------------------------------------
    T    = molP % Temp
    Ptot = molP % Ptot
    RT   = T0/T
!-----------------------------------------
!
! FIRST: create a zero matrix for 
    do j=1, nLines
      do k=1,nLines
        W_jk(j,k) = 0.0_dp
      enddo  
    enddo  
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
        !!! + (xH2O*(HWT0_H2O(i)*(RT**BHW_H2O(i))))
        ! where:
        ! HWT0_H2O = halfwidth broadeing by water vapour (at 296K)
        ! xH2O  = H2O atmospheric percentage 
        ! BHW_H2O = water vapour temperature exponent
        !
          W_jk(j,k) = 2.0*Ptot*HWT !+ i*(-0.008)  
          !stop
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
          if ((W_jk(j,k) .eq. 0.0_dp) .or. (W_jk(k,j) .eq. 0.0_dp)) then 
              if (molP%M .eq. 7) then
                ! O2
                W_jk(jBIG,jSMALL) = K_jkO2(jBIG,jSMALL,dta1,nLines,molP,PerM) 
              else 
                W_jk(jBIG,jSMALL) = K_jkCalc(jBIG,jSMALL,dta1,nLines,molP,PerM)  
              endif
              !stop      
              if ( isnan( W_jk(jBIG,jSMALL) ) .or. isinf( W_jk(jBIG,jSMALL) ) ) then 
                W_jk(jBIG,jSMALL) = 0.0_dp  
                print*, "Wij NaN!", jBIG,jSMALL
                stop    
              endif                          
              ! so downwards transition is (k->j)
              ! so downwards transition is (jBIG -> jSMALL)
              ! pjSMALL·<<jBIG|W|jSMALL>> = pjBIG·<<jSMALL|W|jBIG>>; 
              ! where:
              ! pjBIG = dta1%PopuT(jBIG); pjSMALL = dta1%PopuT(jSMALL)
              !
              r_kj = dta1%PopuT(jBIG)/dta1%PopuT(jSMALL) !pjBIG/pjSMALL
              W_jk(jSMALL,jBIG) = r_kj*W_jk(jBIG,jSMALL)  
              !
              !print*, "downwards transition"
              !write(*,'(a7,i3,a3,i3,a2,2x,E12.4)'),'W_jk(j=',jBIG,';k=',jSMALL,')', W_jk(jBIG,jSMALL) 
              !print*, "upwards transition" 
              !write(*,'(a7,i3,a3,i3,a2,2x,E12.4)'),'W_jk(j=',jSMALL,';k=',jBIG,')', W_jk(jSMALL,jBIG)
          else
            !next transition
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
    !CALL sumRule(nLines,indexI,dta1%D0(1:nLines),Wtest,0.5)
    !
    Return
  END SUBROUTINE WelCAL
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE RN_Wmat(nLines, dta1, Wmat, W_rnO)
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
! Called Routines: "isnan" 
! ---------------  "isinf" 
!     
! Called By: "main program"
! ---------
!     
!     Double Precision Version
!     
! T. Mendaza last change 01 Abr 2016
!--------------------------------------------------------------------------------------------------------------------
!
    use module_common_var
    use module_maths
    Implicit None
    integer*8,        intent(in)    :: nLines
    Double Precision, intent(in)  :: Wmat(nWmax,nWmax)
    Double Precision, intent(out) :: W_rnO(nWmax,nWmax)
    type (dta_SDF), intent(in)    :: dta1
    !Local variables
    integer*8                    :: i, j, k, n
    integer*8                    :: indexS(nLines),indexI(nLines)
    double Precision              :: str(nLines)
    Double Precision              :: Sup, Slow, S_UL_rate
    Double Precision              :: pn, pk
    Double Precision              :: auxW, auxHW
    Double Precision              :: W_rn(nLines,nLines)
    !Other Auxiliar Constants
    integer*8                    :: count1, count2
    double Precision              :: Kaux, Saux, Waux
    double Precision              :: Supaux, Slowaux
!---------
! Reordering by line strength
!
       do k=1,nLines
         !str(k) = dta1%STR(k)
         !str(k) = dta1%sig(k)*dta1%popuT(k)*dta1%D0(k)**2 
         !str(k) = dta1%sig(k)*dta1%popuT(k)*dta1%DipoT0(k)**2
         str(k) = dta1%popuT(k) 
       enddo
       ! The rotational components are sorted according to their intensities in decreasing order
       ! (the first one is the most intense). Hence, << 1 | W | 1 >> is the broadening of the most 
       ! intense line, and << 2 | W | 1 >> (or << 1 | W | 2 >>) are the terms coupling this transition
       ! to the second most intense one.
       ! SORTENING by intensity.
       call bubble_index(nLines,str,indexS,'d')
       ! Here we perform the 'pullback' of indexS. 
       call ibubble_index(nLines,indexS,indexI,'a')
       !do k=1,nLines
       !  write(*,'(a5,E10.2,a5,F7.3,a5,E10.2,a5,F4.0)'), &
       !  'str',str(indexS(k)),'D0',dta1%D0(indexS(k)),'Po',dta1%PopuT(indexS(k)),&
       !  'J',dta1%J(indexS(k),1)
       !enddo
       do i=1,nLines
          do j=1,nLines
            if (i.eq.j) then
              W_rn(i,i)=Wmat(indexS(i), indexS(i))
            else
              W_rn(i,j)=-dabs(Wmat(indexS(i), indexS(j)))
            endif
          enddo
       enddo
       !stop
       !
       ! Then, each column "k" of the matrix, starting from the first one is treated as follows:
       DO n = 1, nLines
         Sup = 0.0; Slow = 0.0
         count1 = 0; count2 = 0
         pn = dta1%popuT( indexS(n) )
         DO k = 1, nLines
          Saux = 0.0
          ! remove assymetry
          if(dta1%iso.gt.2 .AND. dta1%iso.ne.7 &
            .AND.&           
            mod(abs(int(dta1%J(n,1))-int(dta1%J(k,1))),2).ne.0)cycle
          if (k .le. n) then 
            Sup  = Sup  + dabs( dta1%D0( indexS(k) ) )* W_rn( n , k )
            print*, "sup",dabs( dta1%D0( indexS(k) ) ),W_rn( n , k )
            if(isnan(Sup).or.isinf(Sup))stop
            !write(*,'(a2,i3,a1,i3,a2,E10.2)'),"W_rn(",n,",",k,")=", W_rn(n,k) 
            !write(*,'(a4,E10.2)'), "Sup=",Sup
          else !if (k .gt. n) then
            Slow = Slow + dabs(dta1%D0( indexS(k) )) * W_rn( n , k ) 
            print*, "slo",dabs( dta1%D0( indexS(k) ) ),W_rn( n , k )
            if(isnan(Slow).or.isinf(Slow))stop
            !write(*,'(a2,i3,a1,i3,a2,E10.2)'), "W_rn(",n,",",k,")=" , W_rn( n , k ) 
            !write(*,'(a4,E10.2)'), "Slow=",Slow          
          endif
         ENDDO
         S_UL_rate=Sup/Slow
         !print*, S_UL_rate
         !print*, "counting: Sup NaN#", count1, "; Slow NaN#", count2, " out of", nLines
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
              ! K(n,k) = << k | W | n >> 
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
              !print*, W_rn(n,k), W_rn(k,n) 
            endif
            !stop
          endif
         ENDDO  
         !
       ENDDO
       !
       ! use '1.0' if Wii = line-width 
       ! use '2.0' if Wii = half-width 
       CALL sumRule(nLines,indexS,dta1%D0(1:nLines),W_rn,1.0)
       ! 
       ! Reordered by wno
       !
       do i=1,nLines
         do j=1,nLines
           W_rnO(i,j)  =  W_rn( indexI(i) , indexI(j) )
         enddo
       enddo
  END SUBROUTINE RN_Wmat
!--------------------------------------------------------------------------------------------------------------------