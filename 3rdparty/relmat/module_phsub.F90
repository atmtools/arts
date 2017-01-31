!--------------------------------------------------------------------------------------------------------------------
module module_phsub
!--------------------------------------------------------------------------------------------------------------------
! This module contains all subroutines related to 
!
    interface

        subroutine DipCAL(dta1,nLines,molP)
            use module_common_var
            use module_molecSp
            use module_maths
            use module_error
            implicit none
            integer*8, intent(in)             :: nLines
            type  (dta_SDF), intent(inout)    :: dta1
            type  (dta_MOL), intent(in)       :: molP
        end subroutine DipCAL
        
        double precision function pureHund(caseHund, J_i, N_i, J_f, N_f)
            implicit none
            integer (kind=8),intent(in)       :: N_i, N_f
            real*8          , intent(in)      :: J_i, J_f
            character     , intent(in)        :: caseHund
        end function pureHund
        
        subroutine PopuCAL(dta1, nLines, molP)
            use module_common_var
            implicit none
            integer*8,intent(in)            :: nLines
            type (dta_MOL), intent(in)      :: molP
            type (dta_SDF), intent(inout)   :: dta1
        end subroutine PopuCAL
        
        double precision function K_jkCalc(j,k,h,nLines,molP,PerM) 
            use module_common_var
            use module_molecSp
            use module_maths
            implicit none
            integer*8, intent(in)             :: j,k,nLines
            type (dta_SDF) , intent(in)       :: h
            type (dta_MOL) , intent(in)       :: molP, PerM
        end function K_jkCalc

        double precision function K_jkO2(j,k,h,nLines,molP,PerM)
            use module_common_var
            use module_molecSp
            use module_maths
            implicit none
            integer*8, intent(in)             :: j,k,nLines
            type (dta_SDF), intent(in)        :: h
            type (dta_MOL), intent(in)        :: molP, PerM

        end function K_jkO2
        
        double precision function Ql_mol_X(molP,L)
            use module_common_var
            use module_molecSp
            use module_maths
            implicit none
            integer*8, intent(in)             :: L
            type(dta_MOL)  , intent(in)       :: molP
        end function Ql_mol_X
        
        double precision function AFmol_X(molP, PerM, L, step)
            use module_common_var
            implicit none
            integer*8, intent(in)             :: step
            real*8, intent(in)                :: L
            type (dta_MOL), intent(in)        :: molP, PerM
        end function AFmol_X

        subroutine sumRule(nLines,indexS,dipole,Wmat,dfact)
            use module_error
            implicit none
            integer (kind=8), intent(in)    :: nLines
            integer (kind=8), intent(in)    :: indexS(nLines)
            real*8, intent(in)              :: dfact
            Double Precision, intent(in)    :: dipole(nLines)
            Double Precision, intent(in)    :: Wmat(nLines,nLines)
        end subroutine sumRule

        subroutine VarInit(dta1,dta2,molP)
            use module_common_var
            implicit none
            type (dta_SDF), intent(inout)    :: dta1
            type (dta_RMF), intent(inout)    :: dta2
            type (dta_MOL), intent(inout)    :: molP
        end subroutine VarInit

    end interface

END module module_phsub
!--------------------------------------------------------------------------------------------------------------------
! RELAXATION MATRIX SUBROUTINES -------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE DipCAL(dta1,nLines,molP)
!--------------------------------------------------------------------------------------------------------------------
! "DipCAL": DIPole CALculation
! 
! Detailed Description:
! ---------------------
! This subroutine computes the dipole element d_k for a certain level 
! obtaining parameters from Subroutine "Readline". 
!
! Variables:
!
! Input/Output Parameters of Routine (Arguments or Common)
! ----------------------------------
!	E         : lower state energy (cm-1)
!	Sig       : WaveNumber of the Computation (cm-1, Input)
!	PopuT0    : Populations of the Lower state of the Lines
!	            at Temperature T0 (input)
!
! Other important Input Constant (through "module_common_var")
! ------------------------------
!	T0	      : Reference Temperature in Kelvin, 296K (Input).
!	c2        : second radiation constant = hc/k = 1.4388 cm·K
!pure_b_case: --
!
! Output Quantities (through Common Statements)
! -----------------
!	DipoT0	  : Dipole transition Moments of the Line at T0.
!
!
! Accessed Files:  None
! --------------
!
! Called Routines: 'PopuCAL' (POPUlation of the lower state CALculation)
! ---------------  
!
! Called By: Main Program
! ---------
!
!	Double Precision Version
!
!	T.Mendaza, last change 11 March 2016
!--------------------------------------------------------------------------------------------------------------------
!
    use module_common_var
    use module_molecSp
    use module_maths
    use module_error
    Implicit None
    integer*8, intent(in)       :: nLines
    type (dta_SDF),intent(inout):: dta1
    type (dta_MOL),intent(in)   :: molP
    integer*8                   :: N_i, N_f, J_i, J_f
    integer*8                   :: i, j, k
    real*8                      :: s_i, s_f
    double precision            :: cte,fA,fB, w3j
    double precision            :: pure_b_case, purehund, d0
    double precision            :: signDIP, Ia_hit
!
!----------
!
! Check that Arrays for Results are Large Enough.
       if ( nLines .gt. nLmx) then
        call sizeError("1000", nLines, nLmx)
       endif
!
! Compute the Dipole moment at every line 
! (Diatomic approximation... 
!
! INIT. VAR:
!-----------
    DO k = 1, nLines
    !print*, k, nLines
      J_i  = int(dta1%J(k,1)) !J_i  = J LOWER L.
      J_f  = int(dta1%J(k,2)) !J_f  = J UPPER L.
      !
      ! REDUCE MATRIX ELEMENT, D0
      !
        if ( mod((J_f+dta1%lv2(2)),2) .eq. 0) then
        ! j = (-1)**(J_f+lf)
          j = 1
        else
          j = -1
        endif
        ! K_t = marks the dimension of the tensor.
        ! isotropic Raman = 0 
        !   IR ----> K_t = 1
        !   Raman -> K_t = 2
        ! wigner simbols are implemented for -> IR case.
        w3j = wigner3j( dta1%J(k,1), real(K_t,dp), dta1%J(k,2), &
                        real(dta1%lv2(1), dp), real( dta1%lv2(2)-dta1%lv2(1),dp), &
                        real(-dta1%lv2(2), dp))
        d0 = j*dsqrt(2.*dta1%J(k,2) + 1.)*w3j
        dta1%D0(k)    = d0
      ! print*,dta1%D0(k)
!      !
!      ! RIGID ROTOR DIPOLE
!      !
!        if ( mod(J_f,2) .eq. 0) then
!        ! j = (-1)**(J_f)
!          j = 1
!        else
!          j = -1
!        endif
!        w3j = wigner3j( dta1%J(k,1), real(K_t,dp), dta1%J(k,2), &
!                        0.0_dp, 0.0_dp, 0.0_dp)
!        dta1%Drigrotor(k) = j*dsqrt(2.*dta1%J(k,2) + 1.)*w3j
!      !
!      ! DIPOLE MOMENT 
!      !
!      ! Weighted transition moment Squared (|R)
!      ! FORMULA:
!      ! |R = (1/g00)*|R12|^2 ; 
!      !
!      ! where:
!      ! * R12 is the transition moment 
!      !     1) Electric dipole transition (edt)
!      !         g2·A_21 = (16·pi^3/3·h·ε0)·wno^3·|R12|^2.
!      !     units: C2·m2
!      !
!      !     2) magnetic dipole transition (mdt)
!      !         g2·A_21 = (16·µ0·pi^3/3·h)·wno^3·|R12|^2.
!      !     units: A2·m4
!      !
!      !     3) electric-quadrupole transitions (eqt)
!      !         g2·A_21 = (8·pi^5/5·h·ε0)·wno^5·|R12|^2.
!      !     units: C2·m4
!      ! * ε0 is the vacuum permittivity fundamental physical constant
!      ! * µ0 is the vacuum permeability fundamental physical constant
!      !
!      ! NOTE: for conversion to centimetre–gram–second system of units (cgs) 
!      !       electrostatic units, ε0 should be replaced by 1/4pi.
!      ! See:
!      !* Tatum JB. The interpretation of intensities in diatomic molecular spectra. 
!      !  Astrophys J Suppl Ser 1967;14:21–56;
!      !* Tatum JB. Erratum: interpretation of intensities in diatomic molecular spectra. 
!      !  Astrophys J Suppl Ser 1971;22:388.
!      !
!      ! NOTE: the majority of the transitions in HITRAN are 
!      ! electric-dipole.
!      !
!      ! NOTE 2: 
!      ! To calculate the dipole moment you can use the Einstein 
!      ! coefficients or the line intensity.
!      ! We have followed [Simenckova et al. 2006] for this fomulation.
!      !
!      ! UNITS:
!      ! [Debye^2 = 10E-36 ergs·cm3] where [1erg = 1E-07 J]
!      !
!      if (tdcal .eq. 'S') then
!      ! Using LINE-INTENSITY:
!      !intensity factor in HITRAN:
!        Ia_hit = molP % IAb(dta1%iso)
!        !print*, Ia_hit, molP % IAb(dta1%iso)
!        !stop
!      ! Parts of the formula:
!        fA  = -c2*dta1%Sig(k)/T0
!        fB  = 1. - exp(fA)
!        cte = 1.0E036*3.d0*hp_cgs*c/(8.d0*Pi**3)
!        dta1%DipoT0(k)= cte*v_permit*dta1%Str(k)/ &
!                  (Ia_hit*dta1%Sig(k)*&
!                  dta1%PopuT0(k)*fB   & !-> ergs·cm3
!                  )!*1E-07! => J·cm3
!        !
!        ! DIPOLE MOMENT (R12)
!        !
!        dta1%DipoT0(k)= DSQRT(dta1%DipoT0(k)*dta1%swei00(k))
!        
!       else if (tdcal .eq. 'A') then
!       ! THIS PART SHOULD BE CHECKED (UNITS)
!            if (tmt .eq. 'edt') then
!                ! because hp [J·s] and J = Kg m2/s2 = 1E7 g·cm2/s2
!                ! and our v_permit is in c·g·s?
!                ! then we have to transform J into that system?
!                ! conversion C2m2 to J·cm3 ??
!                cte = 3.*1E07*hplank*v_permit/(16.*pi**3)
!                dta1%DipoT0(k) = cte*dta1%A21(k)*&
!                                (dta1%swei0(k)/&
!                                 dta1%Sig(k)**3)
!            else if (tmt .eq. 'mdt') then
!                cte = 3.*hplank/(16.*v_permea*pi**3)
!                dta1%DipoT0(k) = dta1%A21(k)*&
!                                (dta1%swei0(k)/&
!                                 dta1%Sig(k)**3)
!            else if (tmt .eq. 'eqt') then
!                cte = 5.*1E07*hplank*v_permit/(8.*pi**5)
!                dta1%DipoT0(k) = dta1%A21(k)*&
!                                (dta1%swei0(k)/&
!                                 dta1%Sig(k)**5)
!            else
!                print*, 'No transition moment type selected.'
!            endif
!       else
!        print*, 'No dipole moment calc. selected.'
!       endif
!      !
!      ! Dipole sign
!      ! 
!      if ( molP%M .eq. 7 ) then
!        ! if the program is dealing with a diatomic molecule
!        ! then HUND's pure cases help to the calcualtion of 
!        ! the sign.
!        ! DATA:
!        ! N_i  = dta1%N(k,1)
!        ! N_f  = dta1%N(k,2)
!        ! J_i  = dta1%J(k,1) !LOWER L.
!        ! J_f  = dta1%J(k,2) !UPPER L.
!        ! Call the function:
!        !pure_b_case= pureHund(caseHund,      J_i,         N_i,         J_f,         N_f)
!        pure_b_case = pureHund(caseHund, dta1%J(k,1), dta1%N(k,1), dta1%J(k,2), dta1%N(k,2))
!        signDIP =pure_b_case/ABS(pure_b_case)    
!        if (isnan(signDIP) .or. isInf(signDIP)) then
!          signDIP = dta1%DipoT0(k-1)/abs(dta1%DipoT0(k-1))
!        endif
!      else
!        signDIP =d0/ABS(d0) 
!        if (isnan(signDIP) .or. isinf(signDIP)) then
!          signDIP = dta1%DipoT0(k-1)/abs(dta1%DipoT0(k-1))
!        endif
!      endif
!      dta1%DipoT0(k)= dta1%DipoT0(k)*signDIP
!      !print*, dta1%Qlow(k), dta1%D0(k), dta1%popuT0(k)
!      !print*, "----------------------------------"
!      !if (k .eq. 4) stop    
    ENDDO
! 
      Return
  END SUBROUTINE DipCAL
!--------------------------------------------------------------------------------------------------------------------
  double precision function pureHund(caseHund, J_i, N_i, J_f, N_f)
!--------------------------------------------------------------------------------------------------------------------
! " pureHund " = Pure Hund's case
!
! Detailed description:
! ---------------------
! In rotational-vibrational and electronic spectroscopy of DIATOMIC MOLECULES, 
! Hund's coupling cases are idealized cases where specific terms appearing in 
! the molecular Hamiltonian and involving couplings between angular momenta 
! are assumed to dominate over all other terms. There are five cases, 
! traditionally notated with the letters (a) through (e). 
! Most diatomic molecules are somewhere between the idealized cases (a) and (b).
!
! Variables:
! J_i	  = Resultant total angular momentum quantum number L + S, 
!         excluding nuclear spins.(Subindex "i" = initial state) 
! J_f	  = Resultant total angular momentum quantum number L + S,
!         excluding nuclear spins. (Subindex "f" =final   state)
! N_i/f	= the total angular momentum minus the electron spin J - S (i/f)
!         Rotational angular momentum quantum number, 
!         excluding electron and nuclear spins, 
!         in the case where electron spin is present.
! auxiliary var:
! L	= the electronic orbital angular momentum
! S	= the electronic spin angular momentum
!
! Hund's case	Electrostatic	Spin-orbit	    Rotational
! -----------------------------------------------------------
! (a)		    strong	    	intermediate	   weak
! (b)		    strong		    weak		         intermediate
! (c)		    intermediate	strong		       weak
! (d)		    intermediate	weak		         strong
! (e)		    weak		      intermediate	   strong
!				    strong		    intermediate
!
      IMPLICIT NONE
      integer  , parameter         :: dp = kind(1.0D0) !double precission
      integer (kind=8) , intent(in):: N_i, N_f
      double precision, intent(in) :: J_i, J_f
      character, intent(in)        :: caseHund
      double precision             :: cte1
! -----------------
    if (caseHund .eq. 'a') then
        !write(*,'(a10)'), 'Case (a)'
        pureHund = 1.0_dp
    else if (caseHund .eq. 'b') then
        !write(*,'(a10)'), 'Case (b)'
        if ( (N_f .eq. N_i+1 ) .and. ( J_f .eq. J_i+1.0_dp ) ) then
         cte1 = dsqrt(J_i)
      else if ( (N_f .eq. N_i+1 ) .and. ( J_f .eq. J_i ) ) then
         cte1 = (-1.0_dp)*dsqrt(J_i + 1.0_dp)
      else if ( (N_f .eq. N_i-1 ) .and. ( J_f .eq. J_i ) ) then
         cte1 = (-1.0_dp)*dsqrt(J_i)
      else if ( (N_f .eq. N_i-1 ) .and. ( J_f .eq. J_i-1.0_dp ) ) then
         cte1 = dsqrt(J_i + 1.0_dp)
      else ! N_f = N_i !Q band
         cte1 = 1.0_dp
      endif
      pureHund = dsqrt(1.0_dp/(2.0_dp*J_i +1.0_dp)) * cte1
        !print*, cte1, sqrt(1./(2.*J_i +1.)), pureHund
    else if (caseHund .eq. 'c') then
        !write(*,'(a10)'), 'Case (c)'
      pureHund = 1.0_dp
    else if (caseHund .eq. 'd') then
        !write(*,'(a10)'), 'Case (d)'
      pureHund = 1.0_dp
    else if (caseHund .eq. 'e') then
        !write(*,'(a10)'), 'Case (e)'
      pureHund = 1.0_dp
    else
        !write(*,'(a100)'), 'Case selected does not match with any (a,b,c,d,e)'
        pureHund = 1.0_dp
    endif      
    RETURN
  END function pureHund
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE PopuCAL(dta1, nLines, molP)
!--------------------------------------------------------------------------------------------------------------------
!"PopuCAL": POPUlation CALculation
!
! Detailed Description:
! ---------------------
! It calcules the Populations of the Lower state of the 
! Line/transition k at Temperature T0. Elements comes 
! from Subroutine "Readline".	       
!
! Variables:
!
! Input/Output Parameters of Routine (Arguments or Common)
! ----------------------------------
!	k       : marks the transition in study (Input).
!	E       : lower state energy (cm-1).
!	SWeig00 : State Statistical Weight (also called g00 in HITRAN). UNITLESS.
! molN    : = Na·Nmcon; molecular number density (molecules·cm-3)
! PFmol   : internal partition function. UNITLESS.
! ptype   : this parameter marks the result Population calculation 
!           in other words:
!           ('hit') Acording to HITRAN96 
!           ('tra') Acording to Hartmann et. al 2006 
!
! Other important Input Quantities (through "module_common_var")
! --------------------------------
!	c2      : second radiation constant = h·c/kb = 1.4388 cm·K
!	T0	    : Reference Temperature in Kelvin, 296 (K, Input).
!
! Output Quantities (through Common Statements)
! -----------------
!	PopuL0  : Populations of the Lower state of the Lines
!	          at Temperature T0                 (Output)
!           NOTE: 
!             1) with Ha Tran Code, the Population is UNITLESS.
!             2) HITRAN definition:
!                PopuT0(k) = molN*g00(k)*exp(-c2*E(k)/T0)/(PFmol(iso,T0))
!                give the units of the molecular number density (molN)
!
! Accessed Files:  None
! --------------
!
! Called Routines: 'none'
! ---------------  
!
! Called By: 'DipCAL' (Dipole of the lower state CALculation)
! ---------
!
!	Double Precision Version
!
!	T.Mendaza, last change 25 Jan 2017
!--------------------------------------------------------------------------------------------------------------------
!
    use module_common_var
    use module_molecSp
    Implicit None
    integer*8   , intent(in)    :: nLines
    type (dta_MOL), intent(in)    :: molP
    type (dta_SDF), intent(inout) :: dta1
    integer*8                   :: iso
    integer*8                   :: i, j, k
    Double Precision              :: cte1,cte2, Erot_red, T
    Double Precision              :: molN, PFmol, PFr
!-----------------------------------------
    T = molP % Temp
!-----------------------------------------
!
! Check that Arrays for Results are Large Enough, Initialize
    If ( nLines .GT. nLmx) Then
        call sizeError("1000",nLines,nLmx)
        Stop
    End If
!
! Partition function ratio:
!
    !old 
    !PFr=PFmol(molP,T0)/PFmol(molP,T)
    ! QT0 = PFmol(molP,T0)
    ! QT  = PFmol(molP,T )
    ! Now the partition function came as an input from "arts_interface"
    PFr=molP%QT0 / molP%QT
!
    DO k=1,nLines
        !Erot_red = dta1%J(k,1)*( dta1%J(k,1) + 1. )
        !cte1 = -c2*molP%B0*Erot_red/T0
        cte1 = -c2*dta1%E(k)/T0
        cte2 = -c2*dta1%E(k)*(1.d0/T-1.d0/T0)
        !
        dta1%PopuT0(k) = dta1%swei00(k)*exp(cte1)/molP%QT0
        dta1%PopuT(k)  = dta1%PopuT0(k)*PFr*dexp(cte2)
    ENDDO
!
    IF (ptype .eq. 'tra') then
        ! the Population remains as it has been calculated.
    ELSEIF (ptype .eq. 'hit') then
        ! Here we follow the HITRAN 1996 definition of 
        ! Population of the lower level using the molecular number density:
        ! molN = Na·Nmcon
        ! CH4 example (Nmcon = Nch4):
        ! Na   = 6.022·10^23 molecules·mol-1; Avogadro constant
        ! Nch4 = 41.245 mol/m3 = 41.245 mol/(1E006 *cm3) 
        !      = 41.245*1E-06  mol/cm3 (value used in common_var_module)
        ! NIST website:
        ! http://webbook.nist.gov/
        !
        molN= Na*molP%Nmcon ! molecule/cm3
        !molN= Nmcon    ! mol/cm3
        !
        DO k=1,nLines
            dta1%PopuT0(k) = molN*dta1%PopuT0(k)
            dta1%PopuT(k)  = molN*dta1%PopuT(k)
        ENDDO
    ELSE
        print*, "PopuCAL:Non-Valid Population calculation type:", ptype
    ENDIF
!
    Return
  END SUBROUTINE PopuCAL
!--------------------------------------------------------------------------------------------------------------------
  double precision function K_jkCalc(j,k,h,nLines,molP,PerM)
!--------------------------------------------------------------------------------------------------------------------
! "K_jkCalc": OFF-diagonal matrix elements (Wk'k ∞ Q_k'k)
!
! Detailed Description:
! ---------------------
! this function gives the theoretic formulation of relaxation matrix elements:
! 			<< k' | W | k >> 
! They have the same functional form as the rotational state-to-state 
! cross sections within a single vibrational state.
! NOTE: this depends on Temperature but through the use of "Module_common_var"
!
! Accessed Files:  None
! --------------
!
! Called Routines: 
! ---------------  
!
! Called Functions: 'wigner3j' ; 'wigner6j' ; 'AFmol_X' ; 'Ql_mol_X'
! ---------------- 
!
! Called By: 'WelCAL' 
! ---------
!
! Double Precision Version
!
! T.Mendaza, last change 26 January 2017
!--------------------------------------------------------------------------------------------------------------------
!
      use module_common_var
      use module_molecSp
      use module_maths
      IMPLICIT none
      integer*8   , intent(in):: j,k,nLines
      type (dta_SDF), intent(in):: h
      type (dta_MOL), intent(in):: molP, PerM
      !internal variables:
      integer*8       :: L, incr, step
      integer*8       :: i, iniL,endL
      double precision:: Ji, Ji_p, Jf, Jf_p !, jmax
      !double precision :: Ja(nLmx)
      integer*8       :: pos_L_2, aux_j_L, pos_Ji_2
      integer*8       :: li,lf
      double precision      :: suma, cte1, cte2
      double precision      :: w3j1,w3j2,w6j,w3jaux
      double precision      :: AF1,AF2, Qaux
      double precision      :: Afmol_X, Ql_mol_X 
!
!...co2  -- 
! This expresion is used for "downward" transitions (k"->k') only.
! That means that, given two transitions: 
! k := j_f , l_f , j_i , l_i
! k':= j'_f, l'_f, j'_i, l'_i
! so "downward" refears to that ones which (j'_i < j_i)
! the upwards ones can be obtained from detailed balance: 
!           rho_k'·<<k"|W|k'>> = rho_k"·<<k'|W|k">>
! NOTE: for further info about this formulation [Niro et al. 2004]
!
! STEP = 2 -> Because the quadrupole-quadrupole energy potential is dominant
!      in linear molecules.
  !
    !
      step = 2 
      ! j -> Ji(j) > Ji(k)=Ji'
      Ji   = h%J(j,1)
      Jf   = h%J(j,2)
      ! k-> 
      Ji_p = h%J(k,1)
      Jf_p = h%J(k,2)
      !
      ! 
      iniL=max(abs(Ji-Ji_p),abs(Jf-Jf_p))
      If( mod(iniL,step) .ne. 0 )iniL=iniL+1
      endL=min((Ji+Ji_p),(Jf+Jf_p))
      !
!
        ! li and lf are the vibrational angular momentum of the vibrational
        ! available for Linear molecules in HITRAN; 
        ! (due to its symmetry though).
!
       li = h%lv2(1)
       lf = h%lv2(2) 
!
      cte1 = (2.*Ji_p + 1.)*dsqrt((2.*Jf + 1.)*(2.*Jf_p + 1.))
      !cte2 = (-1)**(li+lf)
      ! A faster computation would be... (case (li+lf) has integer value)
      if ( mod( (li+lf+K_t+1),2 ) .eq. 0) then
        !(li+lf) is even
        cte2 = 1
      else !if (modulo( (li+lf),2) .eq. 1) then! 
        !(li+lf) is odd  
        cte2 = -1
      endif  
      suma = 0.d0
      !
      !do L = 1,130
      do L = iniL,endL,step
        ! Since the molecules to be considered are linear -> simmetryc 
        ! and the quadrupole - quadrupole interaction for the sistems mol-X
        ! are meant to be dominant the interaction potential is likely 
        ! dominated by EVEN rank contributions. 
        !
        ! Calculation of the basis rates "Q_mol-X":
        !if (mod(L,2) .eq. 0) then 
        ! L is even
          !print*, "Q ch4-X:"
          Qaux = Ql_mol_X(molP,L)
          !print*,"Qch4-x=", Qaux
        !else 
        ! L is odd
        ! because of that, the basis rates are expected to be small -> 0.
        !  Qaux = 0.0_dp
        !endif
        if ( Qaux .ne. 0.0_dp ) then
          !
          ! wigner 3j-Symbol 1:
          ! ( Ji  L  Ji' )
          ! ( li  0  -li )
          !http://www-stone.ch.cam.ac.uk/wigner.shtml
          w3j1 = wigner3j( Ji_p, real(L, dp), Ji , real(li, dp), 0.0_dp, real(-li, dp) )
          ! CASE wig3j0(M,J1,J2,J)
          !                       m (j1 j2 j)  
          !                    (-)  (m  -m 0)  
          !
          ! wigner 3j-Symbol 2:
          ! ( Jf  L  Jf' )
          ! ( lf  0  -lf )
          w3j2  = wigner3j( Jf_p, real(L, dp), Jf, real(-lf, dp), 0.0_dp, real(lf, dp) )
          !
          ! wigner 6j-Symbol:
          ! { Ji   Jf   1 }
          ! { Jf'  Ji'  L }
          w6j   = wigner6j(Ji,Jf, real(K_t,dp), Jf_p, Ji_p, real(L, dp))
          !
          ! Adiabatic factor 2:
          ! it depends on the temperature too but it is included through common module.
          AF2   = AFmol_X(molP, PerM, real(L,dp), step)
          !
          if ( .not.( isnan(w3j1  ) .or. isinf(w3j1  )) .and. &
               .not.( isnan(w3j2  ) .or. isinf(w3j2  )) .and. &
               .not.( isnan(w6j   ) .or. isinf(w6j   )) .and. &
               .not.( isnan(AF2) .or. isinf(AF2)) ) then
              suma = suma + w3j1 * w3j2* w6j * (2.*L + 1.) * &
                     Qaux / AF2
          else
            call offdiagEL_IN_ERROR(AF2, Qaux, w3j1, w3j2, w6j, 0.0_dp, 0.0_dp)
          endif
        else 
        endif
      enddo
      ! 
      ! Adiabatic factor 1:
      AF1= AFmol_X(molP, PerM, Ji, step)
      K_jkCalc = cte1*cte2*AF1*suma
      !K_jkCalc = cte1*cte2*suma
      if ( isnan( K_jkCalc ) .or. isinf( K_jkCalc  ) ) then
        call offdiagEL_ERROR(cte1, cte2, AF1, suma)
      endif
    RETURN
  END function K_jkCalc
!--------------------------------------------------------------------------------------------------------------------
  double precision function K_jkO2(j,k,h,nLines,molP,PerM)
!--------------------------------------------------------------------------------------------------------------------
! "K_jkO2": OFF-diagonal matrix elements for O2 systems(Wk'k ∞ Q_k'k)
!
! Detailed Description:
! ---------------------
! Same as K_jkcalc
!
! Accessed Files:  None
! --------------
!
! Called Routines: 
! ---------------  
!
! Called Functions: 'wigner3j' ; 'wigner6j' ; 'AFmol_X' ; 'Ql_mol_X'
! ---------------- 
!
! Called By: 'WelCAL' 
! ---------
!
! Double Precision Version
!
! T.Mendaza, last change 26 January 2017
!--------------------------------------------------------------------------------------------------------------------
!
      use module_common_var
      use module_molecSp
      use module_maths
      IMPLICIT none
      integer*8   , intent(in):: j,k,nLines
      type (dta_SDF), intent(in):: h
      type (dta_MOL), intent(in):: molP, PerM
      !internal variables:
      integer*8             :: L, incr, step
      integer*8             :: i, iniL,endL
      double precision      :: Ji, Ji_p, Jf, Jf_p !, jmax
      integer*8             :: Ni, Ni_p, Nf, Nf_p !, jmax
      real*8                :: Si, Sf,n
      integer*8             :: pos_L_2, aux_j_L, pos_Ji_2
      integer*8             :: li,lf
      real*8                :: suma, cte1, cte2
      double precision      :: w3j1,w3j2
      double precision      :: w6j1, w6j2, w6j3
      double precision      :: AF1,AF2, Qaux
      double precision      :: Afmol_X, Ql_mol_X 
!
!...O2  -- 
! This expresion is used for "downward" transitions (k"->k') only.
! That means that, given two transitions: 
! k := jf , Nf , ji , Ni, Si
! k':= jf', Nf', ji', Ni'
! so "downward" refears to that ones which (j'_i < j_i)
! the upwards ones can be obtained from detailed balance: 
!           rho_k'·<<k"|W|k'>> = rho_k"·<<k'|W|k">>
! NOTE: for further info about this formulation [Niro et al. 2004]
!
! STEP = 2 -> Because the quadrupole-quadrupole energy potential is dominant
!      in linear molecules.
  !
    !
      step = 2 
      ! j -> Ji(j) > Ji(k)=Ji'
      Ji   = h%J(j,1); Ni = h%N(j,1); Si = h%espin(j,1)
      Jf   = h%J(j,2); Nf = h%N(j,2); Sf = h%espin(j,2)
      ! k-> 
      Ji_p = h%J(k,1); Ni_p = h%N(k,1)
      Jf_p = h%J(k,2); Nf_p = h%N(k,2)
      !
      ! 
      iniL=max(abs(Ni-Ni_p),abs(Nf-Nf_p))
      If( mod(iniL,step) .ne. 0 )iniL=iniL+1
      endL=min((Ni+Ni_p),(Nf+Nf_p))
      !
      !
      cte1 = (2.*Ni + 1.)*(2.*Ni_p + 1.)*(2.*Nf + 1.)*(2.*Nf_p + 1.)
      cte1 = cte1*(2.*Jf + 1.)*(2.*Jf_p + 1.)*(2.*Ji_p + 1.)**2
      !
      if ( mod( int(Ji_p+Ji+K_t),2 ) .eq. 0) then
        cte2 = 1
      else 
        cte2 = -1
      endif  
      suma = 0.d0
      !
      !do L = 1,130
      do L = iniL,endL,step
        ! Since the molecules to be considered are linear -> simmetryc 
        ! and the quadrupole - quadrupole interaction for the sistems mol-X
        ! are meant to be dominant the interaction potential is likely 
        ! dominated by EVEN rank contributions. 
        !
        ! Calculation of the basis rates "Q_mol-X":
          Qaux = Ql_mol_X(molP,L)
        !
        if ( Qaux .ne. 0.0_dp ) then
          !
          ! wigner 3j-Symbol
          w3j1 = wigner3j( real(Ni_p, dp), real(Ni, dp), real(L, dp), 0.0_dp, 0.0_dp, 0.0_dp )
          w3j2 = wigner3j( real(Nf_p, dp), real(Nf, dp), real(L, dp), 0.0_dp, 0.0_dp, 0.0_dp )
          !
          ! wigner 6j-Symbol
          w6j1 = wigner6j(real(L, dp),Ji, real(Ni_p, dp), real(Si,dp), Ji_p, real(Ni, dp))
          w6j2 = wigner6j(real(L, dp),Jf, Jf_p, real(Sf,dp), real(Nf_p, dp), real(Nf, dp))
          w6j3 = wigner6j(real(L, dp),Ji, Ji_p, real(K_t,dp), Jf_p, Jf)
          !
          ! Adiabatic factor 2:
          ! it depends on the temperature too but it is included through common module.
          AF2  = AFmol_X(molP, PerM, real(L,dp), step)
          !
          if ( .not.( isnan(w3j1  ) .or. isinf(w3j1  )) .and. &
               .not.( isnan(w3j2  ) .or. isinf(w3j2  )) .and. &
               .not.( isnan(w6j1  ) .or. isinf(w6j1  )) .and. &
               .not.( isnan(w6j2  ) .or. isinf(w6j2  )) .and. &
               .not.( isnan(w6j3  ) .or. isinf(w6j3  )) .and. &
               .not.( isnan(AF2) .or. isinf(AF2)) ) then
              suma = suma + w3j1 * w3j2* w6j1 * w6j2 * w6j3 * &
                     (2.*L + 1.) * Qaux / AF2
          else
            call offdiagEL_IN_ERROR(AF2, Qaux, w3j1, w3j2, w6j1, w6j2, w6j3)
          endif
        else 
        endif
      enddo
      ! 
      ! Adiabatic factor 1:
      AF1= AFmol_X(molP, PerM, real(Ni,dp), step)
      K_jkO2 = cte1*cte2*AF1*suma
      if ( isnan( K_jkO2 ) .or. isinf( K_jkO2  ) ) then
        call offdiagEL_ERROR(cte1, cte2, AF1, suma)
      endif
    RETURN
  END function K_jkO2
!--------------------------------------------------------------------------------------------------------------------
  double precision function Ql_mol_X(molP,L)
!--------------------------------------------------------------------------------------------------------------------
! " Ql_mol_X": OFF-diagonal matrix elements (Wk'k ∞ Q_k'k)
! 
! Detailed description:
! ---------------------
! this function gives the theoretic formulation of relaxation matrix elements:
! 			<< k' | W | k >> 
! They have the same functional form as the rotational state-to-state 
! cross sections within a single vibrational state.
! Formulation developed by T.Mendaza developed from the work of [Niro et al. 2004].
!
      use module_common_var
      use module_molecSp
      use module_maths
      IMPLICIT none
      ! a1, a2, a3, dc were declared in module_common_var:
      integer*8        :: L
      real*8           :: E_l, T
      type(dta_MOL)    :: molP
!
! This expresion is used for "downward" transitions (k->k') only.
! That means that, given two transitions: 
! k := j_f l_f j_i l_i
! k':= j'_f l'_f j'_i l'_i
! so "downward" refears to that ones which (j'_i < j_i)
! the upwards ones can be obtained from detailed balance: rho_j·<<k|W|k'>> = rho_k·<<k'|W|k>>
! NOTE: for further info about this notation [Niro et al. 2004]
!
  !
    T = molP % Temp
    !reduce rotational energy from the level L * 1/B0
    E_l = real((L*L + L), dp) 
    !
    if (molP%M .eq. 7) then
      ! mol == O2
      if (L .eq. 0) then
        Ql_mol_X=0.0_dp
      else
      !Ql_mol_X = a1*( ( E_l )**(-a2) )*( exp(-a3*c*hplank*molP%B0*E_l/(T*kb))  )
        Ql_mol_X = molP%a1* &
               ( (  E_l )**(-molP%a2) )* &
               ( ( T/T0 )**(-molP%a3) )
               ![Tran et al., 2006]
      endif
    else 
      if (L .eq. 0) then
        Ql_mol_X=0.0_dp
      else
      !Ql_mol_X = a1*( ( E_l )**(-a2) )*( exp(-a3*c*hplank*molP%B0*E_l/(T*kb))  )
        Ql_mol_X = molP%a1* &
               ( ( E_l )**(-molP%a2) )* &
               ( exp(-molP%a3*c2*molP%B0*E_l/T)  )
               ![Niro et al., 2004]
      endif
    endif
    !write(*,*) 'Ql_mol_X=',Ql_mol_X
    RETURN
  END function Ql_mol_X
!--------------------------------------------------------------------------------------------------------------------
  double precision function AFmol_X(molP, PerM, L, step)
!--------------------------------------------------------------------------------------------------------------------
! "AFmol_X": Adiabaticity-Factor for the chemical reaction CH4-X 
!
! Description:
! The original expresion for this correction factor (a.k.a. Adiabaticity-Factor) 
! is given by Rodriguez et al. 1999.
! AF(J,T) = { 1 + 1/24 * [ (w_j_j2 * dc)/(v_mol-X) ]^2  }^-2
! UNITLESS
!
! where:
! w_j_jstep := frequency spacing 2*pi*c* |E(J) - E(J-step)|; E = Ev + Erot(J)
!           where the Vibrational energy (Ev) is the same for both J, and J-step
!           but not for the rotational part (Erot).
!           This leads to the reduced rotational energy Erot = B0*( J(J+1) )
!           |Erot(J) - Erot(J-step)| = step*B0[2J+1 - step]
! dc : scaling length (fitted parameter)
! v_mol-X := relative mean speed of the system. This value is obtained from the
!             Maxwell–Boltzmann distribution "most expected value":
!             <v> = ( (8*kb*T*mu)/pi)^1/2
!             where: kb is the Boltzmann cte. and mu is the reduced mass
!                    given by: mu = (1/m1 + 1/m2)
!
      use module_common_var
      IMPLICIT none
      integer*8, intent(in)      :: step
      double precision,intent(in):: L
      type(dta_MOL),intent(in)   :: molP, PerM
      double precision,Parameter :: cAF = 0.0006983_dp ! Adiabaticy Factor cte
                                 ! General constant non-molecule dependent included in the 
                                 ! Adiabaticy Factor. It is built as follows:
                                 ! * Velocity term 
                                 !   cte1 = 1./(Na·8·kb / (1E-03·pi))
                                 ! * Freq. spacing term 
                                 !   cte2 = (2·pi·c)^2 
                                 ! * AF correction form self-cte: 
                                 !   1/24
                                 ! * Unit conversion for dc (Å) to m 
                                 !                       dc = 2.2 Å = 2.2E-10 m
                                 !   (to cancel velocity units m/s):
                                 !   1 Å = 1E-10 m 
                                 !   NOTE: dc is squared in the formula
                                 !   cAF = (1./24.)*cte1*cte2*1E-20
                                 ! 
      double precision           :: v_mol_X, w_j_jstep
      double precision           :: mu, T
      !-----------------------------------------
            T = molP % Temp
      !-----------------------------------------
      !
      !
            if (L .eq. 0) then
              AFmol_X = 1.0_dp
            else
              mu = ( 1./molP%mms + 1./PerM%mms ) ! 1/(g/mol)
            ! Velocity Term:
            ! Since it is a velocity its units are m/s. 
            ! Example: N2 at T=300K 
            ! N2_M = 28.006148 g/mol -> 1E-03 · 28.006148 / Na Kg
            ! mu = 1/(1E-03·N2_M/Na) Kg^-1
            ! kb is in J/K where J = Kg·m2/s2
            ! vp =sqrt(2·kb·T·mu) ~ 420 m/s = 420·1E02 cm/s
            ! vm = sqrt(8·kb·T·mu/pi) ~ 473 m/s = 473 · 1E02 cm/s
            ! cte1 = dsqrt(Na*8*kb / (1E-03*pi)) (J/(K·mol))^1/2 = (Kg·m2/(K·s2·mol))^1/2
            ! vm = cte1·sqrt(T/N2_M) ~ 473 m/s
            ! 
              v_mol_X  = 1./T*mu ! ----> squared and inverted velocity term 
            !print*, "v=", v_mol_X
            !
            ! Frequency spacing:
              w_j_jstep  = (molP%B0*( L + L + 1 - step )*step)**2 ! ---> squared freq. spacing
            !print*, "w_j_j-2=", w_j_jstep
            ! Since we have to give this quantity in (s-1) then
            ! we need to multiply the: 
            ! |Erot(J) - Erot(J-step)| = step*B0[2J+1 - step]
            ! by cte2 = 2*pi*c <- for the units shake
            ! Regarding that B0 is in cm-1 then c must be in cm/s (which is the cte by default in this program)
            ! NOTE: careful if you change the units of the constants given in "module_common_var"
            !       lots of calculations in this program are fitted to the units given.
            ! w_j_jstep  = cte2*(molP%B0*( L + L + 1 - step )*step)**2 
            ! UNITS: s-1
            !
              AFmol_X = 1./( 1. + cAF*( v_mol_X )*( w_j_jstep ) * (molP%dc * molP%dc) )**2
            !
            endif
      RETURN
  END function AFmol_X
!--------------------------------------------------------------------------------------------------------------------
  subroutine sumRule(nLines,indexS,dipole,Wmat,dfact)
!--------------------------------------------------------------------------------------------------------------------
        use module_error
        implicit none
        integer (kind=8), intent(in)    :: nLines
        integer (kind=8), intent(in)    :: indexS(nLines)
        real*8, intent(in)              :: dfact
        Double Precision, intent(in)    :: dipole(nLines)
        Double Precision, intent(in)    :: Wmat(nLines,nLines)
        !local variables
        integer (kind=8)                :: i,j
        Double Precision                :: Saux, DipAux, Wij, Wii
        Double Precision, Parameter     :: TOL = 1.0E-02!
        logical                         :: isinf, isnan, testOK
!--------------------------------------------------------------------
        testOK = .true.
        DO i = 1, nLines
          Saux = 0.0D0
          DO j = 1, nLines
            if (j .ne. i) then
              DipAux = dipole( indexS(j) )/dipole( indexS(i) )
              Wij = Wmat(i,j)
              if ( isnan(   Wij) .or. isinf(   Wij) ) stop!Wij = 0.0D0
              if ( isnan(DipAux) .or. isinf(DipAux) ) stop!DipAux = 0.0D0
              Saux = Saux + DipAux*Wij
            else
              Wii = Wmat(i,i)*dfact ! dfact = diagonal factor
            endif
          ENDDO
          if ( (abs(Wii+Saux) .gt. TOL) .and. (i .ne. nLines) ) then
            write(*,'(a4,i3,a14)'), "row#",i,"NO SUM-RULE"
            print*, "half-Width", Wii , "?= SUMij{Wmat}", -Saux
            testOK = .false.
            !stop
          else
            !write(*,'(a4,i3,a14)'), "row#",i,"NO SUM-RULE"
            !print*, "half-Width", Wii , "?= SUMij{Wmat}", -Saux
          endif
        ENDDO
      
        print*, "SUM-RULE TEST FINISHED"
        if ( .not.(testOK) ) then
          call SumRuleError()
        else
          print*, "sumRule: The calculation correctly verifies the sum rule!"
        endif

  end subroutine sumRule
!--------------------------------------------------------------------------------------------------------------------
  SUBROUTINE VarInit(dta1,dta2,molP)
! SUBROUTINE VarInit(dta1,dta2,dta3,dta4,molP)
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
! dta1    : dta type "dta_SDF". Spectrocopic parameters.
! dta2    : dta type "dta_RMF". Relaxation Matrix parameters.
! molP    : dta type "dta_MOL". Molecular structure parameters.
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
! T.Mendaza, last change 11 February 2016
!--------------------------------------------------------------------------------------------------------------------
!
    use module_common_var
    implicit none
    type (dta_SDF), intent(inout)    :: dta1
    type (dta_RMF), intent(inout)    :: dta2
    type (dta_MOL), intent(inout)    :: molP
    integer*8                  :: i,j,k

!----------
print*, "Init. dta_SDF..."
!
! Init. SDF data
      dta1%M   = 0
      dta1%iso = 0 
      ! Local Q.
      !dta1%nu(6,1 or 2) ! lower = 1 or 'low'
      !dta1%nu = reshape( (/ 0.0, 0.0,&
      !                      0.0, 0.0,&
      !                      0.0, 0.0,&
      !                      0.0, 0.0,&
      !                      0.0, 0.0,&
      !                      0.0, 0.0/),(/6,2/),order=(/2,1/))
      
      !
      ! which is the same as 
      !
      !dta1%nu(1,1) = 0.0_dp ; dta1%nu(1,2) = 0.0_dp     
      !dta1%nu(2,1) = 0.0_dp ; dta1%nu(2,2) = 0.0_dp     
      !dta1%nu(3,1) = 0.0_dp ; dta1%nu(3,2) = 0.0_dp        
      !dta1%nu(4,1) = 0.0_dp ; dta1%nu(4,2) = 0.0_dp        
      !dta1%nu(5,1) = 0.0_dp ; dta1%nu(5,2) = 0.0_dp 
      !dta1%nu(6,1) = 0.0_dp ; dta1%nu(6,2) = 0.0_dp  
      !
      dta1%lv2 = (/0,0/)

      do  i=1,nLmx
          !Integer kind
           dta1%J(i,1)    = 0 ; dta1%J(i,2)   = 0 
          !Double precision
           dta1%Sig(i)    = 0.0_dp
           dta1%Str(i)    = 0.0_dp
           dta1%E(i)      = 0.0_dp
           dta1%HWT0(i)   = 0.0_dp
           dta1%BHW(i)    = 0.0_dp
           dta1%SHIFT(i)  = 0.0_dp
           dta1%swei0(i)  = 0.0_dp
           dta1%swei00(i) = 0.0_dp
           dta1%PopuT0(i) = 0.0_dp
           dta1%PopuT(i)  = 0.0_dp
           dta1%D0(i)     = 0.0_dp 
           dta1%Drigrotor(i) = 0.0_dp 
          ! Global Q.
          !integer
          dta1%N(i,1)    = 0 ; dta1%N(i,2)    = 0  
          !real
          dta1%nspin(i,1)= 0.0; dta1%nspin(i,2) = 0.0
          dta1%espin(i,1)= 0.0; dta1%espin(i,2) = 0.0
          !dp
          dta1%F(i,1) = 0.0_dp ; dta1%F(i,2) = 0.0_dp
!         dta1%J_g6(i)= 0.0_dp, dta1%Ju_g6(i)= 0.0_dp
          !character
           dta1%br(i)   = "" ; dta1%br_N(i) = ""
      enddo
!
print*, "Init. dta_RMF..."
!
! Init. RMF data
      do  i=1,nMmx
          !Double precision
           dta2%WT0(i)    = 0.0_dp
           dta2%BTW(i)    = 0.0_dp
      enddo

print*, "Init. Molecular data type..."
!
! Init. MOLECULE data
      molP%M = 0
      do  i=1,mMax
          !Integer kind
           molP%Aco = 0 
          !Double precision
           molP%mms    = 0.0_dp
      enddo
! 
      Return
  END SUBROUTINE VarInit
!--------------------------------------------------------------------------------------------------------------------